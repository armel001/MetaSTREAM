rule mge_detection:
    """Run mobileOG-pl pipeline on assembled contigs."""
    input:
        assembly = "results/{sample}/assembly.fasta",
        db       = config["mobileog"]["db"],
        metadata = config["mobileog"]["metadata"]
    output:
        summary  = "results/{sample}/mge/{sample}.assembly.fasta.summary.csv",
        hits     = "results/{sample}/mge/{sample}.assembly.fasta.mobileOG.Alignment.Out.csv",
        tsv      = "results/{sample}/mge/{sample}.assembly.fasta.tsv",
        proteins = "results/{sample}/mge/{sample}.assembly.fasta.faa"
    params:
        k       = config["mobileog"]["k"],
        evalue  = config["mobileog"]["evalue"],
        pident  = config["mobileog"]["min_identity"],
        qcov    = config["mobileog"]["min_coverage"],
        outdir  = "results/{sample}/mge",
        wrapper = "workflow/scripts/run_mge.sh",
        pypath  = config["mobileog"]["kyanite_py"]
    log:
        "logs/mge/{sample}_mge.log"
    conda:
        "../envs/mobileog.yaml"
    threads: config["resources"]["threads"]
    shell:
        """
        bash {params.wrapper} \
            -i {input.assembly} \
            -d {input.db} \
            -m {input.metadata} \
            -k {params.k} \
            -e {params.evalue} \
            -p {params.pident} \
            -q {params.qcov} \
            --outdir {params.outdir} \
            --pypath {params.pypath} \
            > {log} 2>&1

        mv {params.outdir}/assembly.fasta.faa \
           {params.outdir}/{wildcards.sample}.assembly.fasta.faa
        mv {params.outdir}/assembly.fasta.tsv \
           {params.outdir}/{wildcards.sample}.assembly.fasta.tsv
        mv {params.outdir}/assembly.fasta.mobileOG.Alignment.Out.csv \
           {params.outdir}/{wildcards.sample}.assembly.fasta.mobileOG.Alignment.Out.csv
        mv {params.outdir}/assembly.fasta.summary.csv \
           {params.outdir}/{wildcards.sample}.assembly.fasta.summary.csv
        """

rule mge_aggregate:
    """Aggregate MGE summary files across all samples."""
    input:
        summaries = expand(
            "results/{sample}/mge/{sample}.assembly.fasta.summary.csv",
            sample=config["samples_id"]
        )
    output:
        per_contig = "results/summary/mge_all_samples.tsv",
        per_sample = "results/summary/mge_per_sample.tsv"
    log:
        "logs/mge/aggregate.log"
    conda:
        "../envs/mobileog.yaml"
    threads: 8
    run:
        import pandas as pd

        log_f = open(log[0], 'w')

        MGE_CLASSES = [
            "Bacteriophages",
            "Insertion sequences",
            "Integrative elements",
            "Multiple",
            "Plasmids"
        ]

        dfs          = []
        summary_rows = []

        for f, sample in zip(input.summaries, config["samples_id"]):
            try:
                df = pd.read_csv(f, index_col=0)
                df.insert(0, 'sample_id', sample)
                dfs.append(df)

                row = {'sample_id': sample}
                for cls in MGE_CLASSES:
                    if cls in df.columns:
                        row[cls] = df[cls].sum()
                row['Total_hits'] = df['Total Number of Hits'].sum() \
                                    if 'Total Number of Hits' in df.columns else 0
                row['Total_contigs_with_MGE'] = len(df)
                row['Total_unique_ORFs'] = df['Amount of Unique ORFs'].sum() \
                                           if 'Amount of Unique ORFs' in df.columns else 0

                total = row['Total_hits']
                for cls in MGE_CLASSES:
                    pct_col = f"Percent_{cls.replace(' ', '_')}"
                    row[pct_col] = (row[cls] / total * 100) if total > 0 else 0

                summary_rows.append(row)
                print(f"  ✓ {sample}: {len(df)} contigs, "
                      f"{int(row['Total_hits'])} total hits", file=log_f)

            except Exception as e:
                print(f"  ⚠ {sample}: {e}", file=log_f)

        if dfs:
            all_contigs = pd.concat(dfs, ignore_index=True)
            all_contigs.to_csv(output.per_contig, sep='\t', index=False)

        if summary_rows:
            per_sample = pd.DataFrame(summary_rows)
            per_sample.to_csv(output.per_sample, sep='\t', index=False)

        log_f.close()


rule mge_r_analysis:
    """Matrices R MGE — hits filtrés, normalisation, Pident, présence/absence."""
    input:
        hits_files = expand(
            "results/{sample}/mge/{sample}.assembly.fasta.mobileOG.Alignment.Out.csv",
            sample=config["samples_id"]
        ),
        stats = "results/stats/sequencing_stats.tsv",
    output:
        hits_raw    = "results/r_analysis/mge/01_mge_hits_raw.tsv",
        normalized  = "results/r_analysis/mge/07_mge_normalized.tsv",
        norm_matrix = "results/r_analysis/mge/08_mge_normalized_matrix.tsv",
        pident      = "results/r_analysis/mge/12_mge_pident_summary.tsv",
        presence    = "results/r_analysis/mge/14_mge_presence_absence.tsv",
    params:
        outdir     = "results/r_analysis/mge",
        min_pident = config.get("mge_min_pident", 90),
        min_cov    = config.get("mge_min_cov",    90),
    conda:
        "../envs/r_analysis.yaml"
    log:
        "logs/r_analysis/mge_analysis.log"
    script:
        "../scripts/mge_analysis.R"


rule mge_report:
    """Generate MGE HTML report from R analysis outputs."""
    input:
        normalized = "results/r_analysis/mge/07_mge_normalized.tsv",
        hits_raw   = "results/r_analysis/mge/01_mge_hits_raw.tsv",
    output:
        report = "results/figures/mge/mge_report.html"
    params:
        dpi = config.get("fig_dpi", 150)
    conda:
        "../envs/python_viz.yaml"
    log:
        "logs/figures/mge_report.log"
    script:
        "../scripts/plot_mge_report.py"
