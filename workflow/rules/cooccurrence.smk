# ══════════════════════════════════════════════════════════════════════════════
# Co-occurrence ARG × MGE — jointure contig-level
# RGI main (assembly) × mobileOG-db
# ══════════════════════════════════════════════════════════════════════════════

rule cooccurrence_arg_mge:
    """
    Détection des co-occurrences ARG × MGE au niveau contig.
    Jointure : RGI main (Contig) × mobileOG (Specific Contig).
    Outputs : paires, résumé drug_class × MGE, edge list réseau, stats/sample.
    """
    input:
        rgi_files = expand("results/{sample}/rgi/{sample}.txt",
                           sample=config["samples_id"]),
        mge_files = expand("results/{sample}/mge/assembly.fasta.mobileOG.Alignment.Out.csv",
                           sample=config["samples_id"]),
    output:
        pairs   = "results/r_analysis/cooccurrence/01_cooccurrence_pairs.tsv",
        summary = "results/r_analysis/cooccurrence/02_cooccurrence_summary.tsv",
        network = "results/r_analysis/cooccurrence/03_cooccurrence_network.tsv",
        stats   = "results/r_analysis/cooccurrence/04_cooccurrence_stats.tsv",
    conda:
        "../envs/r_analysis.yaml"
    log:
        "logs/r_analysis/cooccurrence_arg_mge.log"
    script:
        "../scripts/cooccurrence_arg_mge.R"


rule cooccurrence_report:
    """Rapport HTML co-occurrence ARG × MGE."""
    input:
        pairs   = "results/r_analysis/cooccurrence/01_cooccurrence_pairs.tsv",
        network = "results/r_analysis/cooccurrence/03_cooccurrence_network.tsv",
    output:
        report  = "results/figures/cooccurrence/cooccurrence_report.html",
    params:
        dpi = config.get("fig_dpi", 150),
    conda:
        "../envs/python_viz.yaml"
    log:
        "logs/figures/cooccurrence_report.log"
    script:
        "../scripts/plot_cooccurrence_report.py"


rule cooccurrence_all:
    """Target rule — co-occurrence ARG × MGE complète."""
    input:
        rules.cooccurrence_arg_mge.output,
        rules.cooccurrence_report.output,


