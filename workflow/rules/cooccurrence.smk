# ══════════════════════════════════════════════════════════════════════════════
# ARG × MGE co-occurrence — contig-level join
# RGI main (assembly) × mobileOG-db
# ══════════════════════════════════════════════════════════════════════════════

rule cooccurrence_arg_mge:
    """
    Detect ARG × MGE co-occurrence at the contig level.
    Join: RGI main (Contig) × mobileOG (Specific Contig).
    Outputs: pairs, drug_class × MGE summary, network edge list, per-sample stats.
    """
    input:
        rgi_files = expand("results/{sample}/rgi/{sample}.txt",
                           sample=config["samples_id"]),
        mge_files = expand("results/{sample}/mge/{sample}.assembly.fasta.mobileOG.Alignment.Out.csv",
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
    """ARG × MGE co-occurrence HTML report."""
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
    """Target rule — complete ARG × MGE co-occurrence module."""
    input:
        rules.cooccurrence_arg_mge.output,
        rules.cooccurrence_report.output,
