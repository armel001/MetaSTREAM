rule generate_arg_visualizations:
    input:
        arg_matrix_normalized = "results/r_analysis/{tool}/arg_matrix_normalized.tsv",
        arg_matrix_counts     = "results/r_analysis/{tool}/arg_matrix_counts.tsv",
        arg_matrix_presence   = "results/r_analysis/{tool}/arg_matrix_presence.tsv",
        drug_class            = "results/r_analysis/{tool}/drug_class_abundance.tsv",
        mechanism             = "results/r_analysis/{tool}/mechanism_abundance.tsv",
        family                = "results/r_analysis/{tool}/family_abundance.tsv",
        arg_counts            = "results/r_analysis/{tool}/arg_counts.tsv",
        family_prefix         = "results/r_analysis/{tool}/family_prefix_abundance.tsv",
        arg_normalized = "results/r_analysis/{tool}/arg_normalized.tsv",
        config_file           = "config/config.yaml"
    output:
        report = "results/figures/rgi/{tool}_report.html"
    params:
        tool         = lambda w: w.tool,
        min_identity = lambda w: config[w.tool]["min_identity"],
        min_coverage = lambda w: config[w.tool]["min_coverage"]
    log:
        "logs/python_visualization/{tool}_report.log"
    conda:
        "../envs/python_viz.yaml"
    wildcard_constraints:
        tool = "rgi|rgi_bwt"
    threads: 1
    script:
        "../scripts/generate_arg_visualizations.py"
