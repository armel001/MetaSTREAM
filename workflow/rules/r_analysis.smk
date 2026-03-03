rule normalize_and_matrices:
    input:
        rgi_aggregated   = "results/summary/{tool}_all_samples.tsv",
        sequencing_stats = "results/stats/sequencing_stats.tsv"
    output:
        arg_counts            = "results/r_analysis/{tool}/arg_counts.tsv",
        arg_relative          = "results/r_analysis/{tool}/arg_relative.tsv",
        arg_normalized        = "results/r_analysis/{tool}/arg_normalized.tsv",
        arg_presence          = "results/r_analysis/{tool}/arg_presence.tsv",
        arg_matrix_counts     = "results/r_analysis/{tool}/arg_matrix_counts.tsv",
        arg_matrix_relative   = "results/r_analysis/{tool}/arg_matrix_relative.tsv",
        arg_matrix_normalized = "results/r_analysis/{tool}/arg_matrix_normalized.tsv",
        arg_matrix_presence   = "results/r_analysis/{tool}/arg_matrix_presence.tsv",
        drug_class_abundance  = "results/r_analysis/{tool}/drug_class_abundance.tsv",
        mechanism_abundance   = "results/r_analysis/{tool}/mechanism_abundance.tsv",
        family_abundance      = "results/r_analysis/{tool}/family_abundance.tsv",
        family_prefix_abundance = "results/r_analysis/{tool}/family_prefix_abundance.tsv"  # ← nouveau
    params:
        tool = lambda w: w.tool
    log:
        "logs/r_analysis/{tool}/normalize_and_matrices.log"
    conda:
        "../envs/r_analysis.yaml"
    wildcard_constraints:
        tool = "rgi|rgi_bwt"
    threads: 1
    script:
        "../scripts/normalize_and_matrices.R"
