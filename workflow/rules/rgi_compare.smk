rule compare_rgi_methods:
    input:
        rgi     = "results/summary/rgi_all_samples.tsv",
        rgi_bwt = "results/summary/rgi_bwt_all_samples.tsv"
    output:
        tsv     = "results/summary/rgi_comparison.tsv",
        html    = "results/summary/rgi_comparison_report.html"
    log:
        "logs/rgi_compare/compare.log"
    conda:
        "../envs/rgi_preprocessing.yaml"
    threads: 1
    resources:
        mem_mb = 4000
    script:
        "../scripts/compare_rgi_methods.py"
