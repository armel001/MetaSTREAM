rule reads_stats_raw:
    input:
        raw = "resources/reads/{sample}.fastq.gz"
    output:
        stats = "results/{sample}/stats/reads_stats_raw.txt"
    log:
        "logs/reads_stats/{sample}_raw.log"
    conda:
        "../envs/python_filter.yaml"
    threads: 1
    script:
        "../scripts/reads_stats.py"


rule reads_stats_noh:
    input:
        fastq = "results/{sample}/clean/{sample}_noh.fq.gz"
    output:
        stats = "results/{sample}/stats/reads_stats_noh.txt"
    log:
        "logs/reads_stats/{sample}_noh.log"
    conda:
        "../envs/python_filter.yaml"
    threads: 1
    script:
        "../scripts/reads_stats.py"


rule reads_stats_filtered:
    input:
        fastq = "results/{sample}/fastplong/fp_trimmed.fastq.gz"
    output:
        stats = "results/{sample}/stats/reads_stats_filtered.txt"
    log:
        "logs/reads_stats/{sample}_filtered.log"
    conda:
        "../envs/python_filter.yaml"
    threads: 1
    script:
        "../scripts/reads_stats.py"


rule aggregate_reads_stats:
    input:
        raw = expand("results/{sample}/stats/reads_stats_raw.txt",
                    sample=config["samples_id"]),
        noh = expand("results/{sample}/stats/reads_stats_noh.txt",
                    sample=config["samples_id"]),
        filtered = expand("results/{sample}/stats/reads_stats_filtered.txt",
                         sample=config["samples_id"])
    output:
        aggregated = "results/stats/reads_stats_all.tsv",
        summary = "results/stats/reads_stats_summary.txt"
    params:
        samples = config["samples_id"]
    log:
        "logs/reads_stats/aggregate.log"
    conda:
        "../envs/python_filter.yaml"
    threads: 1
    script:
        "../scripts/aggregate_reads_stats.py"
