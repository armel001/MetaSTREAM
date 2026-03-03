# ─── helpers ────────────────────────────────────────────────────────────────

def get_tool_config(wildcards, param):
    return config[wildcards.tool][param]


# ─── rules ──────────────────────────────────────────────────────────────────

rule filter_rgi_quality:
    input:
        rgi = "results/{sample}/{tool}/{sample}.txt"
    output:
        filtered = "results/{sample}/{tool}_preprocessing/{sample}_filtered.txt"
    params:
        min_identity = lambda w: config[w.tool]["min_identity"],
        min_coverage  = lambda w: config[w.tool]["min_coverage"],
        model_types   = lambda w: config[w.tool]["model_types"]
    log:
        "logs/{tool}_filter/{sample}.log"
    conda:
        "../envs/rgi_preprocessing.yaml"
    threads: 1
    resources:
        mem_mb = 2000
    script:
        "../scripts/filter_rgi_results.py"


rule select_rgi_columns:
    input:
        filtered = "results/{sample}/{tool}_preprocessing/{sample}_filtered.txt"
    output:
        selected = "results/{sample}/{tool}_preprocessing/{sample}_selected.txt"
    params:
        columns = lambda w: config[w.tool]["columns_to_keep"]
    log:
        "logs/{tool}_select_columns/{sample}.log"
    conda:
        "../envs/rgi_preprocessing.yaml"
    threads: 1
    script:
        "../scripts/select_rgi_columns.py"


rule aggregate_rgi_samples:
    input:
        selected = expand(
            "results/{sample}/{{tool}}_preprocessing/{sample}_selected.txt",
            sample = config["samples_id"]
        )
    output:
        aggregated = "results/summary/{tool}_all_samples.tsv",
        summary    = "results/summary/{tool}_aggr_summary.txt"
    params:
        samples = config["samples_id"],
        tool    = lambda w: w.tool
    log:
        "logs/{tool}_aggregate/aggregate.log"
    conda:
        "../envs/rgi_preprocessing.yaml"
    threads: 8
    resources:
        mem_mb = 80000
    wildcard_constraints:
        tool = "rgi|rgi_bwt"
    script:
        "../scripts/aggregate_rgi_samples.py"
