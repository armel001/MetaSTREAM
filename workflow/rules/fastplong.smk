rule fastplong:
    input:
        "/data/armel/mythesis-data/data/{sample}.fastq.gz"
    output:
        fq   = "results/{sample}/fastplong/fp_trimmed.fastq.gz",
        html = "results/{sample}/fastplong/fastplong.html",
        json = "results/{sample}/fastplong/fastplong.json"
    conda:
        "../envs/fastplong.yaml"
    threads: 16
    resources:
        mem_mb  = 16000,
        disk_mb = 60000
    log:
        "logs/fastplong/{sample}.log"
    benchmark:
        "benchmarks/fastplong/{sample}.txt"
    shell:
        r"""
        fastplong \
            -i {input[0]} \
            -o {output.fq} \
            -q 8 \
            -t {threads} \
            -h {output.html} \
            -j {output.json} \
            &> {log}
        """
