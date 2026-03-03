rule rename_reads:
    input:
        reads = "results/{sample}/clean/{sample}_noh.fq.gz"
    output:
        renamed = temp("results/{sample}/clean/{sample}.renamed.fastq.gz")
    log:
        "logs/rename_reads/{sample}.log"
    benchmark:
        "benchmarks/rename_reads/{sample}.txt"
    threads: 4
    conda:
        "../envs/python_filter.yaml"
    shell:
        """
        zcat {input.reads} \
            | awk 'NR%4==1 {{printf("@{wildcards.sample}_%d\\n", ++i)}} NR%4!=1 {{print}}' \
            | gzip -c \
            > {output.renamed} \
            2> {log}
        """

