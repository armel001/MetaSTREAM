rule host_depletion:
    input:
        reads   = "results/{sample}/fastplong/fp_trimmed.fastq.gz",
        host_ref = ancient(config["host_genome"])
    output:
        clean_reads = "results/{sample}/clean/{sample}_noh.fq.gz"
    log:
        "logs/host_depletion/{sample}.log"
    benchmark:
        "benchmarks/host_depletion/{sample}.txt"
    threads: 40
    resources:
        mem_mb  = 64000,    # 64 Gb RAM pour minimap2
        disk_mb = 50000     # 50 Gb espace disque temporaire
    conda:
        "../envs/host_depletion.yaml"
    shell:
        """
        echo "[$(date)] Starting host depletion for {wildcards.sample}" > {log}
        echo "[$(date)] Input reads: {input.reads}" >> {log}

        minimap2 \
            -ax map-ont \
            -t {threads} \
            {input.host_ref} \
            {input.reads} \
            2>> {log} | \
        samtools view \
            -b \
            -f 4 \
            -@ {threads} \
            2>> {log} | \
        samtools fastq \
            -@ {threads} \
            - \
            2>> {log} | \
        bgzip \
            -l 9 \
            -@ {threads} \
            > {output.clean_reads} \
            2>> {log}

        # Stats post-depletion
        echo "[$(date)] Counting reads after host depletion..." >> {log}
        echo "Reads after host depletion: $(zcat {output.clean_reads} | \
             awk 'NR%4==1' | wc -l)" >> {log}

        echo "[$(date)] Host depletion completed for {wildcards.sample}" >> {log}
        """
