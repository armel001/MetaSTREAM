rule kraken2:
    input:
        fastq = "results/{sample}/clean/{sample}_noh.fq.gz"
    output:
        report = "results/{sample}/kraken/{sample}_report.txt",
        hits = "results/{sample}/kraken/{sample}_output.txt"
    params:
        db = config["kraken_db"],
        confidence = config["kraken_confidence"]
    threads: config["kraken_threads"]
    resources:
        mem_mb = 100000 # Availablee 190 Go
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        kraken2 --db {params.db} \
                --threads {threads} \
                --confidence {params.confidence} \
                --report {output.report} \
                --output {output.hits} \
                {input.fastq}
        """
