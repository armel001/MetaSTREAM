rule flye_assembly:
    input:
        "results/{sample}/clean/{sample}_noh.fq.gz"
    output:
        "results/{sample}/assembly.fasta"
    conda:
        "../envs/flye.yaml"
    log:
        "logs/flye_assembly/{sample}.log"
    benchmark:
        "benchmarks/flye_assembly/{sample}.txt"
    threads: 22
    resources:
        mem_mb = 128000
    shell:
        """
        flye \
            --nano-hq {input} \
            --out-dir results/{wildcards.sample}/flye_tmp \
            --threads {threads} \
            --meta \
            2> {log}

        mv results/{wildcards.sample}/flye_tmp/assembly.fasta {output}

        rm -rf results/{wildcards.sample}/flye_tmp
        """
