rule medaka:
    input:
        assembly = "results/{sample}/assembly.fasta",
        reads    = "results/{sample}/clean/{sample}.renamed.fastq.gz"
    output:
        consensus = "results/{sample}/medaka/consensus.fasta",
        final     = "results/{sample}/medaka/{sample}_polished.fasta"
    params:
        outdir = "results/{sample}/medaka",
        model  = config["medaka_model"]  # ex: "r1041_e82_400bps_hac_v4.3.0"
    conda:
        "../envs/medaka.yaml"
    threads: 12
    log:
        "logs/medaka/{sample}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}

        echo "==> Polissage Medaka de {wildcards.sample} (modèle: {params.model})"

        medaka_consensus \
            -i {input.reads} \
            -d {input.assembly} \
            -o {params.outdir} \
            -t {threads} \
            -m {params.model} \
            2>> {log}

        cp {output.consensus} {output.final}
        """
