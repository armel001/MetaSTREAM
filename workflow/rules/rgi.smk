rule rgi:
    input:
        assembly = "results/{sample}/medaka/{sample}_polished.fasta"
    output:
        txt = "results/{sample}/rgi/{sample}.txt",
        json = "results/{sample}/rgi/{sample}.json",
        report_dir = directory("results/{sample}/rgi")
    conda:
        "../envs/rgi.yaml"
    threads: 20
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.report_dir}

        echo "==> Running RGI against configured local DB (output in {output.report_dir})."

        rgi main \
            --input_sequence {input.assembly} \
            -o {output.report_dir}/{wildcards.sample} \
            --input_type contig \
            --num_threads {threads} \
            --local \
            --clean
        """

