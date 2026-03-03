rule rgi_bwt:
    input:
        reads = "results/{sample}/clean/{sample}.renamed.fastq.gz"
    output:
        report_dir = directory("results/{sample}/rgi_bwt"),
        txt        = "results/{sample}/rgi_bwt/{sample}.txt",          # standardisé
        allele_txt = "results/{sample}/rgi_bwt/{sample}.allele_mapping_data.txt",
        gene_txt   = "results/{sample}/rgi_bwt/{sample}.gene_mapping_data.txt"
    conda:
        "../envs/rgi.yaml"
    threads: 20
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.report_dir}

        echo "==> Running RGI BWT on ONT reads (output in {output.report_dir})."

        rgi bwt \
            --read_one {input.reads} \
            --output_file {output.report_dir}/{wildcards.sample} \
            --threads {threads} \
            --local \
            --include_other_models \
            --include_wildcard

        # Normalisation : on expose gene_mapping_data comme .txt standard
        cp {output.report_dir}/{wildcards.sample}.gene_mapping_data.txt \
           {output.report_dir}/{wildcards.sample}.txt
        """
