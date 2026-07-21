# ─── Integrative and Conjugative Elements — ICEfinder2 ───────────────────────
# Pipeline: assemblage poli (Medaka) → ICEfinder2 mode Metagenome (Kraken2
#   pour la classification taxonomique par contig avant scan HMM/BLAST)
# Reference: Wang M. et al. 2024, Nucleic Acids Research.
#            "ICEberg 3.0: functional categorization and analysis of the
#            integrative and conjugative elements in bacteria"
# ═══════════════════════════════════════════════════════════════════════════════


rule icefinder2:
    """Detect ICEs/IMEs on the Medaka-polished metagenomic assembly."""
    input:
        assembly = "results/{sample}/medaka/{sample}_polished.fasta"
    output:
        results_dir = directory("results/{sample}/icefinder2")
    params:
        icefinder_dir = config["icefinder2"]["install_dir"],
        mode          = config["icefinder2"].get("mode", "Metagenome"),
        kraken_db     = config["kraken_db"]
    log:
        "logs/icefinder2/{sample}.log"
    benchmark:
        "benchmarks/icefinder2/{sample}.txt"
    conda:
        "../envs/icefinder2.yaml"
    threads: config["resources"]["threads"]
    shell:
        """
        mkdir -p {output.results_dir}
        echo "[$(date)] ICEfinder2 - {wildcards.sample}" > {log}
        echo "[$(date)] Mode: {params.mode}" >> {log}
        echo "[$(date)] Kraken2 db: {params.kraken_db}" >> {log}

        cd {params.icefinder_dir}
        python ICEfinder2.py \
            -i {input.assembly} \
            -t {params.mode} \
            >> {log} 2>&1

        mv {params.icefinder_dir}/result/* {output.results_dir}/ 2>> {log} || true

        echo "[$(date)] ICEfinder2 termine pour {wildcards.sample}" >> {log}
        """
