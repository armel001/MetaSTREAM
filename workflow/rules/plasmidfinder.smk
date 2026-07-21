# ─── Plasmid Replicon Typing — PlasmidFinder ─────────────────────────────────
# Pipeline: assemblage poli (Medaka) → BLASTn contre plasmidfinder_db
# Input: assemblage poli (post-Medaka), PAS le brouillon Flye brut
#   → les motifs rep sont courts (~200-500 pb) et sensibles aux erreurs
#     residuelles d'indels non corrigees par le polissage consensus.
# Reference: Carattoli A. et al. 2014, Antimicrob Agents Chemother.
#            "PlasmidFinder and pMLST: in silico detection and typing
#            of plasmids"
# ═══════════════════════════════════════════════════════════════════════════════


rule plasmidfinder:
    """Type plasmid replicons (Inc groups) on the Medaka-polished assembly."""
    input:
        assembly = "results/{sample}/medaka/{sample}_polished.fasta"
    output:
        results_tab = "results/{sample}/plasmidfinder/{sample}_results_tab.tsv",
        results_txt = "results/{sample}/plasmidfinder/{sample}_results.txt",
        json        = "results/{sample}/plasmidfinder/{sample}_data.json"
    params:
        outdir    = "results/{sample}/plasmidfinder",
        tmp_dir   = "results/{sample}/plasmidfinder/tmp",
        db        = config["plasmidfinder"]["db"],
        databases = config["plasmidfinder"].get("databases", "enterobacteriaceae,gram_positive"),
        method    = config["plasmidfinder"].get("method_path", "blastn"),
        min_cov   = config["plasmidfinder"].get("min_cov", 0.60),
        threshold = config["plasmidfinder"].get("threshold", 0.95),
        extended  = "-x"
    log:
        "logs/plasmidfinder/{sample}.log"
    benchmark:
        "benchmarks/plasmidfinder/{sample}.txt"
    conda:
        "../envs/plasmidfinder.yaml"
    threads: config["resources"]["threads"]
    shell:
        """
        mkdir -p {params.outdir} {params.tmp_dir}
        echo "[$(date)] PlasmidFinder — {wildcards.sample}" > {log}
        echo "[$(date)] Input (polished): {input.assembly}" >> {log}
        echo "[$(date)] DB path: {params.db}" >> {log}
        echo "[$(date)] Databases: {params.databases}" >> {log}
        echo "[$(date)] Method: {params.method}, min_cov={params.min_cov}, threshold={params.threshold}" >> {log}

        plasmidfinder.py \
            -i {input.assembly} \
            -o {params.outdir} \
            -tmp {params.tmp_dir} \
            -p {params.db} \
            -d {params.databases} \
            -mp {params.method} \
            -l {params.min_cov} \
            -t {params.threshold} \
            {params.extended} \
            2>> {log}

        mv {params.outdir}/results_tab.tsv  {output.results_tab}
        mv {params.outdir}/results.txt      {output.results_txt}
        mv {params.outdir}/data.json        {output.json}

        NHITS=$(tail -n +2 {output.results_tab} | wc -l || echo 0)
        echo "[$(date)] Replicons plasmidiques detectes : ${{NHITS}}" >> {log}

        rm -rf {params.tmp_dir}
        echo "[$(date)] PlasmidFinder termine pour {wildcards.sample}" >> {log}
        """
