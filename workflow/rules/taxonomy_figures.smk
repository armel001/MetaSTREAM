rule plot_taxonomy_figures:
    """Generate taxonomic analysis HTML report."""
    input:
        genus_matrix   = fig_genus_matrix,
        species_matrix = fig_species_matrix,
        alpha_div      = fig_alpha_div
    output:
        report = "results/figures/taxonomy/taxonomy_report.html"
    params:
        top_genera  = config.get("fig_top_genera",    15),
        top_species = config.get("fig_top_species",   20),
        top_heatmap = config.get("fig_heatmap_top_n", 30),
        dpi         = config.get("fig_dpi",           150),
        pathogens   = config.get("priority_pathogens", [
            "Escherichia coli",
            "Klebsiella pneumoniae",
            "Pseudomonas aeruginosa",
            "Clostridioides difficile",
            "Clostridium perfringens",
            "Acinetobacter baumannii",
            "Enterococcus faecium"
        ])
    conda:
        "../envs/python_viz.yaml"
    log:
        "logs/figures/taxonomy_report.log"
    script:
        "../scripts/plot_taxonomy_report.py"
