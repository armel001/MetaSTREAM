# ============================================================================
# Taxonomic analysis rules for multi-level classification
# ============================================================================

TAXONOMIC_LEVELS = ["S", "G", "F", "P"]

# ============================================================================
# Combine samples into abundance matrices
# ============================================================================

rule combine_bracken_samples:
    input:
        expand("results/{sample}/bracken/{sample}_{{level}}.txt", 
               sample=config["samples_id"])
    output:
        abundance_matrix = "results/taxonomy/abundance_matrix_{level}.tsv",
        metadata = "results/taxonomy/sample_metadata_{level}.tsv"
    params:
        samples = config["samples_id"]
    wildcard_constraints:
        level = "S|G|F|P"
    conda:
        "../envs/taxo_analysis.yaml"
    log:
        "logs/combine/combine_{level}.log"
    script:
        "../scripts/combine_bracken.py"

# ============================================================================
# Standard filtering (homogeneous samples)
# ============================================================================

rule filter_low_abundance:
    input:
        matrix = "results/taxonomy/abundance_matrix_{level}.tsv"
    output:
        filtered = "results/taxonomy/abundance_matrix_{level}_filtered.tsv",
        removed = "results/taxonomy/removed_taxa_{level}.txt"
    params:
        min_reads = config.get("min_reads_per_taxon", 100),
        min_samples = config.get("min_samples_present", 2),
        min_abundance = config.get("min_relative_abundance", 0.001)
    wildcard_constraints:
        level = "S|G|F|P"
    conda:
        "../envs/taxo_analysis.yaml"
    log:
        "logs/filter/filter_{level}.log"
    script:
        "../scripts/filter_abundance.py"

# ============================================================================
# Independent filtering (heterogeneous samples)
# ============================================================================

rule filter_low_abundance_independent:
    input:
        matrix = "results/taxonomy/abundance_matrix_{level}.tsv"
    output:
        filtered = "results/taxonomy/abundance_matrix_{level}_filtered_independent.tsv",
        removed = "results/taxonomy/removed_taxa_{level}_independent.txt",
        stats = "results/taxonomy/filtering_stats_{level}_independent.txt"
    params:
        min_reads_per_sample = config.get("min_reads_per_sample", 50),
        min_rel_abundance_per_sample = config.get("min_rel_abundance_per_sample", 0.001),
        min_samples_present = config.get("min_samples_independent", 1),
        absolute_min_reads = config.get("absolute_min_total_reads", 10)
    wildcard_constraints:
        level = "S|G|F|P"
    conda:
        "../envs/taxo_analysis.yaml"
    log:
        "logs/filter/filter_independent_{level}.log"
    script:
        "../scripts/filter_abundance_independent.py"

# ============================================================================
# Normalize abundance (standard)
# ============================================================================

rule normalize_abundance:
    input:
        matrix = "results/taxonomy/abundance_matrix_{level}_filtered.tsv"
    output:
        relative = "results/taxonomy/abundance_matrix_{level}_relative.tsv",
        clr = "results/taxonomy/abundance_matrix_{level}_clr.tsv"
    wildcard_constraints:
        level = "S|G|F|P"
    conda:
        "../envs/taxo_analysis.yaml"
    log:
        "logs/normalize/normalize_{level}.log"
    script:
        "../scripts/normalize_abundance.py"

# ============================================================================
# Normalize abundance (independent) - with sequencing stats
# ============================================================================

rule normalize_abundance_independent:
    input:
        matrix = "results/taxonomy/abundance_matrix_{level}_filtered_independent.tsv",
        stats = "results/stats/sequencing_stats.tsv"  # ← Ajouter dépendance
    output:
        relative = "results/taxonomy/abundance_matrix_{level}_relative_independent.tsv",
        clr = "results/taxonomy/abundance_matrix_{level}_clr_independent.tsv"
    wildcard_constraints:
        level = "S|G|F|P"
    conda:
        "../envs/taxo_analysis.yaml"
    log:
        "logs/normalize/normalize_independent_{level}.log"
    script:
        "../scripts/normalize_abundance.py"


# ============================================================================
# Alpha diversity (standard)
# ============================================================================

rule calculate_alpha_diversity:
    input:
        matrix = "results/taxonomy/abundance_matrix_S.tsv"
    output:
        diversity = "results/taxonomy/alpha_diversity.tsv"
    conda:
        "../envs/r_diversity.yaml"
    log:
        "logs/diversity/alpha_diversity.log"
    script:
        "../scripts/alpha_diversity.R"

# ============================================================================
# Alpha diversity (independent)
# ============================================================================

rule calculate_alpha_diversity_independent:
    input:
        matrix = "results/taxonomy/abundance_matrix_S_filtered_independent.tsv"
    output:
        diversity = "results/taxonomy/alpha_diversity_independent.tsv"
    conda:
        "../envs/r_diversity.yaml"
    log:
        "logs/diversity/alpha_diversity_independent.log"
    script:
        "../scripts/alpha_diversity.R"

# ============================================================================
# Summary statistics (standard)
# ============================================================================

rule summary_statistics:
    input:
        kraken_reports = expand("results/{sample}/kraken/{sample}_report.txt", 
                               sample=config["samples_id"]),
        diversity = "results/taxonomy/alpha_diversity.tsv",
        matrices = expand("results/taxonomy/abundance_matrix_{level}_filtered.tsv",
                         level=TAXONOMIC_LEVELS)
    output:
        summary = "results/summary/taxonomy_summary.txt"
    conda:
        "../envs/taxo_analysis.yaml"
    log:
        "logs/summary/taxonomy_summary.log"
    script:
        "../scripts/taxo_summary.py"

# ============================================================================
# Summary statistics (independent)
# ============================================================================

rule summary_statistics_independent:
    input:
        kraken_reports = expand("results/{sample}/kraken/{sample}_report.txt", 
                               sample=config["samples_id"]),
        diversity = "results/taxonomy/alpha_diversity_independent.tsv",
        matrices = expand("results/taxonomy/abundance_matrix_{level}_filtered_independent.tsv",
                         level=TAXONOMIC_LEVELS),
        filtering_stats = expand("results/taxonomy/filtering_stats_{level}_independent.txt",
                                level=TAXONOMIC_LEVELS)
    output:
        summary = "results/summary/taxonomy_summary_independent.txt"
    conda:
        "../envs/taxo_analysis.yaml"
    log:
        "logs/summary/taxonomy_summary_independent.log"
    script:
        "../scripts/taxo_summary_independent.py"
