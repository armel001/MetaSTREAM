#!/usr/bin/env Rscript

library(tidyverse)

# ─── Snakemake params ────────────────────────────────────────────────────────

input_rgi    <- snakemake@input$rgi_aggregated
input_stats  <- snakemake@input$sequencing_stats
tool         <- snakemake@params$tool
log_file     <- snakemake@log[[1]]

output_counts              <- snakemake@output$arg_counts
output_relative            <- snakemake@output$arg_relative
output_normalized          <- snakemake@output$arg_normalized
output_presence            <- snakemake@output$arg_presence
output_matrix_counts       <- snakemake@output$arg_matrix_counts
output_matrix_relative     <- snakemake@output$arg_matrix_relative
output_matrix_normalized   <- snakemake@output$arg_matrix_normalized
output_matrix_presence     <- snakemake@output$arg_matrix_presence
output_drug_class          <- snakemake@output$drug_class_abundance
output_mechanism           <- snakemake@output$mechanism_abundance
output_family              <- snakemake@output$family_abundance
output_family_prefix       <- snakemake@output$family_prefix_abundance

# ─── Column mapping ──────────────────────────────────────────────────────────

if (tool == "rgi") {
  col_aro    <- "best_hit_aro"
  col_drug   <- "drug_class"
  col_mech   <- "resistance_mechanism"
  col_family <- "amr_gene_family"
  col_reads  <- NULL
} else {
  col_aro    <- "aro_term"
  col_drug   <- "drug_class"
  col_mech   <- "resistance_mechanism"
  col_family <- "amr_gene_family"
  col_reads  <- "all_mapped_reads"
}

# ─── Logger ──────────────────────────────────────────────────────────────────

log <- file(log_file, open = "wt")
sink(log)
sink(log, type = "message")

cat(rep("=", 70), "\n", sep = "")
cat("R Analysis: Normalization and Matrices  [", tool, "]\n", sep = "")
cat(rep("=", 70), "\n", sep = "")

# ─── 1. Load ─────────────────────────────────────────────────────────────────

cat("\n[1/8] Loading data...\n")

rgi_all          <- read_tsv(input_rgi,   col_types = cols(), show_col_types = FALSE)
sequencing_stats <- read_tsv(input_stats, col_types = cols(), show_col_types = FALSE)

cat("  ✓ ", tool, ": ", nrow(rgi_all), " observations\n", sep = "")
cat("  ✓ Stats: ", nrow(sequencing_stats), " samples\n", sep = "")
cat("  ✓ ARO column: ", col_aro, "\n", sep = "")

if (!col_aro %in% colnames(rgi_all)) {
  stop("ERROR: ARO column '", col_aro, "' not found. ",
       "Available: ", paste(colnames(rgi_all), collapse = ", "))
}

# ─── 2. Raw counts ───────────────────────────────────────────────────────────

cat("\n[2/8] Computing raw counts...\n")

if (!is.null(col_reads) && col_reads %in% colnames(rgi_all)) {
  arg_counts <- rgi_all %>%
    group_by(sample_id, aro_term = .data[[col_aro]]) %>%
    summarise(
      count        = n(),
      mapped_reads = sum(.data[[col_reads]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(best_hit_aro = aro_term)
  cat("  ✓ Using all_mapped_reads as abundance metric\n")
} else {
  arg_counts <- rgi_all %>%
    count(sample_id, !!sym(col_aro), name = "count") %>%
    rename(best_hit_aro = !!sym(col_aro))
}

# Joindre drug_class
if (col_drug %in% colnames(rgi_all)) {
  drug_lookup <- rgi_all %>%
    select(sample_id,
           best_hit_aro = !!sym(col_aro),
           drug_class   = !!sym(col_drug)) %>%
    distinct()

  arg_counts <- arg_counts %>%
    left_join(drug_lookup, by = c("sample_id", "best_hit_aro"))

  cat("  ✓ drug_class joined to arg_counts\n")
}

# Sauvegarder après le join   ← déplacé ici
write_tsv(arg_counts, output_counts)
cat("  ✓ Saved: ", basename(output_counts),
    " (", nrow(arg_counts), " rows)\n", sep = "")


# ─── 3. Relative abundances ──────────────────────────────────────────────────

cat("\n[3/8] Computing relative abundances...\n")

abundance_col <- if (!is.null(col_reads) && "mapped_reads" %in% colnames(arg_counts)) {
  "mapped_reads"
} else {
  "count"
}

arg_relative <- arg_counts %>%
  group_by(sample_id) %>%
  mutate(
    total_args         = sum(.data[[abundance_col]]),
    relative_abundance = (.data[[abundance_col]] / total_args) * 100
  ) %>%
  ungroup()

write_tsv(arg_relative, output_relative)
cat("  ✓ Saved: ", basename(output_relative),
    " (metric: ", abundance_col, ")\n", sep = "")

# ─── 4. Normalized copy number (per Gb) ──────────────────────────────────────

cat("\n[4/8] Computing normalized copy number (per Gb)...\n")

sequencing_stats_gb <- sequencing_stats %>%
  mutate(total_bases_gb = total_bases / 1e9)

cat("  Dataset sizes (Gb):\n")
for (i in seq_len(nrow(sequencing_stats_gb))) {
  cat(sprintf("    %s: %.3f Gb\n",
              sequencing_stats_gb$sample_id[i],
              sequencing_stats_gb$total_bases_gb[i]))
}

if (!is.null(col_reads) && "mapped_reads" %in% colnames(arg_counts)) {
  arg_normalized <- arg_counts %>%
    left_join(sequencing_stats_gb %>% select(sample_id, total_bases_gb),
              by = "sample_id") %>%
    mutate(normalized_copy_number = mapped_reads / total_bases_gb)
  cat("  ✓ Formula: mapped_reads / total_bases_gb\n")
} else {
  arg_normalized <- arg_counts %>%
    left_join(sequencing_stats_gb %>% select(sample_id, total_bases_gb),
              by = "sample_id") %>%
    mutate(normalized_copy_number = count / total_bases_gb)
  cat("  ✓ Formula: count / total_bases_gb\n")
}

write_tsv(arg_normalized, output_normalized)
cat("  ✓ Saved: ", basename(output_normalized), "\n", sep = "")

# ─── 5. Presence/absence ─────────────────────────────────────────────────────

cat("\n[5/8] Creating presence/absence table...\n")

arg_presence <- arg_counts %>%
  mutate(presence = 1) %>%
  select(sample_id, best_hit_aro, presence)

write_tsv(arg_presence, output_presence)
cat("  ✓ Saved: ", basename(output_presence), "\n", sep = "")

# ─── 6. Matrices ─────────────────────────────────────────────────────────────

cat("\n[6/8] Creating matrices...\n")

arg_matrix_counts <- arg_counts %>%
  select(sample_id, best_hit_aro, count) %>%
  pivot_wider(names_from = best_hit_aro, values_from = count,
              values_fill = 0)
write_tsv(arg_matrix_counts, output_matrix_counts)
cat("  ✓ arg_matrix_counts: ", nrow(arg_matrix_counts), " × ",
    ncol(arg_matrix_counts) - 1, "\n", sep = "")

arg_matrix_relative <- arg_relative %>%
  select(sample_id, best_hit_aro, relative_abundance) %>%
  pivot_wider(names_from = best_hit_aro, values_from = relative_abundance,
              values_fill = 0)
write_tsv(arg_matrix_relative, output_matrix_relative)
cat("  ✓ arg_matrix_relative\n")

arg_matrix_normalized <- arg_normalized %>%
  select(sample_id, best_hit_aro, normalized_copy_number) %>%
  pivot_wider(names_from = best_hit_aro, values_from = normalized_copy_number,
              values_fill = 0)
write_tsv(arg_matrix_normalized, output_matrix_normalized)
cat("  ✓ arg_matrix_normalized (copy number per Gb)\n")

arg_matrix_presence <- arg_presence %>%
  pivot_wider(names_from = best_hit_aro, values_from = presence,
              values_fill = 0)
write_tsv(arg_matrix_presence, output_matrix_presence)
cat("  ✓ arg_matrix_presence\n")

# ─── 7. Functional aggregations ──────────────────────────────────────────────

cat("\n[7/8] Creating functional aggregations...\n")

# Drug class
if (col_drug %in% colnames(rgi_all)) {
  drug_class_abundance <- rgi_all %>%
    separate_rows(!!sym(col_drug), sep = "; ") %>%
    mutate(drug_class = str_trim(.data[[col_drug]])) %>%
    count(sample_id, drug_class, name = "count") %>%
    group_by(sample_id) %>%
    mutate(total = sum(count),
           relative_abundance = (count / total) * 100) %>%
    ungroup()
  write_tsv(drug_class_abundance, output_drug_class)
  cat("  ✓ drug_class_abundance (",
      n_distinct(drug_class_abundance$drug_class), " classes)\n", sep = "")
} else {
  cat("  ⚠ drug_class not found — empty file\n")
  write_tsv(tibble(), output_drug_class)
}

# Mechanism
if (col_mech %in% colnames(rgi_all)) {
  mechanism_abundance <- rgi_all %>%
    separate_rows(!!sym(col_mech), sep = "; ") %>%
    mutate(resistance_mechanism = str_trim(.data[[col_mech]])) %>%
    count(sample_id, resistance_mechanism, name = "count") %>%
    group_by(sample_id) %>%
    mutate(total = sum(count),
           relative_abundance = (count / total) * 100) %>%
    ungroup()
  write_tsv(mechanism_abundance, output_mechanism)
  cat("  ✓ mechanism_abundance (",
      n_distinct(mechanism_abundance$resistance_mechanism),
      " mechanisms)\n", sep = "")
} else {
  cat("  ⚠ resistance_mechanism not found — empty file\n")
  write_tsv(tibble(), output_mechanism)
}

# AMR gene family
if (col_family %in% colnames(rgi_all)) {
  family_abundance <- rgi_all %>%
    separate_rows(!!sym(col_family), sep = "; ") %>%
    mutate(amr_gene_family = str_trim(.data[[col_family]])) %>%
    count(sample_id, amr_gene_family, name = "count") %>%
    group_by(sample_id) %>%
    mutate(total = sum(count),
           relative_abundance = (count / total) * 100) %>%
    ungroup()
  write_tsv(family_abundance, output_family)
  cat("  ✓ family_abundance (",
      n_distinct(family_abundance$amr_gene_family), " families)\n", sep = "")
} else {
  cat("  ⚠ amr_gene_family not found — empty file\n")
  write_tsv(tibble(), output_family)
}

# ─── 8. Family prefix aggregation ────────────────────────────────────────────

cat("\n[8/8] Creating family prefix aggregation...\n")

# Extraire le préfixe de famille
# CTX-M-15 → CTX-M | OXA-101 → OXA | AAC(6')-Ib7 → AAC(6')
extract_family_prefix <- function(name) {
  known_prefixes <- c(
    "CTX-M", "OXA", "AAC", "ANT", "APH", "CMY", "SHV", "TEM",
    "NDM", "KPC", "VIM", "IMP", "GES", "PER", "VEB", "CFX",
    "CfxA", "CfiA", "CARB", "AER", "CAE", "sul", "tet", "erm",
    "msr", "nim", "aad", "aph", "dfr", "qnr", "mcr", "van",
    "bla", "mph", "mef", "erm", "lnu", "vga", "eat", "sal",
    "optrA", "cfr", "poxtA"
  )
  for (prefix in known_prefixes) {
    if (startsWith(name, prefix)) {
      return(prefix)
    }
  }
  # Regex générique : extraire avant le premier tiret + chiffre
  m <- regmatches(name,
                  regexpr("^[A-Za-z\\(\\)\\'\"]+(?:-[A-Za-z]+)?",
                          name, perl = TRUE))
  if (length(m) > 0 && nchar(m) > 0) {
    return(sub("-$", "", m))
  }
  # Fallback : premier mot
  return(strsplit(name, " ")[[1]][1])
}

arg_counts_prefix <- arg_counts %>%
  mutate(family_prefix = sapply(best_hit_aro, extract_family_prefix))

# Agréger par sample + famille
if (!is.null(col_reads) && "mapped_reads" %in% colnames(arg_counts)) {
  family_prefix_abundance <- arg_counts_prefix %>%
    group_by(sample_id, family_prefix) %>%
    summarise(
      n_variants           = n_distinct(best_hit_aro),
      count                = sum(count),
      mapped_reads         = sum(mapped_reads),
      .groups = "drop"
    ) %>%
    left_join(sequencing_stats_gb %>% select(sample_id, total_bases_gb),
              by = "sample_id") %>%
    mutate(normalized_abundance = mapped_reads / total_bases_gb)
} else {
  family_prefix_abundance <- arg_counts_prefix %>%
    group_by(sample_id, family_prefix) %>%
    summarise(
      n_variants = n_distinct(best_hit_aro),
      count      = sum(count),
      .groups = "drop"
    ) %>%
    left_join(sequencing_stats_gb %>% select(sample_id, total_bases_gb),
              by = "sample_id") %>%
    mutate(normalized_abundance = count / total_bases_gb)
}

write_tsv(family_prefix_abundance, output_family_prefix)
cat("  ✓ family_prefix_abundance (",
    n_distinct(family_prefix_abundance$family_prefix),
    " families, ", nrow(family_prefix_abundance), " rows)\n", sep = "")

# Preview top families
top_fam <- family_prefix_abundance %>%
  group_by(family_prefix) %>%
  summarise(total = sum(normalized_abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  head(10)

cat("\n  Top 10 families by normalized abundance:\n")
for (i in seq_len(nrow(top_fam))) {
  cat(sprintf("    %2d. %-20s %.2f\n",
              i, top_fam$family_prefix[i], top_fam$total[i]))
}

# ─── Summary ─────────────────────────────────────────────────────────────────

cat("\n", rep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETED  [", tool, "]\n", sep = "")
cat(rep("=", 70), "\n", sep = "")
cat("Output directory: results/r_analysis/", tool, "/\n", sep = "")
cat("  • Counts and normalization : 4 files\n")
cat("  • Matrices                 : 4 files\n")
cat("  • Functional aggregations  : 3 files\n")
cat("  • Family prefix            : 1 file\n")
cat("  • Total                    : 12 files\n")
cat(rep("=", 70), "\n\n", sep = "")

sink(type = "message")
sink()
close(log)
