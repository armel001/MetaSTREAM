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

# ─── 1b. Reference model filtering (RGI BWT only) ────────────────────────────
#
# RGI BWT maps reads against CARD reference sequences using several model
# types: protein homolog model, protein variant model, rRNA gene variant
# model, and protein overexpression model. SNP-confirmation models (variant,
# rRNA, overexpression) require precise mutation calling from read alignment,
# which is unreliable on fragmented, non-assembled metagenomic reads —
# especially with Nanopore error profiles. The CARD/RGI documentation
# recommends restricting BWT-based analyses to "protein homolog model" hits
# for metagenomic short/long-read data lacking high-quality assemblies.
# This filter is applied defensively even when the current dataset already
# satisfies it, to ensure robustness against future CARD updates or new
# samples introducing SNP-based model hits.

if (tool == "rgi_bwt") {
  cat("\n[1b/8] Filtering RGI BWT hits by reference model type...\n")

  col_model <- "reference_model_type"

  if (!col_model %in% colnames(rgi_all)) {
    stop("ERROR: column '", col_model, "' not found in RGI BWT input. ",
         "Available: ", paste(colnames(rgi_all), collapse = ", "))
  }

  n_before <- nrow(rgi_all)
  model_counts_before <- rgi_all %>% count(.data[[col_model]], name = "n")

  cat("  Reference model types found (before filtering):\n")
  for (i in seq_len(nrow(model_counts_before))) {
    cat(sprintf("    %-30s : %d rows\n",
                model_counts_before[[col_model]][i],
                model_counts_before$n[i]))
  }

  rgi_all <- rgi_all %>%
    filter(.data[[col_model]] == "protein homolog model")

  n_after   <- nrow(rgi_all)
  n_removed <- n_before - n_after

  cat(sprintf("  ✓ Kept 'protein homolog model' hits: %d / %d rows (%.1f%%)\n",
              n_after, n_before, 100 * n_after / n_before))
  if (n_removed > 0) {
    cat(sprintf("  ✓ Removed %d rows (SNP-based / non-homolog models)\n",
                n_removed))
  } else {
    cat("  ✓ No rows removed — all hits were already 'protein homolog model'\n")
  }
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

# Joindre drug_class — fusionner les valeurs multiples par (sample_id,
# best_hit_aro) AVANT la jointure. Nécessaire car un renommage manuel de gène
# peut fusionner plusieurs entrées CARD d'origine ayant des drug_class
# différents ; sans cette fusion, distinct() produirait des clés dupliquées
# côté drug_lookup, et le left_join ferait un fan-out qui casse pivot_wider
# en aval (colonne count basculant en liste au lieu de numérique).
if (col_drug %in% colnames(rgi_all)) {
  drug_lookup <- rgi_all %>%
    select(sample_id,
           best_hit_aro = !!sym(col_aro),
           drug_class   = !!sym(col_drug)) %>%
    distinct() %>%
    group_by(sample_id, best_hit_aro) %>%
    summarise(drug_class = paste(sort(unique(drug_class)), collapse = "; "),
              .groups = "drop")

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

# ─── 8. Family prefix aggregation (literature-standard short labels) ────────
#
# Replaces heuristic regex-parsing of the ARO gene name (best_hit_aro), which
# was unreliable: species-prefixed gene names (e.g. "Klebsiella pneumoniae
# KpnE") were misparsed as bacterial genus names, and case-sensitive matching
# fragmented single families across multiple labels (e.g. ErmF/ErmB/ErmX vs
# "erm" in the known-prefix list).
#
# Instead, this maps directly from the official CARD `amr_gene_family`
# classification — which is assigned independently of how the ARO term is
# phrased — to a standardized short label matching common usage in the AMR
# literature (e.g. blaTEM, blaOXA, qnr, sul, dfr, tet, erm, mph, AAC, ANT,
# APH, RND efflux, MFS efflux...).

cat("\n[8/8] Creating family prefix aggregation (CARD family → standard short label)...\n")

# Manual mapping table — defined once, used both for mapping and for
# unmapped-value detection below.
card_family_manual_map <- c(
  "aminoglycoside bifunctional resistance protein"                          = "AAC/APH bifunctional",
  "16S rRNA methyltransferase (G1405)"                                       = "16S-RMTase",
  "Erm 23S ribosomal RNA methyltransferase"                                  = "Erm",
  "non-erm 23S ribosomal RNA methyltransferase (G748)"                      = "non-Erm 23S-MTase",
  "Cfr 23S ribosomal RNA methyltransferase"                                  = "Cfr",
  "macrolide esterase"                                                       = "Ere",
  "macrolide phosphotransferase (MPH)"                                       = "MPH",
  "msr-type ABC-F protein"                                                   = "Msr",
  "lsa-type ABC-F protein"                                                   = "Lsa",
  "lincosamide nucleotidyltransferase (LNU)"                                 = "Lnu",
  "tetracycline inactivation enzyme"                                         = "Tet (enzymatic)",
  "tetracycline-resistant ribosomal protection protein"                      = "Tet (RPP)",
  "sulfonamide resistant sul"                                                = "Sul",
  "trimethoprim resistant dihydrofolate reductase dfr"                       = "Dfr",
  "quinolone resistance protein (qnr)"                                       = "Qnr",
  "ATP-binding cassette (ABC) antibiotic efflux pump"                        = "ABC efflux",
  "major facilitator superfamily (MFS) antibiotic efflux pump"               = "MFS efflux",
  "multidrug and toxic compound extrusion (MATE) transporter"                = "MATE efflux",
  "resistance-nodulation-cell division (RND) antibiotic efflux pump"         = "RND efflux",
  "small multidrug resistance (SMR) antibiotic efflux pump"                  = "SMR efflux",
  "General Bacterial Porin with reduced permeability to beta-lactams"        = "Porin (beta-lactam)",
  "General Bacterial Porin with reduced permeability to peptide antibiotics" = "Porin (peptide)",
  "Outer Membrane Porin (Opr)"                                               = "Opr",
  "intrinsic colistin resistant phosphoethanolamine transferase"             = "Intrinsic PET",
  "MCR phosphoethanolamine transferase"                                      = "MCR",
  "pmr phosphoethanolamine transferase"                                      = "Pmr",
  "glycopeptide resistance gene cluster"                                     = "Van cluster",
  "Van ligase"   = "Van", "vanH" = "Van", "vanT" = "Van",
  "vanU" = "Van", "vanW" = "Van", "vanXY" = "Van", "vanY" = "Van",
  "defensin resistant mprF"                                                  = "MprF",
  "Intrinsic peptide antibiotic resistant Lps"                               = "Lps",
  "rifampin ADP-ribosyltransferase (Arr)"                                    = "Arr",
  "rifamycin-resistant beta-subunit of RNA polymerase (rpoB)"                = "RpoB",
  "Bleomycin resistant protein"                                              = "Ble",
  "fosfomycin thiol transferase"                                             = "Fos",
  "kdpDE"                                                                    = "KdpDE",
  "methicillin resistant PBP2"                                               = "PBP2 (mecA-like)",
  "RbpA bacterial RNA polymerase-binding protein"                            = "RbpA",
  "streptothricin acetyltransferase (SAT)"                                   = "SAT",
  "tunicamycin resistance protein"                                           = "TmrB",
  "undecaprenyl pyrophosphate related proteins"                              = "UppP",
  "chloramphenicol acetyltransferase (CAT)"                                  = "CAT",
  "nitroimidazole reductase"                                                 = "Nim"
)

map_family_to_prefix <- function(family, manual_map = card_family_manual_map) {
  f <- trimws(family)
  if (is.na(f) || f == "" || f == "NA") return("Unclassified")

  # Rule 1 — Beta-lactamases: "XXX[-N-like] beta-lactamase" / "Beta-lactamase"
  if (grepl("beta-lactamase$", f, ignore.case = TRUE)) {
    prefix <- sub("(?i)(-\\d+)?-like\\s+beta-lactamase$|\\s+beta-lactamase$",
                  "", f, perl = TRUE)
    prefix <- trimws(prefix)

    if (tolower(prefix) == "ampc-type") return("AmpC")
    if (nchar(prefix) == 0) return("Other beta-lactamase")

    tokens <- strsplit(prefix, "\\s+")[[1]]
    if (length(tokens) > 2) {
      last <- tokens[length(tokens)]
      if (last == tolower(last)) {
        last <- paste0(toupper(substr(last, 1, 1)), substr(last, 2, nchar(last)))
      }
      return(last)
    }
    return(prefix)
  }

  # Rule 2 — Aminoglycoside-modifying enzymes: AAC/ANT/APH(...)
  m <- regmatches(f, regexpr("^(AAC|ANT|APH)", f))
  if (length(m) > 0 && nchar(m[1]) > 0) return(m[1])

  # Rule 3 — Manual mapping table
  if (f %in% names(manual_map)) return(unname(manual_map[f]))

  # Unmapped — kept as the original CARD string, flagged below for review
  return(f)
}

# Build family_prefix_abundance from rgi_all (one row per hit), splitting
# multi-value amr_gene_family cells exactly as in step 7, then mapping each
# split value to its standardized short label before aggregating.

if (col_family %in% colnames(rgi_all)) {

  rgi_family_split <- rgi_all %>%
    separate_rows(!!sym(col_family), sep = "; ") %>%
    mutate(amr_gene_family_raw = str_trim(.data[[col_family]]),
           family_prefix        = vapply(amr_gene_family_raw,
                                          map_family_to_prefix,
                                          character(1)))

  if (!is.null(col_reads) && col_reads %in% colnames(rgi_family_split)) {
    family_prefix_abundance <- rgi_family_split %>%
      group_by(sample_id, family_prefix) %>%
      summarise(
        n_variants   = n_distinct(.data[[col_aro]]),
        count        = n(),
        mapped_reads = sum(.data[[col_reads]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      left_join(sequencing_stats_gb %>% select(sample_id, total_bases_gb),
                by = "sample_id") %>%
      mutate(normalized_abundance = mapped_reads / total_bases_gb)
  } else {
    family_prefix_abundance <- rgi_family_split %>%
      group_by(sample_id, family_prefix) %>%
      summarise(
        n_variants = n_distinct(.data[[col_aro]]),
        count      = n(),
        .groups = "drop"
      ) %>%
      left_join(sequencing_stats_gb %>% select(sample_id, total_bases_gb),
                by = "sample_id") %>%
      mutate(normalized_abundance = count / total_bases_gb)
  }

  write_tsv(family_prefix_abundance, output_family_prefix)

  # Flag any family value that did not match a known rule and was returned
  # as-is — useful to detect new CARD families introduced by updates or new
  # samples that require adding an entry to card_family_manual_map.
  is_betalactamase <- grepl("beta-lactamase$", rgi_family_split$amr_gene_family_raw,
                            ignore.case = TRUE)
  is_amino_enzyme  <- grepl("^(AAC|ANT|APH)", rgi_family_split$amr_gene_family_raw)
  is_manual_mapped <- rgi_family_split$amr_gene_family_raw %in% names(card_family_manual_map)
  is_na_or_empty   <- is.na(rgi_family_split$amr_gene_family_raw) |
                       trimws(rgi_family_split$amr_gene_family_raw) %in% c("", "NA")

  unmapped <- unique(rgi_family_split$amr_gene_family_raw[
    !is_betalactamase & !is_amino_enzyme & !is_manual_mapped & !is_na_or_empty
  ])

  if (length(unmapped) > 0) {
    cat("  ⚠ Unmapped CARD families (kept as-is, review recommended):\n")
    for (u in unmapped) cat("    -", u, "\n")
  } else {
    cat("  ✓ All CARD families successfully mapped to short labels\n")
  }

  cat("  ✓ family_prefix_abundance (",
      n_distinct(family_prefix_abundance$family_prefix), " standardized families, ",
      nrow(family_prefix_abundance), " rows)\n", sep = "")

} else {
  cat("  ⚠ amr_gene_family not found — empty file\n")
  write_tsv(tibble(), output_family_prefix)
}

# Preview top families
top_fam <- family_prefix_abundance %>%
  group_by(family_prefix) %>%
  summarise(total = sum(normalized_abundance), .groups = "drop") %>%
  arrange(desc(total)) %>%
  head(15)

cat("\n  Top 15 standardized families by normalized abundance:\n")
for (i in seq_len(nrow(top_fam))) {
  cat(sprintf("    %2d. %-25s %.2f\n",
              i, top_fam$family_prefix[i], top_fam$total[i]))
}

# ─── Summary ─────────────────────────────────────────────────────────────────

cat("\n", rep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETED  [", tool, "]\n", sep = "")
cat(rep("=", 70), "\n", sep = "")
cat("Output directory: results/r_analysis/", tool, "/\n", sep = "")
if (tool == "rgi_bwt") {
  cat("  • Reference model filter   : protein homolog model only\n")
}
cat("  • Counts and normalization : 4 files\n")
cat("  • Matrices                 : 4 files\n")
cat("  • Functional aggregations  : 3 files\n")
cat("  • Family prefix            : 1 file\n")
cat("  • Total                    : 12 files\n")
cat(rep("=", 70), "\n\n", sep = "")

sink(type = "message")
sink()
close(log)
