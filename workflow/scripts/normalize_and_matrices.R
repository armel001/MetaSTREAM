#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# normalize_and_matrices.R — ARG normalization and abundance matrices
# Pipeline MetaSTREAM | RGI / RGI-BWT | CARD database
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(tidyverse))

# ── I/O ───────────────────────────────────────────────────────────────────────
input_rgi   <- snakemake@input$rgi_aggregated
input_stats <- snakemake@input$sequencing_stats
tool        <- snakemake@params$tool
log_file    <- snakemake@log[[1]]

output_counts            <- snakemake@output$arg_counts
output_relative          <- snakemake@output$arg_relative
output_normalized        <- snakemake@output$arg_normalized
output_presence          <- snakemake@output$arg_presence
output_matrix_counts     <- snakemake@output$arg_matrix_counts
output_matrix_relative   <- snakemake@output$arg_matrix_relative
output_matrix_normalized <- snakemake@output$arg_matrix_normalized
output_matrix_presence   <- snakemake@output$arg_matrix_presence
output_drug_class        <- snakemake@output$drug_class_abundance
output_mechanism         <- snakemake@output$mechanism_abundance
output_family            <- snakemake@output$family_abundance
output_family_prefix     <- snakemake@output$family_prefix_abundance

log <- file(log_file, open = "wt"); sink(log); sink(log, type = "message")

# ── Column mapping ────────────────────────────────────────────────────────────
if (tool == "rgi") {
  col_aro <- "best_hit_aro"; col_drug <- "drug_class"
  col_mech <- "resistance_mechanism"; col_family <- "amr_gene_family"
  col_reads <- NULL
} else {
  col_aro <- "aro_term"; col_drug <- "drug_class"
  col_mech <- "resistance_mechanism"; col_family <- "amr_gene_family"
  col_reads <- "all_mapped_reads"
}

# ── Helpers ───────────────────────────────────────────────────────────────────
save_tsv <- function(df, path, label = basename(path)) {
  write_tsv(df, path)
  cat(sprintf("  [OK] %-45s (%d rows)\n", label, nrow(df)))
  invisible(path)
}

# Aggregates a multi-value column (sep="; ") into relative abundance per sample
aggregate_functional <- function(df, col, out_col) {
  df %>%
    separate_rows(!!sym(col), sep = "; ") %>%
    mutate(!!out_col := str_trim(.data[[col]])) %>%
    count(sample_id, !!sym(out_col), name = "count") %>%
    group_by(sample_id) %>%
    mutate(total = sum(count), relative_abundance = count / total * 100) %>%
    ungroup()
}

cat(sprintf("══ NORMALIZATION & MATRICES — MetaSTREAM [%s] ══\n\n", tool))

# ── 1. Load ───────────────────────────────────────────────────────────────────
cat("[1/8] Loading data...\n")

rgi_all <- read_tsv(input_rgi,   show_col_types = FALSE)
stats   <- read_tsv(input_stats, show_col_types = FALSE) %>%
  mutate(total_bases_gb = total_bases / 1e9)

cat(sprintf("  %s: %d obs | Stats: %d samples | ARO col: %s\n",
            tool, nrow(rgi_all), nrow(stats), col_aro))

if (!col_aro %in% colnames(rgi_all))
  stop("ARO column '", col_aro, "' not found. Available: ",
       paste(colnames(rgi_all), collapse = ", "))

# ── 1b. RGI-BWT: restrict to protein homolog model ───────────────────────────
# SNP-confirmation models (variant, rRNA, overexpression) require precise
# mutation calling unreliable on fragmented metagenomic reads (esp. Nanopore).
# CARD/RGI documentation recommends protein homolog model only for BWT mode.
if (tool == "rgi_bwt") {
  cat("\n[1b/8] Filtering by reference model type (protein homolog model)...\n")
  col_model <- "reference_model_type"
  if (!col_model %in% colnames(rgi_all))
    stop("Column '", col_model, "' not found.")

  n_before <- nrow(rgi_all)
  cat("  Model types before filter:\n")
  rgi_all %>% count(.data[[col_model]]) %>%
    { cat(sprintf("    %-30s : %d\n", .[[col_model]], .$n)); . }

  rgi_all <- rgi_all %>% filter(.data[[col_model]] == "protein homolog model")
  n_after <- nrow(rgi_all)
  cat(sprintf("  Kept: %d / %d rows (%.1f%%)\n",
              n_after, n_before, 100 * n_after / n_before))
}

# ── 2. Raw counts ─────────────────────────────────────────────────────────────
cat("\n[2/8] Computing raw counts...\n")

if (!is.null(col_reads) && col_reads %in% colnames(rgi_all)) {
  arg_counts <- rgi_all %>%
    group_by(sample_id, best_hit_aro = .data[[col_aro]]) %>%
    summarise(count = n(), mapped_reads = sum(.data[[col_reads]], na.rm = TRUE),
              .groups = "drop")
  cat("  Using all_mapped_reads as abundance metric\n")
} else {
  arg_counts <- rgi_all %>%
    count(sample_id, !!sym(col_aro), name = "count") %>%
    rename(best_hit_aro = !!sym(col_aro))
}

# Join drug_class — merge multi-value entries before join to prevent fan-out
# in downstream pivot_wider (count column would become list instead of numeric)
if (col_drug %in% colnames(rgi_all)) {
  drug_lookup <- rgi_all %>%
    select(sample_id, best_hit_aro = !!sym(col_aro), drug_class = !!sym(col_drug)) %>%
    distinct() %>%
    group_by(sample_id, best_hit_aro) %>%
    summarise(drug_class = paste(sort(unique(drug_class)), collapse = "; "),
              .groups = "drop")
  arg_counts <- left_join(arg_counts, drug_lookup, by = c("sample_id", "best_hit_aro"))
}
save_tsv(arg_counts, output_counts)

# ── 3. Relative abundances ────────────────────────────────────────────────────
cat("\n[3/8] Computing relative abundances...\n")
abundance_col <- if (!is.null(col_reads) && "mapped_reads" %in% colnames(arg_counts))
  "mapped_reads" else "count"

arg_relative <- arg_counts %>%
  group_by(sample_id) %>%
  mutate(total_args = sum(.data[[abundance_col]]),
         relative_abundance = .data[[abundance_col]] / total_args * 100) %>%
  ungroup()
save_tsv(arg_relative, output_relative)

# ── 4. Normalized copy number (per Gb) ───────────────────────────────────────
cat("\n[4/8] Computing normalized copy number (per Gb)...\n")
cat(sprintf("    %s: %.3f Gb\n", stats$sample_id, stats$total_bases_gb))

norm_col <- if (!is.null(col_reads) && "mapped_reads" %in% colnames(arg_counts))
  "mapped_reads" else "count"

arg_normalized <- arg_counts %>%
  left_join(stats %>% select(sample_id, total_bases_gb), by = "sample_id") %>%
  mutate(normalized_copy_number = .data[[norm_col]] / total_bases_gb)
save_tsv(arg_normalized, output_normalized)

# ── 5. Presence / absence ─────────────────────────────────────────────────────
cat("\n[5/8] Creating presence/absence table...\n")
arg_presence <- arg_counts %>% mutate(presence = 1) %>%
  select(sample_id, best_hit_aro, presence)
save_tsv(arg_presence, output_presence)

# ── 6. Matrices ───────────────────────────────────────────────────────────────
cat("\n[6/8] Creating matrices...\n")

pivot_mat <- function(df, val_col, fill = 0)
  pivot_wider(select(df, sample_id, best_hit_aro, !!sym(val_col)),
              names_from = best_hit_aro, values_from = !!sym(val_col),
              values_fill = fill)

save_tsv(pivot_mat(arg_counts,      "count"),                    output_matrix_counts)
save_tsv(pivot_mat(arg_relative,    "relative_abundance"),       output_matrix_relative)
save_tsv(pivot_mat(arg_normalized,  "normalized_copy_number"),   output_matrix_normalized)
save_tsv(arg_presence %>% pivot_wider(names_from = best_hit_aro,
                                       values_from = presence, values_fill = 0),
         output_matrix_presence)

# ── 7. Functional aggregations ────────────────────────────────────────────────
cat("\n[7/8] Creating functional aggregations...\n")

for (lst in list(
  list(col = col_drug,    out = "drug_class",          path = output_drug_class),
  list(col = col_mech,    out = "resistance_mechanism", path = output_mechanism),
  list(col = col_family,  out = "amr_gene_family",     path = output_family)
)) {
  if (lst$col %in% colnames(rgi_all)) {
    agg <- aggregate_functional(rgi_all, lst$col, lst$out)
    save_tsv(agg, lst$path)
  } else {
    cat(sprintf("  [WARN] %s not found — empty file\n", lst$col))
    write_tsv(tibble(), lst$path)
  }
}

# ── 8. Family prefix aggregation (CARD family → standard short label) ─────────
# Maps CARD amr_gene_family strings to standardized short labels used in the
# AMR literature (blaTEM, qnr, sul, dfr, Erm, AAC, RND efflux, etc.).
# Direct CARD family mapping avoids unreliable regex parsing of ARO gene names
# (species-prefixed names like "Klebsiella pneumoniae KpnE" were misparsed).
cat("\n[8/8] Family prefix aggregation (CARD → standard short label)...\n")

CARD_FAMILY_MAP <- c(
  "aminoglycoside bifunctional resistance protein"                          = "AAC/APH bifunctional",
  "16S rRNA methyltransferase (G1405)"                                      = "16S-RMTase",
  "Erm 23S ribosomal RNA methyltransferase"                                 = "Erm",
  "non-erm 23S ribosomal RNA methyltransferase (G748)"                      = "non-Erm 23S-MTase",
  "Cfr 23S ribosomal RNA methyltransferase"                                 = "Cfr",
  "macrolide esterase"                                                      = "Ere",
  "macrolide phosphotransferase (MPH)"                                      = "MPH",
  "msr-type ABC-F protein"                                                  = "Msr",
  "lsa-type ABC-F protein"                                                  = "Lsa",
  "lincosamide nucleotidyltransferase (LNU)"                                = "Lnu",
  "tetracycline inactivation enzyme"                                        = "Tet (enzymatic)",
  "tetracycline-resistant ribosomal protection protein"                     = "Tet (RPP)",
  "sulfonamide resistant sul"                                               = "Sul",
  "trimethoprim resistant dihydrofolate reductase dfr"                      = "Dfr",
  "quinolone resistance protein (qnr)"                                      = "Qnr",
  "ATP-binding cassette (ABC) antibiotic efflux pump"                       = "ABC efflux",
  "major facilitator superfamily (MFS) antibiotic efflux pump"              = "MFS efflux",
  "multidrug and toxic compound extrusion (MATE) transporter"               = "MATE efflux",
  "resistance-nodulation-cell division (RND) antibiotic efflux pump"        = "RND efflux",
  "small multidrug resistance (SMR) antibiotic efflux pump"                 = "SMR efflux",
  "General Bacterial Porin with reduced permeability to beta-lactams"       = "Porin (beta-lactam)",
  "General Bacterial Porin with reduced permeability to peptide antibiotics"= "Porin (peptide)",
  "Outer Membrane Porin (Opr)"                                              = "Opr",
  "intrinsic colistin resistant phosphoethanolamine transferase"            = "Intrinsic PET",
  "MCR phosphoethanolamine transferase"                                     = "MCR",
  "pmr phosphoethanolamine transferase"                                     = "Pmr",
  "glycopeptide resistance gene cluster"                                    = "Van cluster",
  "Van ligase" = "Van", "vanH" = "Van", "vanT" = "Van",
  "vanU" = "Van", "vanW" = "Van", "vanXY" = "Van", "vanY" = "Van",
  "defensin resistant mprF"                                                 = "MprF",
  "Intrinsic peptide antibiotic resistant Lps"                              = "Lps",
  "rifampin ADP-ribosyltransferase (Arr)"                                   = "Arr",
  "rifamycin-resistant beta-subunit of RNA polymerase (rpoB)"               = "RpoB",
  "Bleomycin resistant protein"                                             = "Ble",
  "fosfomycin thiol transferase"                                            = "Fos",
  "kdpDE"                                                                   = "KdpDE",
  "methicillin resistant PBP2"                                              = "PBP2 (mecA-like)",
  "RbpA bacterial RNA polymerase-binding protein"                           = "RbpA",
  "streptothricin acetyltransferase (SAT)"                                  = "SAT",
  "tunicamycin resistance protein"                                          = "TmrB",
  "undecaprenyl pyrophosphate related proteins"                             = "UppP",
  "chloramphenicol acetyltransferase (CAT)"                                 = "CAT",
  "nitroimidazole reductase"                                                = "Nim"
)

map_family <- function(family, map = CARD_FAMILY_MAP) {
  f <- trimws(family)
  if (is.na(f) || f == "" || f == "NA") return("Unclassified")
  # Beta-lactamases: ends with "beta-lactamase"
  if (grepl("beta-lactamase$", f, ignore.case = TRUE)) {
    prefix <- trimws(sub("(?i)(-\\d+)?-like\\s+beta-lactamase$|\\s+beta-lactamase$",
                         "", f, perl = TRUE))
    if (tolower(prefix) == "ampc-type") return("AmpC")
    if (nchar(prefix) == 0) return("Other beta-lactamase")
    tokens <- strsplit(prefix, "\\s+")[[1]]
    if (length(tokens) > 2) {
      last <- tokens[length(tokens)]
      return(paste0(toupper(substr(last,1,1)), substr(last,2,nchar(last))))
    }
    return(prefix)
  }
  # Aminoglycoside-modifying enzymes: AAC / ANT / APH
  m <- regmatches(f, regexpr("^(AAC|ANT|APH)", f))
  if (length(m) > 0 && nchar(m) > 0) return(m)
  # Manual map
  if (f %in% names(map)) return(unname(map[f]))
  f   # unmapped — kept as-is, flagged below
}

if (col_family %in% colnames(rgi_all)) {
  rgi_split <- rgi_all %>%
    separate_rows(!!sym(col_family), sep = "; ") %>%
    mutate(amr_gene_family_raw = str_trim(.data[[col_family]]),
           family_prefix        = vapply(amr_gene_family_raw, map_family, character(1)))

  agg_col <- if (!is.null(col_reads) && col_reads %in% colnames(rgi_split))
    "mapped_reads" else "count"

  fam_prefix <- rgi_split %>%
    group_by(sample_id, family_prefix) %>%
    summarise(n_variants = n_distinct(.data[[col_aro]]),
              count      = n(),
              !!agg_col := if (agg_col == "mapped_reads")
                sum(.data[[col_reads]], na.rm = TRUE) else n(),
              .groups = "drop") %>%
    left_join(stats %>% select(sample_id, total_bases_gb), by = "sample_id") %>%
    mutate(normalized_abundance = .data[[agg_col]] / total_bases_gb)

  save_tsv(fam_prefix, output_family_prefix)

  # Flag unmapped CARD families for review (new CARD versions / new samples)
  unmapped <- rgi_split %>%
    filter(!grepl("beta-lactamase$", amr_gene_family_raw, ignore.case = TRUE),
           !grepl("^(AAC|ANT|APH)", amr_gene_family_raw),
           !amr_gene_family_raw %in% names(CARD_FAMILY_MAP),
           !is.na(amr_gene_family_raw),
           trimws(amr_gene_family_raw) != "") %>%
    pull(amr_gene_family_raw) %>% unique()

  if (length(unmapped) > 0) {
    cat("  [WARN] Unmapped CARD families (review recommended):\n")
    for (u in unmapped) cat("    -", u, "\n")
  } else {
    cat("  All CARD families successfully mapped\n")
  }

  # Top 15 families by normalized abundance
  cat("\n  Top 15 families (normalized abundance):\n")
  fam_prefix %>%
    group_by(family_prefix) %>% summarise(total = sum(normalized_abundance), .groups="drop") %>%
    arrange(desc(total)) %>% head(15) %>%
    { cat(sprintf("    %2d. %-25s %.2f\n", seq_len(nrow(.)), .$family_prefix, .$total)) }

} else {
  cat("  [WARN] amr_gene_family not found — empty file\n")
  write_tsv(tibble(), output_family_prefix)
}

# ── Summary ───────────────────────────────────────────────────────────────────
cat(sprintf(
  "\n══ COMPLETED [%s] ══\n  protein homolog model filter: %s\n  12 output files → results/r_analysis/%s/\n\n",
  tool,
  if (tool == "rgi_bwt") "yes" else "N/A",
  tool
))

sink(type = "message"); sink(); close(log)
