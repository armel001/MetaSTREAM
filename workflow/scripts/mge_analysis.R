#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# mge_analysis.R — MGE analysis matrices
# Pipeline MetaSTREAM | mobileOG-db beatrix-1.6 | Nanopore shotgun metagenomics
#
# Outputs:
#   01  mge_hits_raw.tsv          Filtered hits + sample_id
#   07  mge_normalized.tsv        Copies/Gb per major category × sample
#   08  mge_normalized_matrix.tsv Wide matrix samples × major categories
#   12  mge_pident_summary.tsv    Pident stats per category × sample
#   14  mge_presence_absence.tsv  Presence/absence major categories × sample
#   15  mge_run_summary.log       Run report
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(tidyverse))

# ── I/O ───────────────────────────────────────────────────────────────────────
if (exists("snakemake")) {
  hits_files <- snakemake@input$hits_files
  stats_file <- snakemake@input$stats
  outdir     <- snakemake@params$outdir
  min_pident <- snakemake@params$min_pident
  min_cov    <- snakemake@params$min_cov
  log_file   <- snakemake@log[[1]]
  sink(log_file, append = FALSE, split = TRUE)
} else {
  suppressPackageStartupMessages(library(argparse))
  parser <- ArgumentParser(description = "MGE analysis — standalone mode")
  parser$add_argument("--hits_dir",   required = TRUE)
  parser$add_argument("--stats",      required = TRUE)
  parser$add_argument("--outdir",     required = TRUE)
  parser$add_argument("--min_pident", type = "double", default = 90)
  parser$add_argument("--min_cov",    type = "double", default = 90)
  args       <- parser$parse_args()
  stats_file <- args$stats
  outdir     <- args$outdir
  min_pident <- args$min_pident
  min_cov    <- args$min_cov
  hits_files <- list.files(args$hits_dir,
                           pattern    = "assembly\\.fasta\\.mobileOG\\.Alignment\\.Out\\.csv",
                           recursive  = TRUE, full.names = TRUE)
  if (length(hits_files) == 0)
    stop("No .mobileOG.Alignment.Out.csv files found in: ", args$hits_dir)
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
t_start <- Sys.time()

cat(sprintf(
  "══ MGE ANALYSIS — MetaSTREAM ══\n  Date   : %s\n  Samples: %d\n  Pident : %.0f%%  Cov: %.0f%%\n  Outdir : %s\n\n",
  format(Sys.time(), "%Y-%m-%d %H:%M"), length(hits_files), min_pident, min_cov, outdir
))

# ── Helper ────────────────────────────────────────────────────────────────────
write_tsv_out <- function(df, filename) {
  path <- file.path(outdir, filename)
  write_tsv(df, path)
  cat(sprintf("  [OK] %-45s (%d rows)\n", filename, nrow(df)))
  invisible(path)
}

# Extracts sample_id from path: .../results/<sample>/mge/...
extract_sample_id <- function(path) {
  parts   <- strsplit(normalizePath(path), .Platform$file.sep)[[1]]
  mge_idx <- which(parts == "mge")
  if (length(mge_idx) == 0) return(basename(dirname(path)))
  parts[mge_idx - 1]
}

# ── Step 1 : Load & concatenate hits ─────────────────────────────────────────
cat("[Step 1] Loading hit files...\n")

raw <- lapply(hits_files, function(f) {
  df <- read_csv(f, show_col_types = FALSE, name_repair = "minimal")
  if (names(df)[1] %in% c("", "...1")) df <- df[, -1]
  df$sample_id <- extract_sample_id(f)
  df
}) %>% bind_rows()

cat(sprintf("  Total hits before filtering: %d\n", nrow(raw)))

# Fix duplicate/ambiguous column names (mobileOG-db header bug)
names(raw)[9]  <- "End_Alignment_Subject"
names(raw)[10] <- "Start_Alignment_Query"
names(raw)[11] <- "End_Alignment_Query"
names(raw)     <- names(raw) %>%
  str_replace_all("\\s+", "_") %>%
  str_replace_all("[/\\.]", "_")

raw <- raw %>% rename(
  major_category = Major_mobileOG_Category,
  minor_category = Minor_mobileOG_Category,
  source_db      = Source_Database,
  evidence_type  = Evidence_Type,
  contig         = Specific_Contig,
  orf_name       = Contig_ORF_Name,
  unique_orf     = Unique_ORF,
  gc_content     = GC_Content,
  gene_name      = Gene_Name,
  mobileog_id    = mobileOG_ID
)

# ── Step 2 : Quality filtering ────────────────────────────────────────────────
cat("\n[Step 2] Quality filtering (Pident + query coverage)...\n")

raw <- raw %>%
  mutate(
    align_length = End_Alignment_Query - Start_Alignment_Query + 1,
    query_cov    = align_length / Query_Sequence_Length * 100
  )

n_before <- nrow(raw)
hits     <- raw %>% filter(Pident >= min_pident, query_cov >= min_cov) %>%
  mutate(orf_global_id = paste(sample_id, orf_name, sep = "|"))
n_after  <- nrow(hits)

cat(sprintf("  Before: %d | After: %d (%.1f%% retained)\n",
            n_before, n_after, n_after / n_before * 100))

# Output 01 — filtered raw hits
write_tsv_out(
  hits %>% select(sample_id, contig, orf_name, orf_global_id, mobileog_id,
                  gene_name, major_category, minor_category, source_db,
                  evidence_type, Pident, Bitscore, query_cov, gc_content,
                  ORF_Start, ORF_End, Sense_or_Antisense_Strand, Partial_Tag) %>%
    arrange(sample_id, contig, orf_name),
  "01_mge_hits_raw.tsv"
)

# ── Step 3 : Counts + normalization ──────────────────────────────────────────
cat("\n[Step 3] Counts and normalization per Gb...\n")

stats <- read_tsv(stats_file, show_col_types = FALSE) %>%
  mutate(total_bases_gb = total_bases / 1e9)

missing <- setdiff(unique(hits$sample_id), stats$sample_id)
if (length(missing) > 0)
  warning("Samples missing from sequencing_stats: ", paste(missing, collapse = ", "))

major_counts <- hits %>%
  group_by(sample_id, major_category) %>%
  summarise(total_hits  = n(),
            unique_orfs = n_distinct(orf_global_id),
            n_contigs   = n_distinct(contig), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(relative_abundance = round(total_hits / sum(total_hits) * 100, 4)) %>%
  ungroup()

# Output 07 — normalized copies/Gb
mge_norm <- major_counts %>%
  left_join(stats %>% select(sample_id, total_bases_gb), by = "sample_id") %>%
  mutate(normalized_copies = round(total_hits / total_bases_gb, 4)) %>%
  select(sample_id, major_category, total_hits, unique_orfs, n_contigs,
         relative_abundance, total_bases_gb, normalized_copies)
write_tsv_out(mge_norm, "07_mge_normalized.tsv")

# Output 08 — wide normalized matrix
write_tsv_out(
  mge_norm %>% select(sample_id, major_category, normalized_copies) %>%
    pivot_wider(names_from = major_category, values_from = normalized_copies,
                values_fill = 0),
  "08_mge_normalized_matrix.tsv"
)

# ── Step 4 : Pident stats ────────────────────────────────────────────────────
cat("\n[Step 4] Pident statistics per category...\n")

# pct_100: proportion of hits with perfect identity (100%)
# Indicator of recently disseminated, well-documented MGEs
write_tsv_out(
  hits %>%
    group_by(sample_id, major_category) %>%
    summarise(n_hits        = n(),
              mean_pident   = round(mean(Pident,   na.rm = TRUE), 2),
              sd_pident     = round(sd(Pident,     na.rm = TRUE), 2),
              median_pident = round(median(Pident, na.rm = TRUE), 2),
              min_pident    = round(min(Pident,    na.rm = TRUE), 2),
              max_pident    = round(max(Pident,    na.rm = TRUE), 2),
              pct_100       = round(mean(Pident == 100) * 100, 2),
              .groups = "drop"),
  "12_mge_pident_summary.tsv"
)

# ── Step 5 : Presence / absence matrix ───────────────────────────────────────
cat("\n[Step 5] Presence/absence matrix...\n")

write_tsv_out(
  major_counts %>%
    select(sample_id, major_category, total_hits) %>%
    pivot_wider(names_from = major_category, values_from = total_hits,
                values_fill = 0) %>%
    mutate(across(-sample_id, ~ if_else(. > 0, 1L, 0L))),
  "14_mge_presence_absence.tsv"
)

# ── Step 6 : Run summary ─────────────────────────────────────────────────────
cat("\n[Step 6] Writing run summary...\n")

top_major   <- hits %>% count(major_category, name = "total_hits") %>% arrange(desc(total_hits))
norm_summary <- mge_norm %>%
  group_by(major_category) %>%
  summarise(mean_norm = round(mean(normalized_copies), 1), .groups = "drop") %>%
  arrange(desc(mean_norm))

log_lines <- c(
  "══ MGE ANALYSIS — RUN SUMMARY ══",
  sprintf("  Date          : %s  |  Elapsed: %s s",
          format(Sys.time(), "%Y-%m-%d %H:%M"),
          round(difftime(Sys.time(), t_start, units = "secs"), 1)),
  sprintf("  Pident ≥ %.0f%%  |  Query cov ≥ %.0f%%", min_pident, min_cov),
  sprintf("  Hits before/after filter: %d / %d (%.1f%%)",
          n_before, n_after, n_after / n_before * 100),
  sprintf("  Samples    : %d  (%s)",
          n_distinct(hits$sample_id),
          paste(sort(unique(hits$sample_id)), collapse = ", ")),
  sprintf("  Unique ORFs: %d  |  Carrying contigs: %d  |  Total hits: %d",
          n_distinct(hits$orf_global_id),
          n_distinct(paste(hits$sample_id, hits$contig)),
          nrow(hits)),
  "",
  "── Major categories (all samples) ──",
  sprintf("  %-35s : %6d  (%5.1f%%)",
          top_major$major_category, top_major$total_hits,
          top_major$total_hits / nrow(hits) * 100),
  "",
  "── Normalized abundance (mean copies/Gb) ──",
  sprintf("  %-35s : %8.1f copies/Gb",
          norm_summary$major_category, norm_summary$mean_norm),
  "",
  "── Output files ──",
  sprintf("  %s", list.files(outdir, pattern = "\\.(tsv|log)$")),
  "══════════════════════════════════════════"
)

writeLines(log_lines, file.path(outdir, "15_mge_run_summary.log"))
cat(paste(log_lines, collapse = "\n"), "\n")
cat("\n══ MGE ANALYSIS COMPLETED ══\n\n")
