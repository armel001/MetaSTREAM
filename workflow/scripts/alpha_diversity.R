#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# alpha_diversity.R — Alpha diversity indices
# Pipeline MetaSTREAM | Bracken taxonomic abundance matrix
#
# Metrics (vegan):
#   Observed richness, Shannon H', Simpson 1-D, Inverse Simpson,
#   Pielou J'  → computed on rarefied counts (common depth, seed=42)
#   Chao1      → computed on full counts (bias-corrected estimator)
#
# References:
#   Oksanen et al. (2022) vegan
#   Willis (2019) Front. Microbiol.
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

# ── I/O ───────────────────────────────────────────────────────────────────────
if (exists("snakemake")) {
  input_matrix <- snakemake@input$matrix
  output_file  <- snakemake@output$diversity
  log_file     <- snakemake@log[[1]]
  level        <- snakemake@params$level %||% "Taxa"
  sink(log_file, append = FALSE, split = TRUE)
} else {
  args         <- commandArgs(trailingOnly = TRUE)
  input_matrix <- args[1]
  output_file  <- args[2]
  level        <- if (length(args) >= 3) args[3] else "Taxa"
}

cat(sprintf(
  "══ ALPHA DIVERSITY — MetaSTREAM ══\n  Level : %s\n  Input : %s\n  Output: %s\n  Date  : %s\n\n",
  level, input_matrix, output_file, format(Sys.time(), "%Y-%m-%d %H:%M")
))

# ── Load matrix ───────────────────────────────────────────────────────────────
raw         <- read_tsv(input_matrix, show_col_types = FALSE)
sample_cols <- setdiff(names(raw), c("name", "taxonomy_id"))

cat(sprintf("  %s: %d | Samples: %d (%s)\n",
            level, nrow(raw), length(sample_cols),
            paste(sample_cols, collapse = ", ")))

# ── Build count matrix (taxa × samples) ───────────────────────────────────────
otu_mat           <- as.matrix(raw[sample_cols])
rownames(otu_mat) <- raw$name
otu_t             <- t(otu_mat)   # vegan expects samples × taxa

# Integer coercion — required by estimateR(); Bracken may produce non-integers
otu_t_int <- round(otu_t)
n_diff    <- sum(abs(otu_t - otu_t_int) > 1e-6)
if (n_diff > 0)
  cat(sprintf("  [WARNING] %d non-integer values detected (Bracken) — rounded.\n", n_diff))

# ── Rarefaction to common depth ───────────────────────────────────────────────
depths    <- rowSums(otu_t_int)
min_depth <- min(depths)
cat(sprintf("  Depth : min=%d  max=%d  ratio=%.2f\n",
            min_depth, max(depths), max(depths) / min_depth))

set.seed(42)
otu_rare <- vegan::rrarefy(otu_t_int, sample = min_depth)
cat(sprintf("  Rarefied to %d reads (seed=42)\n", min_depth))

# ── Diversity indices ─────────────────────────────────────────────────────────
shannon     <- vegan::diversity(otu_rare, index = "shannon")
simpson     <- vegan::diversity(otu_rare, index = "simpson")
inv_simpson <- vegan::diversity(otu_rare, index = "invsimpson")
richness    <- vegan::specnumber(otu_rare)
pielou      <- ifelse(richness > 1, shannon / log(richness), 0)
chao1       <- vegan::estimateR(otu_t_int)["S.chao1", ]   # full counts — not rarefied

# ── Print sample-level summary ────────────────────────────────────────────────
cat(sprintf("\n  %-12s %6s %7s %7s %7s %7s\n", "Sample","S","H'","1-D","J'","Chao1"))
cat(strrep("─", 56), "\n")
for (s in rownames(otu_t))
  cat(sprintf("  %-12s %6d %7.3f %7.3f %7.3f %7.1f\n",
              s, richness[s], shannon[s], simpson[s], pielou[s], chao1[s]))

# ── Export ────────────────────────────────────────────────────────────────────
results <- tibble(
  sample            = names(richness),
  observed_richness = as.integer(richness),
  shannon_diversity = round(shannon, 4),
  simpson_diversity = round(simpson, 4),
  inverse_simpson   = round(inv_simpson, 4),
  pielou_evenness   = round(pielou, 4),
  chao1_richness    = round(chao1, 2),
  total_reads       = as.integer(depths)
) %>% arrange(sample)

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write_tsv(results, output_file)
cat(sprintf("\n  [OK] %s  (%d rows)\n", output_file, nrow(results)))

# ── Statistical summary ───────────────────────────────────────────────────────
cat("\n══ SUMMARY ══\n")
for (m in c("observed_richness","shannon_diversity","simpson_diversity",
            "pielou_evenness","chao1_richness")) {
  v <- results[[m]]
  cat(sprintf("  %-22s  mean=%.3f  sd=%.3f  [%.3f – %.3f]\n",
              m, mean(v), sd(v), min(v), max(v)))
}

avg_H <- mean(results$shannon_diversity)
cat(sprintf("\n  Mean Shannon = %.3f → %s diversity\n", avg_H,
            case_when(avg_H < 2 ~ "low", avg_H < 3 ~ "moderate",
                      avg_H < 4 ~ "high", TRUE ~ "very high")))

avg_J <- mean(results$pielou_evenness)
cat(sprintf("  Mean Pielou  = %.3f → %s evenness\n", avg_J,
            case_when(avg_J > 0.7 ~ "high (uniform distribution)",
                      avg_J > 0.5 ~ "moderate (few dominant taxa)",
                      TRUE        ~ "low (strong dominance)")))

cat("\n══ ALPHA DIVERSITY COMPLETED ══\n\n")
