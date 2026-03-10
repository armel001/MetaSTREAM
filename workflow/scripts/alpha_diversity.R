#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# alpha_diversity.R — Indices de diversité alpha
# Pipeline MetagenAMR | Bracken species-level abundance matrix
#
# Métriques calculées via phyloseq + vegan :
#   - Observed richness (S)
#   - Shannon H'        (vegan::diversity, base log naturel)
#   - Simpson 1-D       (vegan::diversity)
#   - Simpson inverse   (vegan::diversity)
#   - Pielou J'         (H' / ln(S))
#   - Chao1             (vegan::estimateR, bias-corrected)
#
# Références :
#   McMurdie & Holmes (2013) PLoS ONE — phyloseq
#   Oksanen et al. (2022) — vegan
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(vegan)
})

# ── Snakemake I/O ─────────────────────────────────────────────────────────────

if (exists("snakemake")) {
  input_matrix <- snakemake@input$matrix
  output_file  <- snakemake@output$diversity
  log_file     <- snakemake@log[[1]]
  sink(log_file, append = FALSE, split = TRUE)
} else {
  # Mode standalone
  args         <- commandArgs(trailingOnly = TRUE)
  input_matrix <- args[1]
  output_file  <- args[2]
}

cat("══════════════════════════════════════════════════════════════════════\n")
cat("  ALPHA DIVERSITY — MetagenAMR\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Input  : %s\n", input_matrix))
cat(sprintf("  Output : %s\n", output_file))
cat(sprintf("  Date   : %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M")))

# ── Chargement de la matrice ──────────────────────────────────────────────────

cat("[Step 1] Chargement de la matrice d'abondance...\n")

raw <- read_tsv(input_matrix, show_col_types = FALSE)

# Extraire les colonnes samples (tout sauf name et taxonomy_id)
meta_cols   <- c("name", "taxonomy_id")
sample_cols <- setdiff(names(raw), meta_cols)

cat(sprintf("  Espèces  : %d\n", nrow(raw)))
cat(sprintf("  Samples  : %d (%s)\n", length(sample_cols),
            paste(sample_cols, collapse = ", ")))

# ── Construction de l'objet phyloseq ─────────────────────────────────────────

cat("\n[Step 2] Construction de l'objet phyloseq...\n")

# OTU table : lignes = taxa, colonnes = samples
otu_mat <- raw %>%
  select(all_of(sample_cols)) %>%
  as.matrix()

rownames(otu_mat) <- raw$name

# Taxonomy table
tax_mat <- raw %>%
  select(name) %>%
  mutate(Species = name) %>%
  as.data.frame()
rownames(tax_mat) <- raw$name
tax_mat <- tax_matrix <- as.matrix(tax_mat)

# phyloseq object
OTU  <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX  <- tax_table(tax_mat)
ps   <- phyloseq(OTU, TAX)

cat(sprintf("  phyloseq : %d taxa × %d samples\n",
            ntaxa(ps), nsamples(ps)))

# ── Calcul des indices via vegan ──────────────────────────────────────────────
# vegan attend : lignes = samples, colonnes = taxa

cat("\n[Step 3] Calcul des indices de diversité alpha (vegan)...\n")

# Matrice transposée pour vegan : samples × taxa
otu_t <- t(otu_mat)

# Shannon H' — base log naturel (compatible phyloseq/vegan standard)
shannon     <- vegan::diversity(otu_t, index = "shannon")

# Simpson 1-D
simpson     <- vegan::diversity(otu_t, index = "simpson")

# Simpson inverse (1/D)
inv_simpson <- vegan::diversity(otu_t, index = "invsimpson")

# Richesse observée S
richness    <- vegan::specnumber(otu_t)

# Pielou J' = H' / ln(S)
pielou      <- ifelse(richness > 1,
                      shannon / log(richness),
                      0)

# Chao1 — estimateur bias-corrected via vegan::estimateR
chao1_res   <- vegan::estimateR(otu_t)
chao1       <- chao1_res["S.chao1", ]

# Total reads par sample
total_reads <- rowSums(otu_t)

cat(sprintf("  %-12s  %6s  %7s  %7s  %7s  %7s\n",
            "Sample", "S", "H'", "1-D", "J'", "Chao1"))
cat(paste(rep("─", 65), collapse = ""), "\n")
for (s in rownames(otu_t)) {
  cat(sprintf("  %-12s  %6d  %7.3f  %7.3f  %7.3f  %7.1f\n",
              s, richness[s], shannon[s], simpson[s],
              pielou[s], chao1[s]))
}

# ── Assemblage du tableau de sortie ──────────────────────────────────────────

cat("\n[Step 4] Export des résultats...\n")

results <- tibble(
  sample           = names(richness),
  observed_richness = as.integer(richness),
  shannon_diversity = round(shannon, 4),
  simpson_diversity = round(simpson, 4),
  inverse_simpson   = round(inv_simpson, 4),
  pielou_evenness   = round(pielou, 4),
  chao1_richness    = round(chao1, 2),
  total_reads       = as.integer(total_reads)
) %>%
  arrange(sample)

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write_tsv(results, output_file)

cat(sprintf("  [OK] %s  (%d lignes)\n", output_file, nrow(results)))

# ── Résumé statistique ────────────────────────────────────────────────────────

cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("  RÉSUMÉ STATISTIQUE\n")
cat("══════════════════════════════════════════════════════════════════════\n")

for (metric in c("observed_richness", "shannon_diversity",
                 "simpson_diversity", "pielou_evenness", "chao1_richness")) {
  vals <- results[[metric]]
  cat(sprintf("  %-22s  mean=%.3f  sd=%.3f  [%.3f – %.3f]\n",
              metric, mean(vals), sd(vals), min(vals), max(vals)))
}

# Interprétation Shannon
avg_H <- mean(results$shannon_diversity)
interp <- case_when(
  avg_H < 2 ~ "faible",
  avg_H < 3 ~ "modérée",
  avg_H < 4 ~ "élevée",
  TRUE      ~ "très élevée"
)
cat(sprintf("\n  Shannon moyen = %.3f → diversité %s\n", avg_H, interp))

avg_J <- mean(results$pielou_evenness)
interp_J <- case_when(
  avg_J > 0.7 ~ "élevée (distribution uniforme)",
  avg_J > 0.5 ~ "modérée (quelques espèces dominantes)",
  TRUE        ~ "faible (forte domination)"
)
cat(sprintf("  Pielou moyen  = %.3f → équitabilité %s\n", avg_J, interp_J))

cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("  ALPHA DIVERSITY COMPLETED\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
