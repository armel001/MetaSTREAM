#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# alpha_diversity.R — Indices de diversité alpha
# Pipeline MetaSTREAM | Bracken taxonomic abundance matrix (niveau paramétrable)
#
# Métriques calculées via vegan :
#   - Observed richness (S)   — calculée après raréfaction à profondeur commune
#   - Shannon H'        — vegan::diversity, base log naturel — après raréfaction
#   - Simpson 1-D       — vegan::diversity — après raréfaction
#   - Simpson inverse   — vegan::diversity — après raréfaction
#   - Pielou J'         — H' / ln(S) — après raréfaction
#   - Chao1             — vegan::estimateR, bias-corrected — sur comptages COMPLETS
#                          (l'estimateur corrige l'échantillonnage incomplet sans
#                          nécessiter de sous-échantillonnage ; le raréfier serait
#                          contre-productif)
#
# Références :
#   Oksanen et al. (2022) — vegan
#   Willis (2019) Front. Microbiol. — Rarefaction, Alpha Diversity, and Statistics
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

# ── Snakemake I/O ─────────────────────────────────────────────────────────────

if (exists("snakemake")) {
  input_matrix <- snakemake@input$matrix
  output_file  <- snakemake@output$diversity
  log_file     <- snakemake@log[[1]]
  level        <- if (!is.null(snakemake@params$level)) snakemake@params$level else "Taxa"
  sink(log_file, append = FALSE, split = TRUE)
} else {
  # Mode standalone
  args         <- commandArgs(trailingOnly = TRUE)
  input_matrix <- args[1]
  output_file  <- args[2]
  level        <- if (length(args) >= 3) args[3] else "Taxa"
}

cat("══════════════════════════════════════════════════════════════════════\n")
cat("  ALPHA DIVERSITY — MetaSTREAM\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Niveau : %s\n", level))
cat(sprintf("  Input  : %s\n", input_matrix))
cat(sprintf("  Output : %s\n", output_file))
cat(sprintf("  Date   : %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M")))

# ── Chargement de la matrice ──────────────────────────────────────────────────

cat("[Step 1] Chargement de la matrice d'abondance...\n")

raw <- read_tsv(input_matrix, show_col_types = FALSE)

# Extraire les colonnes samples (tout sauf name et taxonomy_id)
meta_cols   <- c("name", "taxonomy_id")
sample_cols <- setdiff(names(raw), meta_cols)

cat(sprintf("  %-9s: %d\n", level, nrow(raw)))
cat(sprintf("  Samples  : %d (%s)\n", length(sample_cols),
            paste(sample_cols, collapse = ", ")))

# ── Construction de la matrice de comptage ───────────────────────────────────

cat("\n[Step 2] Construction de la matrice taxa × samples...\n")

otu_mat <- raw %>%
  select(all_of(sample_cols)) %>%
  as.matrix()
rownames(otu_mat) <- raw$name

cat(sprintf("  Matrice  : %d taxa × %d samples\n", nrow(otu_mat), ncol(otu_mat)))

# ── Calcul des indices via vegan ──────────────────────────────────────────────
# vegan attend : lignes = samples, colonnes = taxa

cat("\n[Step 3] Calcul des indices de diversité alpha (vegan)...\n")

# Matrice transposée pour vegan : samples × taxa
otu_t <- t(otu_mat)

# Coercition entière — requise par estimateR() ; Bracken peut produire des
# comptages non entiers (réassignation probabiliste des reads multi-mappés)
otu_t_int <- round(otu_t)
n_diff <- sum(abs(otu_t - otu_t_int) > 1e-6)
if (n_diff > 0) {
  cat(sprintf("  [AVERTISSEMENT] %d valeurs non entières détectées (Bracken) — arrondies.\n", n_diff))
}

# Diagnostic profondeur de séquençage
depths <- rowSums(otu_t_int)
cat(sprintf("  Profondeur   : min=%d  max=%d  ratio max/min=%.2f\n",
            min(depths), max(depths), max(depths) / min(depths)))

# Raréfaction à profondeur commune pour Shannon / Simpson / Richesse / Pielou
# (Chao1 reste calculé sur les comptages complets — voir en-tête)
set.seed(42)
min_depth <- min(depths)
otu_rare  <- vegan::rrarefy(otu_t_int, sample = min_depth)
cat(sprintf("  Raréfaction  : %d reads (profondeur minimale observée, seed=42)\n", min_depth))

# Shannon H' — base log naturel
shannon     <- vegan::diversity(otu_rare, index = "shannon")

# Simpson 1-D
simpson     <- vegan::diversity(otu_rare, index = "simpson")

# Simpson inverse (1/D)
inv_simpson <- vegan::diversity(otu_rare, index = "invsimpson")

# Richesse observée S (post-raréfaction)
richness    <- vegan::specnumber(otu_rare)

# Pielou J' = H' / ln(S)
pielou      <- ifelse(richness > 1,
                      shannon / log(richness),
                      0)

# Chao1 — estimateur bias-corrected, sur comptages COMPLETS (non raréfiés)
chao1_res   <- vegan::estimateR(otu_t_int)
chao1       <- chao1_res["S.chao1", ]

# Total reads par sample — profondeur réelle (avant raréfaction), conservée
# pour traçabilité
total_reads <- depths

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
