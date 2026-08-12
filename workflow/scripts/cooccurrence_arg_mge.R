#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# cooccurrence_arg_mge.R — Co-occurrence ARG × MGE au niveau contig
# Pipeline MetaSTREAM
#
# Stratégie :
#   1. Charger RGI main (assembly) — un ARG par ligne, clé : Contig
#   2. Charger mobileOG hits — un MGE par ligne, clé : Specific Contig
#   3. Jointure inner_join sur contig_id × sample_id
#   4. Outputs :
#      01_cooccurrence_pairs.tsv    — toutes les paires ARG-MGE par contig
#      02_cooccurrence_summary.tsv  — comptages par (drug_class × mge_category)
#      03_cooccurrence_network.tsv  — format edge list pour visualisation
#      04_cooccurrence_stats.tsv    — statistiques par sample
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(tidyverse)
})

# ── Snakemake I/O ─────────────────────────────────────────────────────────────

if (exists("snakemake")) {
  rgi_files    <- snakemake@input$rgi_files
  mge_files    <- snakemake@input$mge_files
  out_pairs    <- snakemake@output$pairs
  out_summary  <- snakemake@output$summary
  out_network  <- snakemake@output$network
  out_stats    <- snakemake@output$stats
  log_file     <- snakemake@log[[1]]
  sink(log_file, append = FALSE, split = TRUE)
} else {
  # Mode standalone — adapter les chemins
  stop("Run via Snakemake or set paths manually in standalone mode")
}

cat("══════════════════════════════════════════════════════════════════════\n")
cat("  CO-OCCURRENCE ARG × MGE — MetaSTREAM\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat(sprintf("  RGI files : %d\n", length(rgi_files)))
cat(sprintf("  MGE files : %d\n", length(mge_files)))
cat(sprintf("  Date      : %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M")))

# ── Step 1 : Chargement RGI main ──────────────────────────────────────────────

cat("[Step 1] Chargement des hits RGI main (assembly)...\n")

rgi_all <- map_dfr(rgi_files, function(f) {
  sample <- basename(dirname(dirname(f)))   # results/HSE_X/rgi/HSE_X.txt
  df <- read_tsv(f, show_col_types = FALSE) %>%
    select(
      orf_id        = ORF_ID,
      contig_id     = Contig,
      cut_off       = Cut_Off,
      best_hit_aro  = Best_Hit_ARO,
      aro           = ARO,
      drug_class    = `Drug Class`,
      mechanism     = `Resistance Mechanism`,
      amr_family    = `AMR Gene Family`,
      pct_ref       = `Percentage Length of Reference Sequence`
    ) %>%
    mutate(sample_id = sample)
  cat(sprintf("  %-10s : %d ARG hits\n", sample, nrow(df)))
  df
})

cat(sprintf("\n  Total RGI hits : %d\n", nrow(rgi_all)))
cat(sprintf("  ARGs uniques   : %d\n", n_distinct(rgi_all$best_hit_aro)))
cat(sprintf("  Contigs ARG    : %d\n", n_distinct(paste(rgi_all$sample_id, rgi_all$contig_id))))

# ── Step 2 : Chargement mobileOG hits ────────────────────────────────────────

cat("\n[Step 2] Chargement des hits mobileOG...\n")

# Colonnes mobileOG — position fixe (header avec espaces)
mge_cols <- c(
  "seq_title", "query_title", "pident", "bitscore",
  "subj_len", "evalue", "query_len",
  "start_aln_subject", "end_aln_subject",
  "start_aln_query", "end_aln_query",
  "mobileog_id", "gene_name", "accession",
  "major_category", "minor_category", "source_db",
  "evidence_type", "contig_orf_name", "orf_start_stop_strand",
  "orf_start", "orf_end", "strand", "prodigal_id",
  "prodigal_contigs", "unique_orf", "partial_tag",
  "start_codon", "rbs_motif", "rbs_spacer", "gc_content",
  "specific_contig", "final_sample"
)

mge_all <- map_dfr(mge_files, function(f) {
  sample <- basename(dirname(dirname(f)))
  df <- read_csv(f, show_col_types = FALSE) %>%
    rename(
      contig_id      = `Specific Contig`,
      orf_name       = `Contig/ORF Name`,
      mobileog_id    = `mobileOG ID`,
      gene_name      = `Gene Name`,
      major_category = `Major mobileOG Category`,
      minor_category = `Minor mobileOG Category`,
      source_db      = `Source Database`,
      pident         = Pident,
      orf_start      = ORF_Start,
      orf_end        = ORF_End
    ) %>%
    select(contig_id, orf_name, mobileog_id, gene_name,
           major_category, minor_category, source_db,
           pident, orf_start, orf_end) %>%
    mutate(
      sample_id = sample,
      contig_id = str_trim(contig_id)
    )
  cat(sprintf("  %-10s : %d MGE hits\n", sample, nrow(df)))
  df
})

cat(sprintf("\n  Total MGE hits : %d\n", nrow(mge_all)))
cat(sprintf("  MGE uniques    : %d\n", n_distinct(mge_all$mobileog_id)))
cat(sprintf("  Contigs MGE    : %d\n", n_distinct(paste(mge_all$sample_id, mge_all$contig_id))))

# ── Step 3 : Jointure ARG × MGE sur contig_id ────────────────────────────────

cat("\n[Step 3] Jointure ARG × MGE sur contig_id × sample_id...\n")

# Contigs portant au moins 1 ARG
rgi_contigs <- rgi_all %>%
  select(sample_id, contig_id, orf_id, best_hit_aro,
         drug_class, mechanism, amr_family, cut_off, pct_ref) %>%
  distinct()

# Contigs portant au moins 1 MGE
mge_contigs <- mge_all %>%
  select(sample_id, contig_id, orf_name, mobileog_id,
         gene_name, major_category, minor_category,
         source_db, pident, orf_start, orf_end) %>%
  distinct()

# Inner join — seuls les contigs portant ARG ET MGE
cooccurrence <- inner_join(
  rgi_contigs,
  mge_contigs,
  by      = c("sample_id", "contig_id"),
  relationship = "many-to-many"
)

cat(sprintf("  Paires ARG-MGE   : %d\n", nrow(cooccurrence)))
cat(sprintf("  Contigs co-occur : %d\n", n_distinct(paste(cooccurrence$sample_id, cooccurrence$contig_id))))
cat(sprintf("  ARGs impliqués   : %d\n", n_distinct(cooccurrence$best_hit_aro)))
cat(sprintf("  MGEs impliqués   : %d\n", n_distinct(cooccurrence$mobileog_id)))

if (nrow(cooccurrence) == 0) {
  warning("Aucune co-occurrence détectée — vérifier les identifiants de contigs")
}

# ── Step 4 : Outputs ──────────────────────────────────────────────────────────

cat("\n[Step 4] Génération des fichiers de sortie...\n")

dir.create(dirname(out_pairs), showWarnings = FALSE, recursive = TRUE)

# 01 — Paires complètes
cooccurrence %>%
  arrange(sample_id, contig_id, best_hit_aro, major_category) %>%
  write_tsv(out_pairs)
cat(sprintf("  [OK] %s\n", out_pairs))

# 02 — Résumé drug_class × mge_category
summary_df <- cooccurrence %>%
  group_by(sample_id, drug_class, major_category) %>%
  summarise(
    n_pairs       = n(),
    n_contigs     = n_distinct(contig_id),
    n_arg_unique  = n_distinct(best_hit_aro),
    n_mge_unique  = n_distinct(mobileog_id),
    .groups = "drop"
  ) %>%
  arrange(sample_id, desc(n_pairs))

summary_df %>% write_tsv(out_summary)
cat(sprintf("  [OK] %s\n", out_summary))

# 03 — Edge list pour réseau (ARG ↔ MGE category)
network_df <- cooccurrence %>%
  group_by(best_hit_aro, drug_class, major_category) %>%
  summarise(
    weight        = n(),
    n_samples     = n_distinct(sample_id),
    n_contigs     = n_distinct(paste(sample_id, contig_id)),
    .groups = "drop"
  ) %>%
  rename(
    source = best_hit_aro,
    target = major_category
  ) %>%
  arrange(desc(weight))

network_df %>% write_tsv(out_network)
cat(sprintf("  [OK] %s\n", out_network))

# 04 — Statistiques par sample
stats_df <- cooccurrence %>%
  group_by(sample_id) %>%
  summarise(
    n_cooccurrence_contigs = n_distinct(contig_id),
    n_arg_types            = n_distinct(best_hit_aro),
    n_mge_types            = n_distinct(mobileog_id),
    n_drug_classes         = n_distinct(drug_class),
    n_mge_categories       = n_distinct(major_category),
    top_drug_class         = names(sort(table(drug_class), decreasing = TRUE))[1],
    top_mge_category       = names(sort(table(major_category), decreasing = TRUE))[1],
    .groups = "drop"
  )

stats_df %>% write_tsv(out_stats)
cat(sprintf("  [OK] %s\n", out_stats))

# ── Résumé final ──────────────────────────────────────────────────────────────

cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("  RÉSUMÉ CO-OCCURRENCE\n")
cat("══════════════════════════════════════════════════════════════════════\n")
print(as.data.frame(stats_df))

cat("\n  Top associations drug_class × MGE category :\n")
cooccurrence %>%
  count(drug_class, major_category, sort = TRUE) %>%
  head(10) %>%
  print()

cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("  CO-OCCURRENCE COMPLETED\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
