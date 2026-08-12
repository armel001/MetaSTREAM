#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
# mge_analysis.R — Matrices d'analyse MGE
# Pipeline MetaSTREAM | mobileOG-db beatrix-1.6 | Nanopore shotgun metagenomics
# ══════════════════════════════════════════════════════════════════════════════
#
# Outputs :
#   01  mge_hits_raw.tsv          Hits filtrés + sample_id (base downstream)
#   07  mge_normalized.tsv        Copies / Gb par catégorie majeure × sample
#   08  mge_normalized_matrix.tsv Matrice wide samples × major categories
#   12  mge_pident_summary.tsv    Statistiques Pident par catégorie × sample
#   14  mge_presence_absence.tsv  Présence/absence major categories × sample
#   15  mge_run_summary.log       Rapport de synthèse
#
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(tidyverse))

# ── Mode Snakemake ou standalone ──────────────────────────────────────────────

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
  parser$add_argument("--hits_dir",   required = TRUE,
    help = "Répertoire results/ contenant les sous-dossiers sample/mge/")
  parser$add_argument("--stats",      required = TRUE,
    help = "sequencing_stats.tsv (colonnes : sample_id, total_bases)")
  parser$add_argument("--outdir",     required = TRUE,
    help = "Répertoire de sortie")
  parser$add_argument("--min_pident", type = "double", default = 90,
    help = "Seuil Pident minimum (défaut : 90)")
  parser$add_argument("--min_cov",    type = "double", default = 90,
    help = "Seuil couverture query minimum (défaut : 90)")
  args       <- parser$parse_args()
  stats_file <- args$stats
  outdir     <- args$outdir
  min_pident <- args$min_pident
  min_cov    <- args$min_cov

  hits_files <- list.files(
    path      = args$hits_dir,
    pattern   = "assembly\\.fasta\\.mobileOG\\.Alignment\\.Out\\.csv",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(hits_files) == 0)
    stop("Aucun fichier .mobileOG.Alignment.Out.csv trouvé dans : ", args$hits_dir)
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

t_start <- Sys.time()
cat("══════════════════════════════════════════════════════════════════════\n")
cat("  MGE ANALYSIS — MetaSTREAM\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat(sprintf("  Date       : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M")))
cat(sprintf("  Samples    : %d fichiers hits\n", length(hits_files)))
cat(sprintf("  Min Pident : %.0f%%\n", min_pident))
cat(sprintf("  Min Cov    : %.0f%%\n", min_cov))
cat(sprintf("  Outdir     : %s\n", outdir))
cat("══════════════════════════════════════════════════════════════════════\n\n")

# ── Helper ────────────────────────────────────────────────────────────────────

write_tsv_out <- function(df, filename) {
  path <- file.path(outdir, filename)
  write_tsv(df, path)
  cat(sprintf("  [OK] %-45s (%d lignes)\n", filename, nrow(df)))
  invisible(path)
}

# ── Step 1 : Chargement et concaténation ──────────────────────────────────────

cat("[Step 1] Chargement des fichiers hits...\n")

# Extraire sample_id depuis le chemin : .../results/<sample>/mge/...
extract_sample_id <- function(path) {
  parts   <- strsplit(normalizePath(path), .Platform$file.sep)[[1]]
  mge_idx <- which(parts == "mge")
  if (length(mge_idx) == 0) return(basename(dirname(path)))
  parts[mge_idx - 1]
}

raw_list <- lapply(hits_files, function(f) {
  sid <- extract_sample_id(f)
  df  <- read_csv(f, show_col_types = FALSE, name_repair = "minimal")
  # Supprimer colonne index automatique (sans nom ou "...1")
  if (names(df)[1] %in% c("", "...1")) df <- df[, -1]
  df$sample_id <- sid
  df
})

raw <- bind_rows(raw_list)
cat(sprintf("  Hits totaux avant filtrage : %d\n", nrow(raw)))

# Normaliser les noms de colonnes
# Renommer par position les colonnes au header ambigü (mobileOG-db bug)
# Col 9  = End of Alignment in Subject (dupliqué avec col 11 dans le header)
# Col 10 = Start of Alignment in Query
# Col 11 = End of Alignment in Query
names(raw)[9]  <- "End_Alignment_Subject"
names(raw)[10] <- "Start_Alignment_Query"
names(raw)[11] <- "End_Alignment_Query"

# Normalisation générale des noms restants
names(raw) <- names(raw) %>%
  str_replace_all("\\s+", "_") %>%
  str_replace_all("[/\\.]", "_")

# Renommer les colonnes clés
raw <- raw %>% rename(
  major_category = Major_mobileOG_Category,
  minor_category = Minor_mobileOG_Category,
  source_db      = Source_Database,
  evidence_type  = Evidence_Type,
  contig         = Specific_Contig,
  orf_name       = Contig_ORF_Name,   # identifiant ORF réel : contig_XXXX_N
  unique_orf     = Unique_ORF,        # index local mobileOG-db (non global)
  gc_content     = GC_Content,
  gene_name      = Gene_Name,
  mobileog_id    = mobileOG_ID
)

# ── Step 2 : Filtrage qualité ─────────────────────────────────────────────────

cat("\n[Step 2] Filtrage qualité (Pident + couverture query)...\n")

# query_cov = fraction de l'ORF (query) couverte par l'alignement
# = longueur alignée / longueur totale du query × 100
raw <- raw %>%
  mutate(
    align_length = End_Alignment_Query - Start_Alignment_Query + 1,
    query_cov    = align_length / Query_Sequence_Length * 100
  )

n_before <- nrow(raw)
hits <- raw %>% filter(Pident >= min_pident, query_cov >= min_cov)
n_after  <- nrow(hits)

cat(sprintf("  Avant filtrage : %d hits\n",  n_before))
cat(sprintf("  Après filtrage : %d hits  (retenus : %.1f%%)\n",
            n_after, n_after / n_before * 100))

# ── Sortie 01 : hits bruts filtrés ────────────────────────────────────────────
# Conserver toutes les colonnes utiles pour les analyses downstream
# (co-occurrence ARG × MGE notamment : contig, ORF_Start, ORF_End)

# Identifiant ORF global : sample_id + orf_name (contig_XXXX_N)
# orf_name est unique au sein d'un sample ; le préfixe sample_id le rend global
hits <- hits %>%
  mutate(orf_global_id = paste(sample_id, orf_name, sep = "|"))

hits_out <- hits %>%
  select(
    sample_id, contig, orf_name, orf_global_id, mobileog_id, gene_name,
    major_category, minor_category, source_db, evidence_type,
    Pident, Bitscore, query_cov, gc_content,
    ORF_Start, ORF_End, Sense_or_Antisense_Strand, Partial_Tag
  ) %>%
  arrange(sample_id, contig, orf_name)

write_tsv_out(hits_out, "01_mge_hits_raw.tsv")

# ── Step 3 : Comptages par catégorie majeure + normalisation ──────────────────

cat("\n[Step 3] Comptages et normalisation par Gb...\n")

stats <- read_tsv(stats_file, show_col_types = FALSE) %>%
  mutate(total_bases_gb = total_bases / 1e9)

missing_samples <- setdiff(unique(hits$sample_id), stats$sample_id)
if (length(missing_samples) > 0)
  warning("Samples absents de sequencing_stats : ",
          paste(missing_samples, collapse = ", "))

major_counts <- hits %>%
  group_by(sample_id, major_category) %>%
  summarise(
    total_hits  = n(),
    unique_orfs = n_distinct(orf_global_id),
    n_contigs   = n_distinct(contig),
    .groups     = "drop"
  ) %>%
  group_by(sample_id) %>%
  mutate(
    total_hits_sample  = sum(total_hits),
    relative_abundance = round(total_hits / total_hits_sample * 100, 4)
  ) %>%
  ungroup()

# ── Sortie 07 : normalisé par Gb ──────────────────────────────────────────────

mge_norm <- major_counts %>%
  left_join(stats %>% select(sample_id, total_bases_gb), by = "sample_id") %>%
  mutate(normalized_copies = round(total_hits / total_bases_gb, 4)) %>%
  select(sample_id, major_category, total_hits, unique_orfs, n_contigs,
         relative_abundance, total_bases_gb, normalized_copies)

write_tsv_out(mge_norm, "07_mge_normalized.tsv")

# ── Sortie 08 : matrice wide normalisée ───────────────────────────────────────

norm_matrix <- mge_norm %>%
  select(sample_id, major_category, normalized_copies) %>%
  pivot_wider(
    names_from  = major_category,
    values_from = normalized_copies,
    values_fill = 0
  )

write_tsv_out(norm_matrix, "08_mge_normalized_matrix.tsv")

# ── Step 4 : Statistiques Pident par catégorie ────────────────────────────────

cat("\n[Step 4] Statistiques Pident par catégorie...\n")

# pct_100 : proportion de hits à identité parfaite (100%)
# Signal d'éléments identiques à des références connues dans mobileOG-db
# → indicateur de dissémination récente et documentée

pident_summary <- hits %>%
  group_by(sample_id, major_category) %>%
  summarise(
    n_hits      = n(),
    mean_pident = round(mean(Pident,        na.rm = TRUE), 2),
    sd_pident   = round(sd(Pident,          na.rm = TRUE), 2),
    median_pident = round(median(Pident,    na.rm = TRUE), 2),
    min_pident  = round(min(Pident,         na.rm = TRUE), 2),
    max_pident  = round(max(Pident,         na.rm = TRUE), 2),
    pct_100     = round(mean(Pident == 100) * 100, 2),
    .groups     = "drop"
  )

write_tsv_out(pident_summary, "12_mge_pident_summary.tsv")

# ── Step 5 : Matrice présence / absence ───────────────────────────────────────

cat("\n[Step 5] Matrice présence/absence...\n")

# Base : matrice wide des hits absolus
major_matrix_wide <- major_counts %>%
  select(sample_id, major_category, total_hits) %>%
  pivot_wider(
    names_from  = major_category,
    values_from = total_hits,
    values_fill = 0
  )

presence <- major_matrix_wide %>%
  mutate(across(-sample_id, ~ if_else(. > 0, 1L, 0L)))

write_tsv_out(presence, "14_mge_presence_absence.tsv")

# ── Step 6 : Résumé ───────────────────────────────────────────────────────────

cat("\n[Step 6] Rédaction du résumé...\n")

t_end   <- Sys.time()
elapsed <- round(difftime(t_end, t_start, units = "secs"), 1)

# Top 5 catégories majeures (tous samples confondus)
top_major <- hits %>%
  count(major_category, name = "total_hits") %>%
  arrange(desc(total_hits))

# Statistiques globales
total_hits_all <- nrow(hits)
n_samples_ok   <- n_distinct(hits$sample_id)
n_orfs_total   <- n_distinct(hits$orf_global_id)
n_contigs_total <- n_distinct(paste(hits$sample_id, hits$contig))

log_lines <- c(
  "══════════════════════════════════════════════════════════════════════",
  "  MGE ANALYSIS — RUN SUMMARY",
  "══════════════════════════════════════════════════════════════════════",
  sprintf("  Date              : %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  sprintf("  Elapsed           : %s s", elapsed),
  "",
  "── Paramètres de filtrage ────────────────────────────────────────────",
  sprintf("  Min Pident        : %.0f%%", min_pident),
  sprintf("  Min query cov     : %.0f%%", min_cov),
  sprintf("  Hits avant filtre : %d", n_before),
  sprintf("  Hits retenus      : %d (%.1f%%)",
          n_after, n_after / n_before * 100),
  "",
  "── Données ───────────────────────────────────────────────────────────",
  sprintf("  Samples           : %d  (%s)",
          n_samples_ok,
          paste(sort(unique(hits$sample_id)), collapse = ", ")),
  sprintf("  ORFs uniques      : %d", n_orfs_total),
  sprintf("  Contigs porteurs  : %d  (unique sample × contig)", n_contigs_total),
  sprintf("  Hits totaux       : %d", total_hits_all),
  "",
  "── Catégories majeures — total hits (tous samples) ──────────────────",
  paste(
    sprintf("  %-35s : %6d  (%5.1f%%)",
            top_major$major_category,
            top_major$total_hits,
            top_major$total_hits / total_hits_all * 100),
    collapse = "\n"
  ),
  "",
  "── Normalisation (copies / Gb, moyenne inter-samples) ───────────────",
  paste(
    mge_norm %>%
      group_by(major_category) %>%
      summarise(mean_norm = round(mean(normalized_copies), 1), .groups = "drop") %>%
      arrange(desc(mean_norm)) %>%
      mutate(line = sprintf("  %-35s : %8.1f copies/Gb",
                            major_category, mean_norm)) %>%
      pull(line),
    collapse = "\n"
  ),
  "",
  "── Fichiers générés ──────────────────────────────────────────────────",
  paste(
    sprintf("  %s", list.files(outdir, pattern = "\\.(tsv|log)$")),
    collapse = "\n"
  ),
  "══════════════════════════════════════════════════════════════════════"
)

writeLines(log_lines, file.path(outdir, "15_mge_run_summary.log"))
cat(paste(log_lines, collapse = "\n"), "\n")

cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("  MGE ANALYSIS COMPLETED\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
