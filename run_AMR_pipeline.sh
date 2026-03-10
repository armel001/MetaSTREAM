#!/usr/bin/env bash

################################################################################
# AMR Pipeline for Metagenomic Analysis
# Author       : Thibaut Armel Chérif GNIMADI
# Affiliation  : CERFIG
# Description  : Metagenomic analysis pipeline with Kraken2/Bracken and RGI
# Version      : 4.0
# Date         : 2026-03-09
################################################################################

set -euo pipefail

CORES=22
CONFIG_FILE="config/config.yaml"
SNAKEFILE="workflow/Snakefile"
RESULTS_DIR="results"
SUMMARY_DIR="${RESULTS_DIR}/summary"
TARGET="all"

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Options:
  -t, --target TARGET   Snakemake target rule to execute (default: all)
  -c, --cores  N        Number of cores to use (default: ${CORES})
  -n, --dry-run         Dry run — show what would be executed without running
  -h, --help            Show this help message

Available targets:

  ── Taxonomy ──────────────────────────────────────────────────────────────────
  quality_control     QC (NanoPlot) + stats + clean reads
  taxonomy_kraken     + Kraken2 + Bracken + matrices d'abondance
  taxonomy_analysis   + filtrage, normalisation, diversité alpha, summaries
  taxonomy_viz        + rapport HTML taxonomique
  taxonomy_all        Pipeline taxonomique complet (= all sans RGI)

  ── Pipeline global ───────────────────────────────────────────────────────────
  all                 Tout : taxonomie + RGI main + RGI BWT + MGE

  ── RGI — analyse R (matrices) ────────────────────────────────────────────────
  r_analysis_rgi      Matrices R — rgi main (assembly)
  r_analysis_bwt      Matrices R — rgi_bwt (reads)
  r_analysis_all      Matrices R — les deux outils

  ── RGI — visualisation ───────────────────────────────────────────────────────
  viz_rgi             Rapport HTML — rgi main
  viz_bwt             Rapport HTML — rgi_bwt
  viz_all             Rapports HTML — les deux outils

  ── RGI — pipelines complets ──────────────────────────────────────────────────
  rgi_main_all        RGI assembly : détection → matrices R → rapport HTML
  rgi_bwt_all         RGI reads    : détection → matrices R → rapport HTML
  rgi_compare_all     RGI main + BWT + rapport comparatif

  ── MGE ───────────────────────────────────────────────────────────────────────
  mge_all             Détection des éléments génétiques mobiles (mobileOG-db)

Examples:
  $(basename "$0")                           # Pipeline complet (all)
  $(basename "$0") -t taxonomy_all           # Taxonomie uniquement
  $(basename "$0") -t taxonomy_viz           # Juste le rapport HTML taxo
  $(basename "$0") -t r_analysis_bwt         # Matrices R BWT uniquement
  $(basename "$0") -t viz_bwt                # Rapport HTML BWT uniquement
  $(basename "$0") -t rgi_bwt_all            # Pipeline BWT complet
  $(basename "$0") -t rgi_compare_all        # RGI main + BWT + comparaison
  $(basename "$0") -t mge_all                # MGE uniquement
  $(basename "$0") -t rgi_main_all -n        # Dry run RGI main
  $(basename "$0") -t all -c 32              # Pipeline complet sur 32 cœurs
EOF
    exit 0
}

DRY_RUN=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -t|--target) TARGET="$2"; shift 2 ;;
        -c|--cores)  CORES="$2";  shift 2 ;;
        -n|--dry-run) DRY_RUN="--dry-run"; shift ;;
        -h|--help) usage ;;
        *) echo "ERROR: Unknown option: $1"; usage ;;
    esac
done

VALID_TARGETS=(
    # Taxonomy
    "quality_control"
    "taxonomy_kraken"
    "taxonomy_analysis"
    "taxonomy_viz"
    "taxonomy_all"
    # Global
    "all"
    # R analysis
    "r_analysis_rgi"
    "r_analysis_bwt"
    "r_analysis_all"
    # Visualisation
    "viz_rgi"
    "viz_bwt"
    "viz_all"
    # RGI pipelines
    "rgi_main_all"
    "rgi_bwt_all"
    "rgi_compare_all"
    # MGE
    "mge_all"
)
VALID=false
for t in "${VALID_TARGETS[@]}"; do
    [[ "$TARGET" == "$t" ]] && VALID=true && break
done

if [[ "$VALID" == false ]]; then
    echo "ERROR: Invalid target '${TARGET}'"
    echo "Run '$(basename "$0") --help' for the full list of available targets."
    exit 1
fi

if [[ ! -f "${CONFIG_FILE}" ]]; then echo "ERROR: Config file not found: ${CONFIG_FILE}"; exit 1; fi
if [[ ! -f "${SNAKEFILE}" ]];   then echo "ERROR: Snakefile not found: ${SNAKEFILE}";     exit 1; fi

echo "======================================"
echo "  AMR PIPELINE — EXECUTION"
echo "======================================"
echo "  Config    : ${CONFIG_FILE}"
echo "  Snakefile : ${SNAKEFILE}"
echo "  Target    : ${TARGET}"
echo "  Cores     : ${CORES}"
[[ -n "${DRY_RUN}" ]] && echo "  Mode      : DRY RUN"
echo "======================================"
echo ""

snakemake \
    --snakefile "${SNAKEFILE}" \
    --configfile "${CONFIG_FILE}" \
    --cores "${CORES}" \
    -j 2 \
    --use-conda \
    --rerun-triggers mtime \
    ${DRY_RUN} \
    -- "${TARGET}"

echo ""
echo "======================================"
[[ -n "${DRY_RUN}" ]] \
    && echo "  Dry run completed — no files written" \
    || echo "  Pipeline completed successfully"
echo "======================================"
echo ""
