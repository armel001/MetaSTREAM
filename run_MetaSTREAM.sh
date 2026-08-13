#!/usr/bin/env bash
# AMR Pipeline for Metagenomic Analysis
# Author       : Thibaut Armel Chérif GNIMADI
# Affiliation  : CERFIG
# Description  : Metagenomic analysis pipeline with Kraken2/Bracken and RGI
# Version      : 4.2
# Date         : 2026-07-21
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
  taxonomy_kraken     + Kraken2 + Bracken + abundance matrices
  taxonomy_analysis   + filtering, normalization, alpha diversity, summaries
  taxonomy_viz        + taxonomic HTML report
  taxonomy_all        Complete taxonomy pipeline (= all without RGI)

  ── Medaka ─────────────────────────────────────────────────────────────────────
  medaka_all          Medaka polishing of all raw Flye assemblies

  ── PlasmidFinder ─────────────────────────────────────────────────────────────
  plasmidfinder_all   Plasmid replicon typing (Inc groups) on Medaka-polished
                       assemblies (depends on medaka_all)

  ── Full pipeline ─────────────────────────────────────────────────────────────
  all                 Everything: taxonomy + Medaka + PlasmidFinder + RGI main
                       + RGI BWT + MGE

  ── RGI — R analysis (matrices) ───────────────────────────────────────────────
  r_analysis_rgi      R matrices — RGI main (assembly)
  r_analysis_bwt      R matrices — RGI BWT (reads)
  r_analysis_all      R matrices — both tools

  ── RGI — visualization ───────────────────────────────────────────────────────
  viz_rgi             HTML report — RGI main
  viz_bwt             HTML report — RGI BWT
  viz_all             HTML reports — both tools

  ── RGI — full pipelines ──────────────────────────────────────────────────────
  rgi_main_all        RGI assembly: detection → R matrices → HTML report
  rgi_bwt_all         RGI reads: detection → R matrices → HTML report
  rgi_compare_all     RGI main + BWT + comparative report

  ── MGE ───────────────────────────────────────────────────────────────────────
  mge_all             Mobile genetic element detection (mobileOG-db)
  cooccurrence_all    ARG × MGE co-occurrence

Examples:
  $(basename "$0")                           # Complete pipeline (all)
  $(basename "$0") -t taxonomy_all           # Taxonomy only
  $(basename "$0") -t taxonomy_viz           # Taxonomy HTML report only
  $(basename "$0") -t medaka_all             # Medaka polishing only
  $(basename "$0") -t plasmidfinder_all      # Plasmid typing only
  $(basename "$0") -t r_analysis_bwt         # BWT R matrices only
  $(basename "$0") -t viz_bwt                # BWT HTML report only
  $(basename "$0") -t rgi_bwt_all            # Complete BWT pipeline
  $(basename "$0") -t rgi_compare_all        # RGI main + BWT + comparison
  $(basename "$0") -t mge_all                # MGE only
  $(basename "$0") -t rgi_main_all -n        # Dry run RGI main
  $(basename "$0") -t all -c 32              # Complete pipeline on 32 cores
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
    # Medaka
    "medaka_all"
    # PlasmidFinder
    "plasmidfinder_all"
    # Global
    "all"
    # R analysis
    "r_analysis_rgi"
    "r_analysis_bwt"
    "r_analysis_all"
    # Visualization
    "viz_rgi"
    "viz_bwt"
    "viz_all"
    # RGI pipelines
    "rgi_main_all"
    "rgi_bwt_all"
    "rgi_compare_all"
    # MGE
    "mge_all"
    # Co-occurrence
    "cooccurrence_all"
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
    --use-conda \
    ${DRY_RUN} \
    -- "${TARGET}"

echo ""
echo "======================================"
[[ -n "${DRY_RUN}" ]] \
    && echo "  Dry run completed — no files written" \
    || echo "  Pipeline completed successfully"
echo "======================================"
echo ""
