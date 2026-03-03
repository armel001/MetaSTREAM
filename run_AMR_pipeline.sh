#!/usr/bin/env bash

################################################################################
# AMR Pipeline for Metagenomic Analysis
# Author       : Thibaut Armel Chérif GNIMADI
# Affiliation  : CERFIG
# Description  : Metagenomic analysis pipeline with Kraken2/Bracken and RGI
# Version      : 3.0
# Date         : 2026-03-03
################################################################################

set -euo pipefail  # Stop on error, undefined variables, or pipe failure

# ==============================================================================
# CONFIGURATION
# ==============================================================================

CORES=22
CONFIG_FILE="config/config.yaml"
SNAKEFILE="workflow/Snakefile"
RESULTS_DIR="results"
SUMMARY_DIR="${RESULTS_DIR}/summary"

# Default target
TARGET="all"

# ==============================================================================
# USAGE
# ==============================================================================

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Options:
  -t, --target TARGET   Snakemake target rule to execute (default: all)
  -c, --cores  N        Number of cores to use (default: ${CORES})
  -n, --dry-run         Dry run — show what would be executed without running
  -h, --help            Show this help message

Available targets:
  all               Base pipeline: QC, taxonomy, stats (no RGI)
  rgi_main_all      RGI main pipeline (assembly-based AMR detection)
  rgi_bwt_all       RGI BWT pipeline (read-based AMR detection)
  rgi_compare_all   Both RGI pipelines + comparison report

Examples:
  $(basename "$0")                        # Run base pipeline
  $(basename "$0") -t rgi_main_all        # Run RGI main only
  $(basename "$0") -t rgi_bwt_all         # Run RGI BWT only
  $(basename "$0") -t rgi_compare_all     # Run both RGI + comparison
  $(basename "$0") -t rgi_main_all -n     # Dry run RGI main
EOF
    exit 0
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

DRY_RUN=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -t|--target)
            TARGET="$2"
            shift 2
            ;;
        -c|--cores)
            CORES="$2"
            shift 2
            ;;
        -n|--dry-run)
            DRY_RUN="--dry-run"
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            usage
            ;;
    esac
done

# ==============================================================================
# VALIDATION
# ==============================================================================

VALID_TARGETS=("all" "rgi_main_all" "rgi_bwt_all" "rgi_compare_all")
VALID=false
for t in "${VALID_TARGETS[@]}"; do
    [[ "$TARGET" == "$t" ]] && VALID=true && break
done

if [[ "$VALID" == false ]]; then
    echo "ERROR: Invalid target '${TARGET}'"
    echo "Valid targets: ${VALID_TARGETS[*]}"
    exit 1
fi

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "ERROR: Config file not found: ${CONFIG_FILE}"
    exit 1
fi

if [[ ! -f "${SNAKEFILE}" ]]; then
    echo "ERROR: Snakefile not found: ${SNAKEFILE}"
    exit 1
fi

# ==============================================================================
# RUN
# ==============================================================================

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
    -j 1 \
    --use-conda \
    --rerun-triggers params \
    ${DRY_RUN} \
    -- "${TARGET}"

echo ""
echo "======================================"
[[ -n "${DRY_RUN}" ]] \
    && echo "  Dry run completed — no files written" \
    || echo "  Pipeline completed successfully"
echo "======================================"
echo ""
