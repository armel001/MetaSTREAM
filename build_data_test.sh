#!/usr/bin/env bash
# =============================================================================
# MetagenAMR — build resources/data_test/ from non-study Nanopore reads
# =============================================================================
# Pools one or more source fastq(.gz) files, then subsamples two small mock
# "test" samples for out-of-the-box pipeline validation (technical smoke
# test — not intended to produce biologically meaningful ARG/taxonomy
# results, only to confirm the DAG executes end-to-end after cloning).
#
# Requires seqkit (not a pipeline runtime dependency — install once for
# this data-prep step only):
#   mamba create -n data_prep -c bioconda -c conda-forge seqkit
#   conda activate data_prep
# =============================================================================
set -euo pipefail

# --- EDIT THESE ---------------------------------------------------------
SOURCE_FASTQ=(
  "/data/armel/mythesis-data/data/EPU02A.fastq.gz"
  "/data/armel/mythesis-data/data/EPU03A.fastq.gz"
)
OUT_DIR="resources/data_test"
N_READS_PER_SAMPLE=3000
SEED=42
# -------------------------------------------------------------------------

mkdir -p "${OUT_DIR}"
POOL="${OUT_DIR}/.pool.fastq"

echo "[1/3] Pooling source fastq..."
zcat "${SOURCE_FASTQ[@]}" > "${POOL}"

echo "[2/3] Subsampling two mock test samples (seqkit sample, n=${N_READS_PER_SAMPLE})..."
# Different seeds -> test_1 and test_2 don't share identical reads,
# closer to two genuinely distinct samples for a meaningful smoke test.
seqkit sample -s "${SEED}"     -n "${N_READS_PER_SAMPLE}" "${POOL}" -o "${OUT_DIR}/test_1.fastq.gz"
seqkit sample -s "$((SEED+1))" -n "${N_READS_PER_SAMPLE}" "${POOL}" -o "${OUT_DIR}/test_2.fastq.gz"

echo "[3/3] Cleaning up pooled intermediate..."
rm "${POOL}"

echo ""
echo "Done — verify these stay well under GitHub's 100 MB/file limit:"
du -sh "${OUT_DIR}"/test_*.fastq.gz
