#!/usr/bin/env bash
# =============================================================================
# run_medaka.sh — external GPU-server Medaka polishing
# =============================================================================
# Medaka polishing requires a GPU; this pipeline's compute node does not have
# one. This script runs medaka_consensus on a remote GPU server for each
# sample, pulling assembly.fasta + cleaned reads from the main pipeline
# server, and pushing the polished consensus back into
# results/{sample}/medaka/ — the same location workflow/rules/medaka.smk
# would produce if run locally with GPU access.
#
# This script is run MANUALLY, outside the Snakemake DAG. It is not a
# Snakemake rule and is not invoked automatically by run_AMR_pipeline.sh.
#
# NOTE (open methodological question — see project docs): rgi.smk and
# mge.smk currently detect ARGs/MGEs on the raw Flye assembly
# (results/{sample}/assembly.fasta), while plasmidfinder.smk uses the
# Medaka-polished assembly produced here. Whether ARG/MGE detection should
# also run on the polished assembly for consistency is still open.
#
# Configuration — set as environment variables, or edit the defaults below:
#   MEDAKA_REMOTE       user@remote-gpu-server
#   MEDAKA_SSH_KEY       path to the SSH key for that server
#   MEDAKA_REMOTE_DIR    path to this repository on the main pipeline server
#   MEDAKA_ENV_PREFIX    path to the conda/venv environment providing
#                         medaka_consensus on the GPU server
#
# Usage:
#   MEDAKA_REMOTE=user@gpu.example.org \
#   MEDAKA_SSH_KEY=$HOME/.ssh/id_ed25519 \
#   MEDAKA_REMOTE_DIR=/path/to/MetaSTREAM \
#   MEDAKA_ENV_PREFIX=$HOME/envs/medaka-conda \
#     ./run_medaka.sh H1_R H2_R H3_R
# =============================================================================
set -euo pipefail

# ─── Adjust these ─────────────────────────────────────────────────────────
SAMPLES=("$@")   # pass sample IDs as arguments, e.g. ./run_medaka.sh H1_R H2_R
REMOTE="${MEDAKA_REMOTE:?Set MEDAKA_REMOTE, e.g. user@remote-gpu-server}"
SSH_KEY="${MEDAKA_SSH_KEY:?Set MEDAKA_SSH_KEY, e.g. \$HOME/.ssh/id_ed25519}"
REMOTE_DIR="${MEDAKA_REMOTE_DIR:?Set MEDAKA_REMOTE_DIR, path to the repo on the main server}"
MODEL="${MEDAKA_MODEL:-r1041_e82_400bps_hac_v5.2.0}"
THREADS="${MEDAKA_THREADS:-12}"
BATCH_SIZE="${MEDAKA_BATCH_SIZE:-100}"   # lower on CUDA OOM (e.g. RTX 3090 24GB: 100 is safe, try 150 to push it)
ENV_PREFIX="${MEDAKA_ENV_PREFIX:?Set MEDAKA_ENV_PREFIX, path to the medaka conda/venv on the GPU server}"
WORK_ROOT="${MEDAKA_WORK_ROOT:-$HOME/medaka_run}"
# ────────────────────────────────────────────────────────────────────────

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "Usage: $0 SAMPLE_ID [SAMPLE_ID ...]"
    exit 1
fi

echo "==> Activating Medaka environment"
export PATH="${ENV_PREFIX}/bin:${PATH}"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"   # reduces CUDA memory fragmentation
command -v medaka_consensus >/dev/null || { echo "ERROR: medaka_consensus not found in ${ENV_PREFIX}/bin"; exit 1; }

FAILED=()

for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "════════════════════════════════════════════════════"
    echo "Sample: ${SAMPLE}"
    echo "════════════════════════════════════════════════════"

    WORKDIR="${WORK_ROOT}/${SAMPLE}"
    mkdir -p "$WORKDIR"
    cd "$WORKDIR"

    # A failure on one sample must not stop the others.
    (
        set -euo pipefail

        echo "==> Fetching input files from ${REMOTE}"
        rsync -avP -e "ssh -i ${SSH_KEY}" "${REMOTE}:${REMOTE_DIR}/results/${SAMPLE}/assembly.fasta" assembly.fasta
        rsync -avP -e "ssh -i ${SSH_KEY}" "${REMOTE}:${REMOTE_DIR}/results/${SAMPLE}/clean/${SAMPLE}.renamed.fastq.gz" reads.fastq.gz

        echo "==> Medaka polishing of ${SAMPLE} (model: ${MODEL})"
        medaka_consensus \
            -i reads.fastq.gz \
            -d assembly.fasta \
            -o medaka \
            -t "${THREADS}" \
            -m "${MODEL}" \
            -b "${BATCH_SIZE}"

        cp medaka/consensus.fasta "medaka/${SAMPLE}_polished.fasta"
        echo "==> Polishing complete: medaka/${SAMPLE}_polished.fasta"

        echo "==> Pushing results back to ${REMOTE}"
        ssh -i "${SSH_KEY}" "${REMOTE}" "mkdir -p ${REMOTE_DIR}/results/${SAMPLE}/medaka"
        rsync -avP -e "ssh -i ${SSH_KEY}" medaka/ "${REMOTE}:${REMOTE_DIR}/results/${SAMPLE}/medaka/"
    ) || { echo "!! FAILED: ${SAMPLE} (see output above)"; FAILED+=("$SAMPLE"); }

done

echo ""
echo "════════════════════════════════════════════════════"
if [[ ${#FAILED[@]} -eq 0 ]]; then
    echo "All samples completed successfully: ${SAMPLES[*]}"
else
    echo "Failed samples: ${FAILED[*]}"
    exit 1
fi
