#!/bin/bash
# run_mge.sh — wrapper mobileOG-pl indépendant du répertoire courant
# Usage: bash run_mge.sh -i assembly.fasta -d db.dmnd -m metadata.csv
#        -k 15 -e 1e-20 -p 90 -q 90 --outdir results/sample/mge

# ─── Parse arguments ─────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)      INPUT="$2";    shift 2 ;;
    -d|--db)         DIAMOND="$2";  shift 2 ;;
    -m|--metadata)   METADATA="$2"; shift 2 ;;
    -k|--kvalue)     KVALUE="$2";   shift 2 ;;
    -e|--escore)     ESCORE="$2";   shift 2 ;;
    -p|--pidentvalue) PIDENTVALUE="$2"; shift 2 ;;
    -q|--queryscore) QUERYSCORE="$2"; shift 2 ;;
    --outdir)        OUTDIR="$2";   shift 2 ;;
    --pypath)        PYPATH="$2";   shift 2 ;;
    *) shift ;;
  esac
done

# ─── Chemins absolus ─────────────────────────────────────────────────────────
INPUT=$(readlink -f "$INPUT")
DIAMOND=$(readlink -f "$DIAMOND")
METADATA=$(readlink -f "$METADATA")
OUTDIR=$(readlink -f "$OUTDIR")
mkdir -p "$OUTDIR"

# Préfixe de base pour les outputs = nom du fichier input sans chemin
BASENAME=$(basename "$INPUT")      # assembly.fasta
PREFIX="${OUTDIR}/${BASENAME}"     # results/HSE_2/mge/assembly.fasta

# ─── Step 1 : Prodigal ───────────────────────────────────────────────────────
echo "[1/3] Prodigal — protein prediction"
prodigal \
    -i "$INPUT" \
    -p meta \
    -a "${PREFIX}.faa" \
    -q

# ─── Step 2 : DIAMOND blastp ─────────────────────────────────────────────────
echo "[2/3] DIAMOND blastp — alignment against mobileOG-db"
diamond blastp \
    -q "${PREFIX}.faa" \
    --db "$DIAMOND" \
    --outfmt 6 stitle qtitle pident bitscore slen evalue qlen sstart send qstart qend \
    -k "$KVALUE" \
    -o "${PREFIX}.tsv" \
    -e "$ESCORE" \
    --query-cover "$QUERYSCORE" \
    --id "$PIDENTVALUE" \
    --quiet

# ─── Step 3 : mobileOGs-pl-kyanite.py ───────────────────────────────────────
echo "[3/3] mobileOGs-pl-kyanite.py — classification"
python "$PYPATH" \
    --o "${PREFIX}" \
    --i "${PREFIX}.tsv" \
    -m "$METADATA"

echo "Done — outputs in ${OUTDIR}"
