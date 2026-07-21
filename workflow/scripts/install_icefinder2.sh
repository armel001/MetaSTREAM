#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════════════
# Installation one-time d'ICEfinder2 (ICEberg 3.0 engine)
# A executer manuellement UNE SEULE FOIS, hors Snakemake
# ═══════════════════════════════════════════════════════════════════════════════
set -euo pipefail

INSTALL_DIR="/data/armel/mythesis-data/tools/icefinder2"
ENV_NAME="icefinder2"

mkdir -p "$(dirname "${INSTALL_DIR}")"

echo "==> 1. Creation de l'environnement conda avec les 8 dependances"
conda create -y --name "${ENV_NAME}" -c bioconda -c conda-forge \
    hmmer">=3.1" \
    blast">=2.10.1" \
    kraken2">=2.0.9" \
    seqkit"=0.12.0" \
    prodigal"=2.6.3" \
    prokka"=1.14.6" \
    macsyfinder \
    defense-finder \
    biopython \
    ete3 \
    python">=3.8"

echo "==> 2. Recuperation du code ICEfinder2 (mirroir EBI-Metagenomics)"
git clone https://github.com/EBI-Metagenomics/icefinder2.git "${INSTALL_DIR}"

echo "==> 3. Installation des bases DefenseFinder"
conda run -n "${ENV_NAME}" defense-finder update

echo "==> 4. Installation des modeles MacSyFinder (CONJscan)"
conda run -n "${ENV_NAME}" macsydata install --org macsy-models CONJScan || \
    echo "AVERTISSEMENT : verifier manuellement l'installation du modele CONJScan"

echo ""
echo "==================================================================="
echo "  Prochaine etape MANUELLE obligatoire :"
echo "  Editer ${INSTALL_DIR}/config.ini"
echo "  et renseigner les chemins vers chaque binaire installe dans"
echo "  l'environnement conda '${ENV_NAME}' (utiliser :"
echo "  conda run -n ${ENV_NAME} which <hmmsearch|blastn|kraken2|seqkit|prodigal|prokka|defense-finder|macsyfinder>"
echo "  pour obtenir chaque chemin exact)."
echo "==================================================================="
