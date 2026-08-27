#!/bin/bash
# =============================================================
# Étape 8 : Génération rapport Quarto
# =============================================================

set -euo pipefail

REPORT_DIR="results/report"
mkdir -p "${REPORT_DIR}"

# Copier les figures dans le dossier rapport
cp -r results/population_genetics/figures "${REPORT_DIR}/"
cp results/annotation/missense_variants.tsv "${REPORT_DIR}/"

# Rendre le rapport
quarto render report_PRJNA1465284.qmd \
    --output-dir "${REPORT_DIR}" \
    --to html

echo "Rapport généré : ${REPORT_DIR}/report_PRJNA1465284.html"
