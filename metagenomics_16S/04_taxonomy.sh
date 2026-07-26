#!/bin/bash
# ============================================================
# 04_taxonomy.sh — Classification taxonomique SILVA 138
# Environnement : qiime2-amplicon-2024.10
# ============================================================

set -euo pipefail

PROJECT_DIR="$HOME/Documents/github_projet/metagenomics_16S"
DADA2_DIR="$PROJECT_DIR/results/03_dada2"
TAX_DIR="$PROJECT_DIR/results/04_taxonomy"
DB_DIR="$PROJECT_DIR/databases"
LOG_DIR="$PROJECT_DIR/logs"
THREADS=12

mkdir -p "$TAX_DIR"

# Classifier pré-entraîné SILVA 138 V4 (515F/806R)
CLASSIFIER="$DB_DIR/silva-138-99-seqs-515-806.qza"

if [[ ! -f "$CLASSIFIER" ]]; then
  echo "[$(date +%T)] Téléchargement classifier SILVA 138 V4..."
  mkdir -p "$DB_DIR"
  wget https://data.qiime2.org/2024.10/common/silva-138-99-seqs-515-806.qza \
       -O "$CLASSIFIER"
fi

echo "[$(date +%T)] Classification taxonomique en cours..."
qiime feature-classifier classify-sklearn \
  --i-classifier "$CLASSIFIER" \
  --i-reads "$DADA2_DIR/rep-seqs.qza" \
  --p-n-jobs "$THREADS" \
  --p-confidence 0.7 \
  --o-classification "$TAX_DIR/taxonomy.qza" \
  --verbose 2>&1 | tee "$LOG_DIR/taxonomy.log"

echo "[$(date +%T)] Classification terminée."

# --- Visualisations -----------------------------------------
qiime metadata tabulate \
  --m-input-file "$TAX_DIR/taxonomy.qza" \
  --o-visualization "$TAX_DIR/taxonomy.qzv"

# Barplot taxonomique
qiime taxa barplot \
  --i-table "$DADA2_DIR/table.qza" \
  --i-taxonomy "$TAX_DIR/taxonomy.qza" \
  --m-metadata-file "$PROJECT_DIR/metadata/metadata_qiime2.tsv" \
  --o-visualization "$TAX_DIR/taxa-barplot.qzv"

# Filtrer mitochondries et chloroplastes
echo "[$(date +%T)] Filtrage mitochondries/chloroplastes..."
qiime taxa filter-table \
  --i-table "$DADA2_DIR/table.qza" \
  --i-taxonomy "$TAX_DIR/taxonomy.qza" \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table "$TAX_DIR/table-filtered.qza"

echo "[$(date +%T)] Script 04 terminé."
echo "→ qiime tools view $TAX_DIR/taxa-barplot.qzv"
