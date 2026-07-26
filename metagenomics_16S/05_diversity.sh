#!/bin/bash
# ============================================================
# 05_diversity.sh — Diversité alpha et bêta
# Environnement : qiime2-amplicon-2024.10
# ============================================================

set -euo pipefail

PROJECT_DIR="$HOME/Documents/github_projet/metagenomics_16S"
TAX_DIR="$PROJECT_DIR/results/04_taxonomy"
DIV_DIR="$PROJECT_DIR/results/05_diversity"
META_DIR="$PROJECT_DIR/metadata"
LOG_DIR="$PROJECT_DIR/logs"
THREADS=12

mkdir -p "$DIV_DIR"

# -----------------------------------------------------------
# PROFONDEUR DE RARÉFACTION — À définir après table.qzv
# Règle : choisir la profondeur qui retient ~95% des samples
# Regarder le graphe "Interactive Sample Detail" dans table.qzv
# -----------------------------------------------------------
SAMPLING_DEPTH=5000  # ajuster après inspection de table.qzv

echo "[$(date +%T)] Arbre phylogénétique (nécessaire pour UniFrac)..."
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "$PROJECT_DIR/results/03_dada2/rep-seqs.qza" \
  --o-alignment "$DIV_DIR/aligned-rep-seqs.qza" \
  --o-masked-alignment "$DIV_DIR/masked-aligned-rep-seqs.qza" \
  --o-tree "$DIV_DIR/unrooted-tree.qza" \
  --o-rooted-tree "$DIV_DIR/rooted-tree.qza" \
  --p-n-threads "$THREADS" \
  2>&1 | tee "$LOG_DIR/phylogeny.log"

echo "[$(date +%T)] Analyses de diversité alpha et bêta..."
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "$DIV_DIR/rooted-tree.qza" \
  --i-table "$TAX_DIR/table-filtered.qza" \
  --p-sampling-depth "$SAMPLING_DEPTH" \
  --m-metadata-file "$META_DIR/metadata_qiime2.tsv" \
  --output-dir "$DIV_DIR/core-metrics-results" \
  --p-n-jobs-or-threads "$THREADS" \
  2>&1 | tee "$LOG_DIR/diversity.log"

echo "[$(date +%T)] Diversité alpha — tests statistiques..."

# Shannon
qiime diversity alpha-group-significance \
  --i-alpha-diversity "$DIV_DIR/core-metrics-results/shannon_vector.qza" \
  --m-metadata-file "$META_DIR/metadata_qiime2.tsv" \
  --o-visualization "$DIV_DIR/shannon-significance.qzv"

# Faith PD (phylogénétique)
qiime diversity alpha-group-significance \
  --i-alpha-diversity "$DIV_DIR/core-metrics-results/faith_pd_vector.qza" \
  --m-metadata-file "$META_DIR/metadata_qiime2.tsv" \
  --o-visualization "$DIV_DIR/faith-pd-significance.qzv"

echo "[$(date +%T)] Diversité bêta — tests PERMANOVA..."

# Bray-Curtis
qiime diversity beta-group-significance \
  --i-distance-matrix "$DIV_DIR/core-metrics-results/bray_curtis_distance_matrix.qza" \
  --m-metadata-file "$META_DIR/metadata_qiime2.tsv" \
  --m-metadata-column "malaria_status" \
  --p-method permanova \
  --p-pairwise \
  --o-visualization "$DIV_DIR/bray-curtis-significance.qzv"

# UniFrac pondéré
qiime diversity beta-group-significance \
  --i-distance-matrix "$DIV_DIR/core-metrics-results/weighted_unifrac_distance_matrix.qza" \
  --m-metadata-file "$META_DIR/metadata_qiime2.tsv" \
  --m-metadata-column "malaria_status" \
  --p-method permanova \
  --p-pairwise \
  --o-visualization "$DIV_DIR/weighted-unifrac-significance.qzv"

echo "[$(date +%T)] Script 05 terminé."
echo "Résultats : $DIV_DIR/core-metrics-results/"
