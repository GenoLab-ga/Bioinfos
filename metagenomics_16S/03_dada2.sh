#!/bin/bash
# ============================================================
# 03_dada2.sh — Débruitage DADA2 + génération table ASVs
# Environnement : qiime2-amplicon-2024.10
# ============================================================

set -euo pipefail

PROJECT_DIR="$HOME/Documents/github_projet/metagenomics_16S"
QIIME_DIR="$PROJECT_DIR/results/02_qiime2"
DADA2_DIR="$PROJECT_DIR/results/03_dada2"
LOG_DIR="$PROJECT_DIR/logs"
THREADS=12

mkdir -p "$DADA2_DIR"

# -----------------------------------------------------------
# PARAMÈTRES DE TRONCATURE — À AJUSTER après visualisation
# de demux_paired.qzv dans qiime tools view
# Région V4 (515F/806R) : lectures 251bp paired-end
# Règle : tronquer là où la qualité médiane descend < Q25
# -----------------------------------------------------------
TRUNC_LEN_F=230   # forward — ajuster après QC
TRUNC_LEN_R=200   # reverse — souvent plus court car qualité dégrade

echo "[$(date +%T)] Démarrage DADA2 denoising..."
echo "[$(date +%T)] Paramètres : trunc_f=$TRUNC_LEN_F trunc_r=$TRUNC_LEN_R threads=$THREADS"

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$QIIME_DIR/demux_paired.qza" \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f "$TRUNC_LEN_F" \
  --p-trunc-len-r "$TRUNC_LEN_R" \
  --p-n-threads "$THREADS" \
  --p-chimera-method consensus \
  --o-table "$DADA2_DIR/table.qza" \
  --o-representative-sequences "$DADA2_DIR/rep-seqs.qza" \
  --o-denoising-stats "$DADA2_DIR/denoising-stats.qza" \
  --verbose 2>&1 | tee "$LOG_DIR/dada2.log"

echo "[$(date +%T)] DADA2 terminé."

# --- Visualisations -----------------------------------------
echo "[$(date +%T)] Génération des visualisations..."

qiime feature-table summarize \
  --i-table "$DADA2_DIR/table.qza" \
  --o-visualization "$DADA2_DIR/table.qzv"

qiime feature-table tabulate-seqs \
  --i-data "$DADA2_DIR/rep-seqs.qza" \
  --o-visualization "$DADA2_DIR/rep-seqs.qzv"

qiime metadata tabulate \
  --m-input-file "$DADA2_DIR/denoising-stats.qza" \
  --o-visualization "$DADA2_DIR/denoising-stats.qzv"

echo "[$(date +%T)] Visualisations générées."
echo "→ qiime tools view $DADA2_DIR/table.qzv"
echo "→ qiime tools view $DADA2_DIR/denoising-stats.qzv"
