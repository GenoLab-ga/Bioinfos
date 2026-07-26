#!/bin/bash
# ============================================================
# 02_qiime2_import.sh — Import FASTQ + visualisation qualité
# ============================================================

set -euo pipefail

PROJECT_DIR="$HOME/Documents/github_projet/metagenomics_16S"
META_DIR="$PROJECT_DIR/metadata"
QIIME_DIR="$PROJECT_DIR/results/02_qiime2"

mkdir -p "$QIIME_DIR"

echo "[$(date +%T)] Import des FASTQ dans QIIME2..."
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$META_DIR/manifest.tsv" \
  --output-path "$QIIME_DIR/demux_paired.qza" \
  --input-format PairedEndFastqManifestPhred33V2

echo "[$(date +%T)] Génération visualisation qualité..."
qiime demux summarize \
  --i-data "$QIIME_DIR/demux_paired.qza" \
  --o-visualization "$QIIME_DIR/demux_paired.qzv"

echo "[$(date +%T)] Import terminé."
echo "Visualise avec : qiime tools view $QIIME_DIR/demux_paired.qzv"
