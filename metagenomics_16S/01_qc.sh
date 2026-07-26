#!/bin/bash
# ============================================================
# 01_qc.sh — Contrôle qualité des lectures brutes
# Auteur : Karl Mounguele / GenoLabGab
# ============================================================

set -euo pipefail

PROJECT_DIR="$HOME/Documents/github_projet/metagenomics_16S"
FASTQ_DIR="$PROJECT_DIR/fastq"
QC_DIR="$PROJECT_DIR/results/01_qc"
FASTP_DIR="$QC_DIR/fastp"
FASTQC_DIR="$QC_DIR/fastqc_raw"
MULTIQC_DIR="$QC_DIR/multiqc"
LOG_DIR="$PROJECT_DIR/logs"
THREADS=4

mkdir -p "$FASTP_DIR" "$FASTQC_DIR" "$MULTIQC_DIR"

# --- Étape 1 : FastQC sur les bruts -------------------------
echo "[$(date +%T)] FastQC sur lectures brutes..."
ls "$FASTQ_DIR"/*_1.fastq.gz | parallel -j 6 \
  "fastqc {} {= s/_1/_2/ =} \
   --outdir $FASTQC_DIR \
   --threads 2 \
   --quiet"
echo "[$(date +%T)] FastQC terminé."

# --- Étape 2 : fastp (trimming + filtrage) ------------------
echo "[$(date +%T)] Trimming avec fastp..."
ls "$FASTQ_DIR"/*_1.fastq.gz | sed 's/_1.fastq.gz//' | \
  xargs -P 3 -I {} bash -c "
    SRR=\$(basename {})
    fastp \
      --in1 {}_1.fastq.gz \
      --in2 {}_2.fastq.gz \
      --out1 $FASTP_DIR/\${SRR}_1.fastp.fastq.gz \
      --out2 $FASTP_DIR/\${SRR}_2.fastp.fastq.gz \
      --json $FASTP_DIR/\${SRR}.json \
      --html $FASTP_DIR/\${SRR}.html \
      --thread $THREADS \
      --length_required 100 \
      --qualified_quality_phred 20 \
      --detect_adapter_for_pe \
      --correction \
      --disable_trim_poly_g \
      2>> $LOG_DIR/fastp_errors.log
  "
echo "[$(date +%T)] fastp terminé."

# --- Étape 3 : MultiQC --------------------------------------
echo "[$(date +%T)] Génération rapport MultiQC..."
multiqc "$FASTQC_DIR" "$FASTP_DIR" \
  --outdir "$MULTIQC_DIR" \
  --filename "multiqc_report" \
  --title "PRJNA417939 — QC 16S Mali Malaria" \
  --quiet
echo "[$(date +%T)] MultiQC terminé. Rapport : $MULTIQC_DIR/multiqc_report.html"
