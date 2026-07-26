#!/bin/bash
# ============================================================
# convert_sra.sh — Conversion .sra → .fastq.gz
# Auteur : Karl Mounguele / GenoLabGab
# ============================================================

set -euo pipefail

PROJECT_DIR="$HOME/Documents/github_projet/metagenomics_16S"
RAW_DIR="$PROJECT_DIR/raw_data"
FASTQ_DIR="$PROJECT_DIR/fastq"
LOG_DIR="$PROJECT_DIR/logs"

mkdir -p "$FASTQ_DIR"

# Liste tous les fichiers .sra trouvés
find "$RAW_DIR" -name "*.sra" > "$LOG_DIR/sra_files.txt"
echo "[$(date +%T)] $(wc -l < "$LOG_DIR/sra_files.txt") fichiers .sra trouvés."

# Conversion en parallèle : 3 jobs × 4 threads = 12 CPUs
echo "[$(date +%T)] Démarrage fasterq-dump..."
cat "$LOG_DIR/sra_files.txt" | parallel -j 3 --bar \
  "fasterq-dump {} \
   --split-files \
   --threads 4 \
   --outdir $FASTQ_DIR \
   2>> $LOG_DIR/fasterq_errors.log"

echo "[$(date +%T)] Conversion terminée. Compression gzip..."

# Compression en parallèle
find "$FASTQ_DIR" -name "*.fastq" | parallel -j 6 "gzip {}"

echo "[$(date +%T)] Compression terminée."

# Vérification finale
N=$(ls "$FASTQ_DIR"/*.fastq.gz 2>/dev/null | wc -l)
echo "[$(date +%T)] Fichiers .fastq.gz présents : $N / 2212 attendus"

# Fichiers suspects (< 100kb)
SUSPECTS=$(find "$FASTQ_DIR" -name "*.fastq.gz" -size -100k)
if [[ -n "$SUSPECTS" ]]; then
  echo "[$(date +%T)] ATTENTION — fichiers suspects :"
  echo "$SUSPECTS" | tee "$LOG_DIR/suspects.txt"
else
  echo "[$(date +%T)] Aucun fichier suspect. Tout est propre."
fi
