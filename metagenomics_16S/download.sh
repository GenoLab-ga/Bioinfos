#!/bin/bash
# ============================================================
# download.sh — Téléchargement complet PRJNA417939
# Auteur : Karl Mounguele / GenoLabGab
# Date   : $(date +%F)
# ============================================================

set -euo pipefail  # arrêt immédiat si erreur

# --- Chemins -------------------------------------------------
PROJECT_DIR="$HOME/Documents/github_projet/metagenomics_16S"
RAW_DIR="$PROJECT_DIR/raw_data"
FASTQ_DIR="$PROJECT_DIR/fastq"
META_DIR="$PROJECT_DIR/metadata"
LOG_DIR="$PROJECT_DIR/logs"
SRR_LIST="$META_DIR/srr_list.txt"

# --- Création des dossiers -----------------------------------
mkdir -p "$RAW_DIR" "$FASTQ_DIR" "$META_DIR" "$LOG_DIR"
echo "[$(date +%T)] Dossiers créés."

# --- Étape 1 : Extraction des accessions SRR -----------------
echo "[$(date +%T)] Extraction des SRR depuis SraRunTable.csv..."
cut -d',' -f1 "$META_DIR/SraRunTable.csv" | tail -n +2 > "$SRR_LIST"
echo "[$(date +%T)] $(wc -l < "$SRR_LIST") accessions extraites."

# --- Étape 2 : prefetch --------------------------------------
echo "[$(date +%T)] Démarrage prefetch (4 jobs parallèles)..."
cat "$SRR_LIST" | parallel -j 4 --bar \
  "prefetch {} --output-directory $RAW_DIR 2>> $LOG_DIR/prefetch_errors.log"
echo "[$(date +%T)] prefetch terminé."

# --- Étape 3 : fasterq-dump + gzip --------------------------
echo "[$(date +%T)] Conversion en FASTQ + compression gzip..."
cat "$SRR_LIST" | parallel -j 3 --bar \
  "fasterq-dump $RAW_DIR/{}/{}.sra \
   --split-files \
   --threads 4 \
   --outdir $FASTQ_DIR \
   2>> $LOG_DIR/fasterq_errors.log && \
   gzip $FASTQ_DIR/{}_1.fastq $FASTQ_DIR/{}_2.fastq"
echo "[$(date +%T)] Conversion terminée."

# --- Étape 4 : Vérification ----------------------------------
echo "[$(date +%T)] Vérification de l'intégrité..."
N_FASTQ=$(ls "$FASTQ_DIR"/*.fastq.gz 2>/dev/null | wc -l)
echo "[$(date +%T)] Fichiers FASTQ.gz présents : $N_FASTQ / 2212 attendus"

SUSPECTS=$(find "$FASTQ_DIR" -name "*.fastq.gz" -size -100k)
if [[ -n "$SUSPECTS" ]]; then
  echo "[$(date +%T)] ATTENTION — fichiers suspects (< 100kb) :"
  echo "$SUSPECTS"
else
  echo "[$(date +%T)] Aucun fichier suspect détecté."
fi

echo "[$(date +%T)] Script terminé."
