#!/bin/bash
################################################################################
# RICE GENOMICS WGS PIPELINE — 04_mapping.sh
# Alignement des reads sur IRGSP v1.0 avec BWA-MEM
# Workflow : BWA-MEM → samtools fixmate → sort → markdup
#
# EXÉCUTION : bash 04_mapping.sh
# DURÉE     : ~60-90 min pour 8 samples
#
# SORTIES :
#   results/03_bam/                       — BAM triés + indexés
#   results/03_bam/marked_duplicates/     — BAM dédupliqués
#   results/03_bam/bam_list.txt           — liste pour variant calling
#   results/alignment_stats.tsv          — statistiques d'alignement
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

log_message "INFO" "=== ÉTAPE 4 : MAPPING ==="
eval "$CONDA_ACTIVATE"

[ -f "$REF_INDEX" ] || error_exit "Index BWA introuvable — lancer 01_setup.sh"
[ -f "$SAMPLE_INFO_TRIM" ] || error_exit "sample_info_trimmed.tsv introuvable — lancer 03_qc_preprocessing.sh"

SAMPLE_COUNT=$(tail -n +2 "$SAMPLE_INFO_TRIM" | grep -v "^$" | wc -l)
log_message "INFO" "$SAMPLE_COUNT samples à aligner"
log_message "INFO" "Référence : $REF_FASTA"

mkdir -p "$BAM_DIR/marked_duplicates"

PROCESSED=0

while IFS=$'\t' read -r sample_acc run_acc country variety R1_TRIM R2_TRIM; do
    [ "$sample_acc" = "sample_accession" ] && continue
    [[ -z "$sample_acc" ]] && continue

    PROCESSED=$((PROCESSED+1))
    log_message "INFO" "[$PROCESSED/$SAMPLE_COUNT] $sample_acc ($country)"

    BAM_FILE="$BAM_DIR/${sample_acc}.bam"
    BAM_DEDUP="$BAM_DIR/marked_duplicates/${sample_acc}.dedup.bam"

    if [ -f "$BAM_DEDUP" ] && [ -f "${BAM_DEDUP}.bai" ]; then
        log_message "INFO" "  ✓ Déjà aligné"
        continue
    fi

    [ -f "$R1_TRIM" ] && [ -f "$R2_TRIM" ] || { log_message "WARN" "  ✗ FASTQ trimmés manquants"; continue; }

    # BWA-MEM → fixmate (requis pour markdup) → sort
    # -R : Read Group obligatoire pour bcftools/GATK
    log_message "INFO" "  → BWA-MEM + fixmate + sort ($BWA_THREADS threads)..."
    bwa mem -t "$BWA_THREADS" \
        -R "@RG\tID:${run_acc}\tSM:${sample_acc}\tPL:ILLUMINA\tLB:${sample_acc}" \
        "$REF_FASTA" "$R1_TRIM" "$R2_TRIM" \
        2>> "$LOG_FILE" \
        | samtools fixmate -@ "$SAMTOOLS_THREADS" -m - - \
        | samtools sort -@ "$SAMTOOLS_THREADS" -o "$BAM_FILE" \
        2>> "$LOG_FILE" || { log_message "WARN" "  ✗ BWA-MEM échoué"; rm -f "$BAM_FILE"; continue; }

    samtools index -@ "$SAMTOOLS_THREADS" "$BAM_FILE"

    # Statistiques d'alignement
    FLAGSTAT=$(samtools flagstat "$BAM_FILE")
    TOTAL=$(echo "$FLAGSTAT" | grep "in total" | awk '{print $1}')
    MAPPED=$(echo "$FLAGSTAT" | grep "mapped (" | head -1 | awk '{print $1}')
    PCT=$(echo "scale=1; $MAPPED*100/$TOTAL" | bc)
    log_message "INFO" "  → Mapping : $MAPPED/$TOTAL reads ($PCT%)"

    # Marquage des duplicats (nécessite fixmate -m en amont)
    log_message "INFO" "  → Marquage des duplicats..."
    samtools markdup -@ "$SAMTOOLS_THREADS" "$BAM_FILE" "$BAM_DEDUP" 2>> "$LOG_FILE"
    samtools index -@ "$SAMTOOLS_THREADS" "$BAM_DEDUP"

    log_message "INFO" "  ✓ $sample_acc terminé"

done < "$SAMPLE_INFO_TRIM"

# --- Créer la liste BAM pour le variant calling ---
BAM_LIST="$BAM_DIR/bam_list.txt"
find "$BAM_DIR/marked_duplicates/" -name "*.dedup.bam" | sort > "$BAM_LIST"
log_message "INFO" "Liste BAM : $BAM_LIST ($(wc -l < $BAM_LIST) fichiers)"

# --- Tableau des statistiques ---
STATS="$RESULTS_DIR/alignment_stats.tsv"
echo -e "sample\tcountry\tvariety\ttotal_reads\tmapped_reads\tpct_mapped" > "$STATS"

while IFS=$'\t' read -r sample_acc run_acc country variety R1 R2; do
    [ "$sample_acc" = "sample_accession" ] && continue
    BAM_DEDUP="$BAM_DIR/marked_duplicates/${sample_acc}.dedup.bam"
    [ -f "$BAM_DEDUP" ] || continue
    FLAGSTAT=$(samtools flagstat "$BAM_DEDUP")
    TOTAL=$(echo "$FLAGSTAT" | grep "in total" | awk '{print $1}')
    MAPPED=$(echo "$FLAGSTAT" | grep "mapped (" | head -1 | awk '{print $1}')
    PCT=$(echo "scale=1; $MAPPED*100/$TOTAL" | bc)
    echo -e "${sample_acc}\t${country}\t${variety}\t${TOTAL}\t${MAPPED}\t${PCT}" >> "$STATS"
done < "$SAMPLE_INFO_TRIM"

BAM_SIZE=$(du -sh "$BAM_DIR" | cut -f1)
log_message "INFO" "=== MAPPING TERMINÉ ==="
log_message "INFO" "Taille BAM     : $BAM_SIZE"
log_message "INFO" "Statistiques   : $STATS"
log_message "INFO" "Étape suivante : bash 05_variant_calling.sh"
