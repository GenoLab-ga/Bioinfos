#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
log_message "INFO" "=== Starting Read Mapping ==="
eval "$CONDA_ACTIVATE"
SAMPLE_INFO_TRIM="$DATA_DIR/sample_info_trimmed.tsv"
SAMPLE_COUNT=$(tail -n +2 "$SAMPLE_INFO_TRIM" | wc -l)
mkdir -p "$BAM_DIR/marked_duplicates"
PROCESSED=0
while IFS=$'\t' read -r sample_acc run_acc country variety R1_TRIM R2_TRIM; do
    [ "$sample_acc" = "sample_accession" ] && continue
    PROCESSED=$((PROCESSED + 1))
    log_message "INFO" "[$PROCESSED/$SAMPLE_COUNT] $sample_acc ($country)"
    BAM_FILE="$BAM_DIR/${sample_acc}.bam"
    BAM_DEDUP="$BAM_DIR/marked_duplicates/${sample_acc}.dedup.bam"
    if [ -f "$BAM_DEDUP" ] && [ -f "${BAM_DEDUP}.bai" ]; then
        log_message "INFO" "  ✓ Déjà fait"
        continue
    fi
    log_message "INFO" "  → BWA-MEM + fixmate + sort..."
    bwa mem -t "$BWA_THREADS" \
        -R "@RG\tID:${run_acc}\tSM:${sample_acc}\tPL:ILLUMINA\tLB:${sample_acc}" \
        "$REF_FASTA" "$R1_TRIM" "$R2_TRIM" 2>> "$LOG_FILE" \
        | samtools fixmate -@ "$SAMTOOLS_THREADS" -m - - \
        | samtools sort -@ "$SAMTOOLS_THREADS" -o "$BAM_FILE" 2>> "$LOG_FILE"
    samtools index -@ "$SAMTOOLS_THREADS" "$BAM_FILE"

    log_message "INFO" "  → Marquage duplicats..."
    samtools markdup -@ "$SAMTOOLS_THREADS" \
        "$BAM_FILE" "$BAM_DEDUP" 2>> "$LOG_FILE"
    samtools index -@ "$SAMTOOLS_THREADS" "$BAM_DEDUP"
    FLAGSTAT=$(samtools flagstat "$BAM_DEDUP")
    TOTAL=$(echo "$FLAGSTAT" | grep "in total" | awk '{print $1}')
    MAPPED=$(echo "$FLAGSTAT" | grep "mapped (" | head -1 | awk '{print $1}')
    PCT=$(echo "scale=1; $MAPPED*100/$TOTAL" | bc)
    log_message "INFO" "  ✓ $MAPPED/$TOTAL reads mappés ($PCT%)"
done < "$SAMPLE_INFO_TRIM"
find "$BAM_DIR/marked_duplicates/" -name "*.dedup.bam" | sort > "$BAM_DIR/bam_list.txt"
log_message "INFO" "=== MAPPING COMPLETE === $(wc -l < $BAM_DIR/bam_list.txt) BAM files"
log_message "INFO" "Étape suivante : bash 05_variant_calling.sh"
