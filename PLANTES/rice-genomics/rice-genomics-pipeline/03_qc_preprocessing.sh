#!/bin/bash
################################################################################
# RICE GENOMICS WGS PIPELINE — 03_qc_preprocessing.sh
# Contrôle qualité et trimming avec fastp
#
# EXÉCUTION : bash 03_qc_preprocessing.sh
# DURÉE     : ~15-20 min pour 8 samples (5.7 Go)
#
# SORTIES :
#   results/01_qc/          — rapports HTML + JSON fastp par sample
#   results/01_qc/fastp_summary.tsv — tableau récapitulatif QC
#   results/01_qc/multiqc/  — rapport MultiQC agrégé
#   results/02_trimmed/     — FASTQ nettoyés
#   data/sample_info_trimmed.tsv — chemins vers FASTQ trimmés
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

log_message "INFO" "=== ÉTAPE 3 : QC ET PREPROCESSING ==="
eval "$CONDA_ACTIVATE"

check_tool "fastp" || error_exit "fastp introuvable — mamba install -c bioconda fastp"
log_message "INFO" "fastp : $(fastp --version 2>&1 | head -1)"

[ -f "$SAMPLE_INFO" ] || error_exit "sample_info.tsv introuvable — lancer 02_download_data.sh"

SAMPLE_COUNT=$(tail -n +2 "$SAMPLE_INFO" | grep -v "^$" | wc -l)
log_message "INFO" "$SAMPLE_COUNT samples à traiter"

mkdir -p "$QC_DIR" "$TRIM_DIR"

PROCESSED=0

while IFS=$'\t' read -r sample_acc run_acc country variety R1 R2; do
    [ "$sample_acc" = "sample_accession" ] && continue
    [[ -z "$sample_acc" ]] && continue

    PROCESSED=$((PROCESSED+1))
    log_message "INFO" "[$PROCESSED/$SAMPLE_COUNT] $sample_acc ($country — $variety)"

    R1_TRIM="$TRIM_DIR/${sample_acc}_R1.fastq.gz"
    R2_TRIM="$TRIM_DIR/${sample_acc}_R2.fastq.gz"
    HTML="$QC_DIR/${sample_acc}_fastp.html"
    JSON="$QC_DIR/${sample_acc}_fastp.json"

    if [ -f "$R1_TRIM" ] && [ -f "$R2_TRIM" ]; then
        log_message "INFO" "  ✓ Déjà trimmé"
        continue
    fi

    [ -f "$R1" ] && [ -f "$R2" ] || { log_message "WARN" "  ✗ FASTQ manquants"; continue; }

    fastp \
        --in1 "$R1" --in2 "$R2" \
        --out1 "$R1_TRIM" --out2 "$R2_TRIM" \
        --html "$HTML" --json "$JSON" \
        --detect_adapter_for_pe \
        --cut_front --cut_tail \
        --cut_mean_quality "$FASTP_MEAN_QUALITY" \
        --length_required "$FASTP_MIN_LENGTH" \
        --thread "$FASTP_THREADS" \
        --compression 4 \
        2>> "$LOG_FILE" || { log_message "WARN" "  ✗ fastp échoué"; rm -f "$R1_TRIM" "$R2_TRIM"; continue; }

    if [ -f "$JSON" ]; then
        R_BEFORE=$(python3 -c "import json; d=json.load(open('$JSON')); print(d['summary']['before_filtering']['total_reads'])" 2>/dev/null || echo "?")
        R_AFTER=$(python3 -c "import json; d=json.load(open('$JSON')); print(d['summary']['after_filtering']['total_reads'])" 2>/dev/null || echo "?")
        Q30=$(python3 -c "import json; d=json.load(open('$JSON')); print(round(d['summary']['after_filtering']['q30_rate']*100,1))" 2>/dev/null || echo "?")
        log_message "INFO" "  ✓ Reads : $R_BEFORE → $R_AFTER | Q30 : $Q30%"
    fi

done < "$SAMPLE_INFO"

# --- Tableau récapitulatif QC ---
log_message "INFO" "Génération du tableau récapitulatif QC..."
STATS="$QC_DIR/fastp_summary.tsv"
echo -e "sample\tcountry\tvariety\treads_before\treads_after\tq20_rate\tq30_rate\tgc_content\tpct_retained" > "$STATS"

while IFS=$'\t' read -r sample_acc run_acc country variety R1 R2; do
    [ "$sample_acc" = "sample_accession" ] && continue
    [[ -z "$sample_acc" ]] && continue
    JSON="$QC_DIR/${sample_acc}_fastp.json"
    [ -f "$JSON" ] || continue
    python3 << PYEOF >> "$STATS"
import json
d = json.load(open("$JSON"))
b = d["summary"]["before_filtering"]
a = d["summary"]["after_filtering"]
pct = round(a["total_reads"]/b["total_reads"]*100, 1) if b["total_reads"] > 0 else 0
print(f"$sample_acc\t$country\t$variety\t{b['total_reads']}\t{a['total_reads']}\t{round(a['q20_rate']*100,2)}\t{round(a['q30_rate']*100,2)}\t{round(a['gc_content']*100,2)}\t{pct}")
PYEOF
done < "$SAMPLE_INFO"

# --- MultiQC ---
if check_tool "multiqc"; then
    log_message "INFO" "Génération du rapport MultiQC..."
    multiqc "$QC_DIR" --outdir "$QC_DIR/multiqc" --filename "multiqc_report" --force 2>> "$LOG_FILE"
    log_message "INFO" "✓ MultiQC : $QC_DIR/multiqc/multiqc_report.html"
fi

# --- sample_info_trimmed.tsv ---
echo -e "sample_accession\trun_accession\tcountry\tvariety_group\tR1_trim\tR2_trim" > "$SAMPLE_INFO_TRIM"
while IFS=$'\t' read -r sample_acc run_acc country variety R1 R2; do
    [ "$sample_acc" = "sample_accession" ] && continue
    [[ -z "$sample_acc" ]] && continue
    R1T="$TRIM_DIR/${sample_acc}_R1.fastq.gz"
    R2T="$TRIM_DIR/${sample_acc}_R2.fastq.gz"
    [ -f "$R1T" ] && [ -f "$R2T" ] && \
        echo -e "${sample_acc}\t${run_acc}\t${country}\t${variety}\t${R1T}\t${R2T}" >> "$SAMPLE_INFO_TRIM"
done < "$SAMPLE_INFO"

TRIM_SIZE=$(du -sh "$TRIM_DIR" | cut -f1)
log_message "INFO" "=== QC TERMINÉ ==="
log_message "INFO" "Taille trimmés : $TRIM_SIZE"
log_message "INFO" "Résumé QC      : $STATS"
log_message "INFO" "Étape suivante : bash 04_mapping.sh"
