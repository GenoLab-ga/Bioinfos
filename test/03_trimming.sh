#!/bin/bash
# =============================================================
# Étape 3 : Trimming & filtrage qualité — fastp
# Projet  : PRJNA1465284 – P. falciparum drug resistance
# Auteur  : GenoLabGab
# Usage   : bash 03_trimming.sh
# =============================================================

set -euo pipefail

# ── Paramètres ───────────────────────────────────────────────
FASTQ_DIR="data/fastq_files"
TRIM_DIR="data/trimmed_files"
REPORT_DIR="results/trimming"
EXCLUDE_FILE="data/excluded_samples.txt"
THREADS=6
LOG="logs/03_trimming.log"

# ── Création des dossiers ─────────────────────────────────────
mkdir -p "${TRIM_DIR}" "${REPORT_DIR}" logs

echo "=======================================" | tee -a "${LOG}"
echo " Trimming démarré : $(date)"             | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"

# ── Vérification outil ────────────────────────────────────────
if ! command -v fastp &>/dev/null; then
    echo "ERREUR : fastp introuvable." | tee -a "${LOG}"
    exit 1
fi

# ── Charger les exclusions ────────────────────────────────────
declare -A EXCLUDED
while IFS= read -r excl || [ -n "$excl" ]; do
    [[ -z "$excl" || "$excl" =~ ^# ]] && continue
    EXCLUDED["$excl"]=1
done < "${EXCLUDE_FILE}"

echo "Échantillons exclus : ${#EXCLUDED[@]}" | tee -a "${LOG}"

# ── Compteurs ─────────────────────────────────────────────────
PROCESSED=0
SKIPPED=0
ERRORS=0

# ── Boucle principale ─────────────────────────────────────────
# Récupère les IDs uniques depuis les R1
for R1 in "${FASTQ_DIR}"/*_1.fastq.gz; do
    SAMPLE=$(basename "${R1}" _1.fastq.gz)
    R2="${FASTQ_DIR}/${SAMPLE}_2.fastq.gz"

    # Vérifier exclusion
    if [[ -v EXCLUDED["${SAMPLE}"] ]]; then
        echo "[SKIP] ${SAMPLE} — exclu (qualité dégradée)" | tee -a "${LOG}"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Vérifier que R2 existe
    if [ ! -f "${R2}" ]; then
        echo "[WARN] ${SAMPLE} — R2 manquant, ignoré" | tee -a "${LOG}"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Vérifier si déjà traité
    if [ -f "${TRIM_DIR}/${SAMPLE}_1.trimmed.fastq.gz" ]; then
        echo "[DONE] ${SAMPLE} — déjà trimmé, skip" | tee -a "${LOG}"
        PROCESSED=$((PROCESSED + 1))
        continue
    fi

    echo "[RUN] ${SAMPLE}..." | tee -a "${LOG}"

    fastp \
        --in1  "${R1}" \
        --in2  "${R2}" \
        --out1 "${TRIM_DIR}/${SAMPLE}_1.trimmed.fastq.gz" \
        --out2 "${TRIM_DIR}/${SAMPLE}_2.trimmed.fastq.gz" \
        \
        --detect_adapter_for_pe \
        \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --n_base_limit 5 \
        \
        --length_required 50 \
        \
        --cut_tail \
        --cut_tail_window_size 4 \
        --cut_tail_mean_quality 20 \
        \
        --thread "${THREADS}" \
        \
        --html "${REPORT_DIR}/${SAMPLE}_fastp.html" \
        --json "${REPORT_DIR}/${SAMPLE}_fastp.json" \
        \
        2>>"${LOG}" \
    || {
        echo "[ERREUR] ${SAMPLE}" | tee -a "${LOG}"
        ERRORS=$((ERRORS + 1))
        continue
    }

    PROCESSED=$((PROCESSED + 1))
    echo "[OK] ${SAMPLE}" | tee -a "${LOG}"

done

# ── Rapport MultiQC sur les JSON fastp ────────────────────────
echo "" | tee -a "${LOG}"
echo "Génération MultiQC post-trimming..." | tee -a "${LOG}"

multiqc "${REPORT_DIR}"/*.json \
    --outdir "${REPORT_DIR}/multiqc_trimming" \
    --filename "multiqc_trimming_PRJNA1465284" \
    --title "P. falciparum – Post-trimming QC" \
    --force \
    2>>"${LOG}"

# ── Résumé final ──────────────────────────────────────────────
echo "" | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"
echo " Trimming terminé : $(date)"             | tee -a "${LOG}"
echo " Traités  : ${PROCESSED}"                | tee -a "${LOG}"
echo " Exclus   : ${SKIPPED}"                  | tee -a "${LOG}"
echo " Erreurs  : ${ERRORS}"                   | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"
