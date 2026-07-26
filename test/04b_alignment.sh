#!/bin/bash
# =============================================================
# Étape 4b : Alignement BWA-MEM + fixmate + sort + markdup
# Projet   : PRJNA1465284 – P. falciparum drug resistance
# Auteur   : GenoLabGab
# =============================================================

set -euo pipefail

# ── Paramètres ───────────────────────────────────────────────
TRIM_DIR="data/trimmed_files"
BAM_DIR="data/bam_files"
REF="data/reference/Pf3D7.fasta"
EXCLUDE_FILE="data/excluded_samples.txt"
THREADS=6
LOG="logs/04b_alignment.log"

# ── Création dossiers ─────────────────────────────────────────
mkdir -p "${BAM_DIR}/tmp" logs

echo "=======================================" | tee -a "${LOG}"
echo " Alignement démarré : $(date)"           | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"

# ── Vérification outils ───────────────────────────────────────
for tool in bwa samtools; do
    if ! command -v "${tool}" &>/dev/null; then
        echo "ERREUR : ${tool} introuvable." | tee -a "${LOG}"
        exit 1
    fi
done

# ── Vérification référence ────────────────────────────────────
if [ ! -f "${REF}.bwt" ]; then
    echo "ERREUR : référence non indexée." | tee -a "${LOG}"
    exit 1
fi

# ── Charger exclusions ────────────────────────────────────────
declare -A EXCLUDED
while IFS= read -r excl || [ -n "$excl" ]; do
    [[ -z "$excl" || "$excl" =~ ^# ]] && continue
    EXCLUDED["$excl"]=1
done < "${EXCLUDE_FILE}"

# ── Compteurs ─────────────────────────────────────────────────
PROCESSED=0
SKIPPED=0
ERRORS=0
TOTAL=$(ls "${TRIM_DIR}"/*_1.trimmed.fastq.gz | wc -l)

echo "Échantillons à traiter : ${TOTAL}" | tee -a "${LOG}"

# ── Nettoyage fichiers incomplets d'une exécution précédente ──
echo "Nettoyage fichiers incomplets..." | tee -a "${LOG}"
for bam in "${BAM_DIR}"/*.sorted.bam \
           "${BAM_DIR}"/*.fixmate.bam; do
    [ -f "$bam" ] && rm -f "$bam"
done

# ── Boucle principale ─────────────────────────────────────────
for R1 in "${TRIM_DIR}"/*_1.trimmed.fastq.gz; do
    SAMPLE=$(basename "${R1}" _1.trimmed.fastq.gz)
    R2="${TRIM_DIR}/${SAMPLE}_2.trimmed.fastq.gz"
    FINAL_BAM="${BAM_DIR}/${SAMPLE}.markdup.bam"

    # ── Exclusions ────────────────────────────────────────────
    if [[ -v EXCLUDED["${SAMPLE}"] ]]; then
        echo "[SKIP] ${SAMPLE} — exclu" | tee -a "${LOG}"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # ── R2 manquant ───────────────────────────────────────────
    if [ ! -f "${R2}" ]; then
        echo "[WARN] ${SAMPLE} — R2 manquant" | tee -a "${LOG}"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # ── Déjà traité ───────────────────────────────────────────
    if [ -f "${FINAL_BAM}" ] && [ -f "${FINAL_BAM}.bai" ]; then
        echo "[DONE] ${SAMPLE} — déjà aligné, skip" | tee -a "${LOG}"
        PROCESSED=$((PROCESSED + 1))
        continue
    fi

    echo "" | tee -a "${LOG}"
    echo "[RUN] ${SAMPLE} ($(( PROCESSED + SKIPPED + 1 ))/${TOTAL})..." \
        | tee -a "${LOG}"

    RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}\tPU:${SAMPLE}"

    # ── Pipeline complet en un seul pipe ──────────────────────
    # bwa mem → fixmate (-m calcule ms tag) → sort → markdup
    set +e
    bwa mem \
        -t "${THREADS}" \
        -R "${RG}" \
        -M \
        "${REF}" \
        "${R1}" "${R2}" \
        2>>"${LOG}" \
    | samtools fixmate \
        -@ "${THREADS}" \
        -m \
        -O bam \
        - - \
        2>>"${LOG}" \
    | samtools sort \
        -@ "${THREADS}" \
        -T "${BAM_DIR}/tmp/${SAMPLE}" \
        -O bam \
        - \
        2>>"${LOG}" \
    | samtools markdup \
        -@ "${THREADS}" \
        - \
        "${FINAL_BAM}" \
        2>>"${LOG}"

    PIPE_STATUS=${PIPESTATUS[@]}
    set -e

    # ── Vérification statut du pipe ───────────────────────────
    PIPE_FAIL=0
    for s in ${PIPE_STATUS}; do
        [ "$s" -ne 0 ] && PIPE_FAIL=1
    done

    if [ "${PIPE_FAIL}" -eq 1 ] || [ ! -s "${FINAL_BAM}" ]; then
        echo "[ERREUR] pipe échoué pour ${SAMPLE} — statuts : ${PIPE_STATUS}" \
            | tee -a "${LOG}"
        rm -f "${FINAL_BAM}"
        ERRORS=$((ERRORS + 1))
        continue
    fi

    # ── Indexation BAM final ───────────────────────────────────
    samtools index "${FINAL_BAM}" 2>>"${LOG}"

    # ── Statistiques alignement ────────────────────────────────
    samtools flagstat "${FINAL_BAM}" \
        > "${BAM_DIR}/${SAMPLE}.flagstat.txt" \
        2>>"${LOG}"

    PCT=$(grep "mapped (" "${BAM_DIR}/${SAMPLE}.flagstat.txt" \
          | head -1 | awk '{print $5}' | tr -d '(')

    echo "[OK] ${SAMPLE} — aligné : ${PCT}" | tee -a "${LOG}"
    PROCESSED=$((PROCESSED + 1))

done

# ── Nettoyage tmp ─────────────────────────────────────────────
rm -rf "${BAM_DIR}/tmp"

# ── Résumé MultiQC ────────────────────────────────────────────
echo "" | tee -a "${LOG}"
echo "Génération MultiQC alignement..." | tee -a "${LOG}"

multiqc "${BAM_DIR}"/*.flagstat.txt \
    --outdir "results/alignment_qc" \
    --filename "multiqc_alignment_PRJNA1465284" \
    --title "P. falciparum – Alignment QC" \
    --force \
    2>>"${LOG}"

# ── Résumé final ──────────────────────────────────────────────
echo "" | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"
echo " Alignement terminé : $(date)"           | tee -a "${LOG}"
echo " Traités  : ${PROCESSED}"                | tee -a "${LOG}"
echo " Exclus   : ${SKIPPED}"                  | tee -a "${LOG}"
echo " Erreurs  : ${ERRORS}"                   | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"

# ── Alerte taux alignement < 70% ──────────────────────────────
echo "" | tee -a "${LOG}"
echo "=== Échantillons alignement < 70% ===" | tee -a "${LOG}"
ALERTS=0
for stat in "${BAM_DIR}"/*.flagstat.txt; do
    SAMPLE=$(basename "${stat}" .flagstat.txt)
    PCT_NUM=$(grep "mapped (" "${stat}" | head -1 \
              | grep -oP '\d+\.\d+(?=%)' | head -1)
    if [ -n "${PCT_NUM}" ] && \
       (( $(echo "${PCT_NUM} < 70" | bc -l) )); then
        echo "  ALERTE : ${SAMPLE} — ${PCT_NUM}%" | tee -a "${LOG}"
        ALERTS=$((ALERTS + 1))
    fi
done
[ "${ALERTS}" -eq 0 ] && \
    echo "  Aucun — tous les échantillons > 70%" | tee -a "${LOG}"
