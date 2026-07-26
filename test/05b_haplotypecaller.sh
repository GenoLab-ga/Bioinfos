#!/bin/bash
# =============================================================
# Étape 5b : GATK HaplotypeCaller – mode GVCF par échantillon
# Projet   : PRJNA1465284 – P. falciparum drug resistance
# Auteur   : GenoLabGab
# =============================================================

set -euo pipefail

# ── Paramètres ───────────────────────────────────────────────
BAM_DIR="data/bam_files"
GVCF_DIR="data/gvcf_files"
REF="data/reference/Pf3D7.fasta"
INTERVALS="data/reference/resistance_loci.interval_list"
EXCLUDE_FILE="data/excluded_samples.txt"
THREADS=4
LOG="logs/05b_haplotypecaller.log"

mkdir -p "${GVCF_DIR}" logs

echo "=======================================" | tee -a "${LOG}"
echo " HaplotypeCaller démarré : $(date)"      | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"

# ── Vérifications ─────────────────────────────────────────────
if ! command -v gatk &>/dev/null; then
    echo "ERREUR : gatk introuvable." | tee -a "${LOG}"
    exit 1
fi

if [ ! -f "${INTERVALS}" ]; then
    echo "ERREUR : intervalles manquants. Lance d'abord 05a." | tee -a "${LOG}"
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
TOTAL=$(ls "${BAM_DIR}"/*.markdup.bam | wc -l)

echo "BAM disponibles : ${TOTAL}" | tee -a "${LOG}"

# ── Boucle HaplotypeCaller ────────────────────────────────────
for BAM in "${BAM_DIR}"/*.markdup.bam; do
    SAMPLE=$(basename "${BAM}" .markdup.bam)
    GVCF="${GVCF_DIR}/${SAMPLE}.g.vcf.gz"

    # Exclusion
    if [[ -v EXCLUDED["${SAMPLE}"] ]]; then
        echo "[SKIP] ${SAMPLE}" | tee -a "${LOG}"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Déjà traité
    if [ -f "${GVCF}" ] && [ -f "${GVCF}.tbi" ]; then
        echo "[DONE] ${SAMPLE} — GVCF existant, skip" | tee -a "${LOG}"
        PROCESSED=$((PROCESSED + 1))
        continue
    fi

    echo "" | tee -a "${LOG}"
    echo "[RUN] ${SAMPLE} ($(( PROCESSED + SKIPPED + 1 ))/${TOTAL})..." \
        | tee -a "${LOG}"

    gatk HaplotypeCaller \
        --input "${BAM}" \
        --output "${GVCF}" \
        --reference "${REF}" \
        --intervals "${INTERVALS}" \
        \
        --emit-ref-confidence GVCF \
        --sample-ploidy 1 \
        \
        --min-base-quality-score 20 \
        --minimum-mapping-quality 30 \
        \
        --native-pair-hmm-threads "${THREADS}" \
        \
        --tmp-dir /tmp \
        --verbosity ERROR \
        2>>"${LOG}" \
    || {
        echo "[ERREUR] HaplotypeCaller ${SAMPLE}" | tee -a "${LOG}"
        rm -f "${GVCF}" "${GVCF}.tbi"
        ERRORS=$((ERRORS + 1))
        continue
    }

    echo "[OK] ${SAMPLE}" | tee -a "${LOG}"
    PROCESSED=$((PROCESSED + 1))

done

# ── Résumé ────────────────────────────────────────────────────
echo "" | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"
echo " HaplotypeCaller terminé : $(date)"      | tee -a "${LOG}"
echo " Traités  : ${PROCESSED}"                | tee -a "${LOG}"
echo " Exclus   : ${SKIPPED}"                  | tee -a "${LOG}"
echo " Erreurs  : ${ERRORS}"                   | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"

# ── Générer la liste des GVCF pour phase 5c ───────────────────
echo "" | tee -a "${LOG}"
echo "Génération de la liste GVCF..." | tee -a "${LOG}"

ls "${GVCF_DIR}"/*.g.vcf.gz > "${GVCF_DIR}/gvcf_list.txt"
echo "GVCF listés : $(wc -l < ${GVCF_DIR}/gvcf_list.txt)" \
    | tee -a "${LOG}"
