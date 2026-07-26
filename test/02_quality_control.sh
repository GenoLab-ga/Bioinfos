#!/bin/bash


# ======================================================== #
	# ÉTAPE 2 : Contrôle qualité FastQC + MulticQC
	# Projet : PRJNA1465284 - P. Falciparum drug resistance
	# Auteur : GenoLabGab
# ======================================================== #

set -euo pipefail

# == Paramètres ========================================== #
FASTQ_DIR="data/fastq_files"
QC_DIR="results/qc"
FASTQC_DIR="${QC_DIR}/fastqc_reports"
MULTIQC_DIR="${QC_DIR}/multiqc_report"
THREADS=6
LOG="logs/02_qc.log"


# == création des dossiers =============================== #

mkdir -p "${FASTQC_DIR}" "${MULTIQC_DIR}" logs

# == ===================================================== #

echo "================================" | tee -a "${LOG}"
echo " QC demarré : $(data)" | tee -a "${LOG}"
echo "================================" | tee -a "${LOG}"


# === Vérification des outils =============================
for tool in fastqc multiqc; do
	if ! command -v "${tool}" &>/dev/null; then
		echo "ERREUR : ${tool} introuvable. Active l'environnement rnaseq_env." | tee -a "${LOG}"
		exit 1
	fi

done

# == Comptage des fichiers à traiter =====================
N_FILES=$(ls "${FASTQ_DIR}"/*.fastq.gz 2>/dev/null | wc -l)

if [ "${N_FILES}" -eq 0 ]; then
    echo "ERREUR : aucun fichier .fastq.gz trouvé dans ${FASTQ_DIR}" | tee -a "${LOG}"
    exit 1
fi

echo "Fichiers trouvés : ${N_FILES}" | tee -a "${LOG}"

# ── FastQC sur tous les fichiers ──────────────────────────────
echo "" | tee -a "${LOG}"
echo "Lancement FastQC sur ${N_FILES} fichiers (${THREADS} threads)..." | tee -a "${LOG}"

fastqc \
    "${FASTQ_DIR}"/*.fastq.gz \
    --outdir "${FASTQC_DIR}" \
    --threads "${THREADS}" \
    --quiet \
    2>>"${LOG}"

echo "FastQC terminé : $(date)" | tee -a "${LOG}"

# ── MultiQC — agrégation de tous les rapports ─────────────────
echo "" | tee -a "${LOG}"
echo "Génération du rapport MultiQC..." | tee -a "${LOG}"

multiqc "${FASTQC_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "multiqc_PRJNA1465284" \
    --title "P. falciparum – Drug Resistance QC – PRJNA1465284" \
    --comment "Amplicon PfSMARRTer – kelch13, crt, mdr1, dhfr, dhps" \
    --force \
    2>>"${LOG}"

echo "MultiQC terminé : $(date)" | tee -a "${LOG}"

# ── Détection automatique des échantillons suspects ───────────
echo "" | tee -a "${LOG}"
echo "--- Détection échantillons suspects ---" | tee -a "${LOG}"

FAIL_COUNT=0
for zip in "${FASTQC_DIR}"/*_fastqc.zip; do
    SAMPLE=$(basename "${zip}" _fastqc.zip)
    # Extraire le statut global depuis le rapport FastQC
    STATUS=$(unzip -p "${zip}" "*/summary.txt" 2>/dev/null | grep "FAIL" | wc -l)
    if [ "${STATUS}" -gt 2 ]; then
        echo "  SUSPECT : ${SAMPLE} (${STATUS} modules FAIL)" | tee -a "${LOG}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
done

echo "" | tee -a "${LOG}"
echo "Échantillons suspects : ${FAIL_COUNT}" | tee -a "${LOG}"
echo "Rapport MultiQC : ${MULTIQC_DIR}/multiqc_PRJNA1465284.html" | tee -a "${LOG}"
echo "" | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"
echo " QC terminé : $(date)" | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"
