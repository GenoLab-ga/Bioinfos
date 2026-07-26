#!/bin/bash
# =============================================================
# Étape 4a : Téléchargement et indexation référence Pf3D7
# Source   : PlasmoDB release 68
# =============================================================

set -euo pipefail

REF_DIR="data/reference"
LOG="logs/04a_reference.log"

mkdir -p "${REF_DIR}" logs

echo "=======================================" | tee -a "${LOG}"
echo " Référence démarrée : $(date)"           | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"

# ── Téléchargement génome Pf3D7 ───────────────────────────────
echo "Téléchargement Pf3D7 depuis PlasmoDB..." | tee -a "${LOG}"

wget -q --show-progress \
    "https://plasmodb.org/common/downloads/release-68/PfalciparumNF54/fasta/data/PlasmoDB-68_PfalciparumNF54_Genome.fasta" \
    -O "${REF_DIR}/Pf3D7.fasta" \
    2>>"${LOG}" \
|| {
    # Fallback NCBI si PlasmoDB indisponible
    echo "PlasmoDB indisponible, tentative NCBI..." | tee -a "${LOG}"
    wget -q --show-progress \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/765/GCF_000002765.6_GCA_000002765/GCF_000002765.6_GCA_000002765_genomic.fna.gz" \
        -O "${REF_DIR}/Pf3D7.fasta.gz" \
        2>>"${LOG}"
    gunzip "${REF_DIR}/Pf3D7.fasta.gz"
}

echo "Taille référence : $(du -sh ${REF_DIR}/Pf3D7.fasta)" | tee -a "${LOG}"

# ── Indexation BWA ────────────────────────────────────────────
echo "Indexation BWA..." | tee -a "${LOG}"
bwa index "${REF_DIR}/Pf3D7.fasta" 2>>"${LOG}"

# ── Index samtools ────────────────────────────────────────────
echo "Index samtools fai..." | tee -a "${LOG}"
samtools faidx "${REF_DIR}/Pf3D7.fasta" 2>>"${LOG}"

# ── Dictionnaire séquences (pour GATK) ───────────────────────
echo "Dictionnaire GATK..." | tee -a "${LOG}"
samtools dict "${REF_DIR}/Pf3D7.fasta" \
    -o "${REF_DIR}/Pf3D7.dict" \
    2>>"${LOG}"

echo "" | tee -a "${LOG}"
echo "Fichiers générés :" | tee -a "${LOG}"
ls -lh "${REF_DIR}/" | tee -a "${LOG}"

echo "=======================================" | tee -a "${LOG}"
echo " Référence prête : $(date)"              | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"
