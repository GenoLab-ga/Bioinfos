#!/bin/bash
################################################################################
# RICE GENOMICS WGS PIPELINE — 01_setup.sh
# Initialisation : arborescence, outils, génome de référence + indexation
#
# EXÉCUTION : bash 01_setup.sh
# DURÉE     : ~10 min (téléchargement + indexation BWA)
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

log_message "INFO" "=== ÉTAPE 1 : SETUP ==="
print_config

# --- Arborescence ---
init_project
check_disk_space 50

# --- Environnement conda ---
log_message "INFO" "Vérification de l'environnement : $MAMBA_ENV"
if eval "$CONDA_ACTIVATE" &> /dev/null; then
    log_message "INFO" "✓ Environnement activé"
else
    error_exit "Impossible d'activer l'environnement : $MAMBA_ENV"
fi
eval "$CONDA_ACTIVATE"

# --- Outils requis ---
log_message "INFO" "Vérification des outils..."
REQUIRED_TOOLS=("bwa" "samtools" "bcftools" "bedtools" "vcftools" "fastp" "multiqc" "plink2" "snpEff")
MISSING=()
for tool in "${REQUIRED_TOOLS[@]}"; do
    if check_tool "$tool"; then
        log_message "INFO" "  ✓ $tool"
    else
        MISSING+=("$tool")
    fi
done

if [ ${#MISSING[@]} -gt 0 ]; then
    log_message "WARN" "Outils manquants : ${MISSING[*]}"
    log_message "WARN" "Installer avec : mamba install -c bioconda ${MISSING[*]}"
fi

# --- Génome de référence ---
log_message "INFO" "Vérification du génome de référence..."

if [ ! -f "$REF_FASTA" ]; then
    log_message "INFO" "Téléchargement IRGSP v1.0 depuis Ensembl Plants..."
    cd "$REF_DIR"
    wget -c "$REF_GENOME_URL" -O IRGSP-1.0_genome.fasta.gz
    gunzip IRGSP-1.0_genome.fasta.gz
    log_message "INFO" "✓ Génome téléchargé : $REF_FASTA"
else
    log_message "INFO" "✓ Génome déjà présent : $REF_FASTA"
fi

# --- Index BWA ---
# IMPORTANT : créer les liens symboliques pour que BWA trouve l'index
if [ ! -f "$REF_INDEX" ]; then
    log_message "INFO" "Indexation BWA (peut prendre ~5 min)..."
    bwa index -p "$REF_DIR/IRGSP-1.0_genome" "$REF_FASTA" 2>> "$LOG_FILE"
    log_message "INFO" "✓ Index BWA créé"
else
    log_message "INFO" "✓ Index BWA déjà présent"
fi

# Créer les liens symboliques .fasta.* pour que BWA trouve l'index via le chemin FASTA
for ext in amb ann bwt pac sa; do
    if [ ! -f "$REF_DIR/IRGSP-1.0_genome.fasta.${ext}" ]; then
        ln -sf "$REF_DIR/IRGSP-1.0_genome.${ext}" \
               "$REF_DIR/IRGSP-1.0_genome.fasta.${ext}"
    fi
done
log_message "INFO" "✓ Liens symboliques BWA créés"

# --- Dictionnaire samtools ---
if [ ! -f "$REF_DICT" ]; then
    log_message "INFO" "Création du dictionnaire samtools..."
    samtools dict "$REF_FASTA" > "$REF_DICT"
    log_message "INFO" "✓ Dictionnaire créé"
fi

# --- Index samtools FAI ---
if [ ! -f "$REF_FASTA.fai" ]; then
    log_message "INFO" "Création de l'index FAI..."
    samtools faidx "$REF_FASTA"
    log_message "INFO" "✓ Index FAI créé"
fi

# --- Base SnpEff ---
log_message "INFO" "Vérification de la base SnpEff Oryza_sativa..."
SNPEFF_DATA_DIR=$(dirname "$SNPEFF_CFG")/data
if [ ! -d "$SNPEFF_DATA_DIR/Oryza_sativa" ]; then
    log_message "INFO" "Téléchargement de la base SnpEff Oryza_sativa..."
    snpEff download "$SNPEFF_DB" -config "$SNPEFF_CFG" 2>> "$LOG_FILE" || \
        log_message "WARN" "Téléchargement SnpEff échoué — relancer manuellement"
else
    log_message "INFO" "✓ Base SnpEff Oryza_sativa présente"
fi

log_message "INFO" "=== SETUP TERMINÉ ==="
log_message "INFO" "Étape suivante : bash 02_download_data.sh"
