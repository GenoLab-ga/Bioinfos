#!/bin/bash
################################################################################
# RICE GENOMICS WGS PIPELINE — config.sh
# Configuration centrale partagée par tous les scripts
# Projet : 3000 Rice Genomes (PRJEB6180) — Génétique des populations
# Auteur : Karl Mounguele | GenoLabGab
################################################################################
#
# USAGE : ce fichier est sourcé automatiquement par chaque script
#   source "$SCRIPT_DIR/config.sh"
#
# SEUL FICHIER À MODIFIER pour adapter le pipeline à une nouvelle machine.
################################################################################

# ==============================================================================
# SECTION 1 — ARBORESCENCE DU PROJET
# Modifier PROJECT_DIR selon votre machine
# ==============================================================================

PROJECT_DIR="$HOME/Documents/github_projet/PLANTES/rice-genomics"

DATA_DIR="$PROJECT_DIR/data"
RESULTS_DIR="$PROJECT_DIR/results"
LOGS_DIR="$PROJECT_DIR/logs"
SCRIPTS_DIR="$PROJECT_DIR/scripts"

REF_DIR="$DATA_DIR/reference"
RAW_FASTQ_DIR="$DATA_DIR/raw_fastq"

QC_DIR="$RESULTS_DIR/01_qc"
TRIM_DIR="$RESULTS_DIR/02_trimmed"
BAM_DIR="$RESULTS_DIR/03_bam"
VCF_DIR="$RESULTS_DIR/04_vcf"
ANNO_DIR="$RESULTS_DIR/05_annotation"
ANALYSIS_DIR="$RESULTS_DIR/06_analysis"


# ==============================================================================
# SECTION 2 — ENVIRONNEMENT CONDA/MAMBA
# ==============================================================================

MAMBA_ENV="genomics_env"
CONDA_ACTIVATE="source ~/miniforge3/etc/profile.d/mamba.sh && mamba activate $MAMBA_ENV"


# ==============================================================================
# SECTION 3 — GÉNOME DE RÉFÉRENCE
# Oryza sativa japonica — IRGSP v1.0
# Téléchargé depuis Ensembl Plants release-60
# ==============================================================================

REF_GENOME_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz"
REF_ANNOTATION_URL="https://rapdb.dna.affrc.go.jp/download/irgsp1/IRGSP-1.0_representative_2024-09.gff3.gz"

REF_FASTA="$REF_DIR/IRGSP-1.0_genome.fasta"
REF_GFF="$REF_DIR/IRGSP-1.0.gff3"
REF_INDEX="$REF_DIR/IRGSP-1.0_genome.bwt"   # index BWA classique
REF_DICT="$REF_DIR/IRGSP-1.0_genome.dict"


# ==============================================================================
# SECTION 4 — FICHIERS DES ÉCHANTILLONS
# ==============================================================================

SAMPLES_FILE="$PROJECT_DIR/samples.tsv"
SAMPLE_INFO="$DATA_DIR/sample_info.tsv"
SAMPLE_INFO_TRIM="$DATA_DIR/sample_info_trimmed.tsv"

# ENA metadata TSV (téléchargé lors du setup)
ENA_TSV="$DATA_DIR/PRJEB6180_full.tsv"


# ==============================================================================
# SECTION 5 — RESSOURCES COMPUTATIONNELLES
# Configuré pour AMD Ryzen 5 PRO — 12 CPUs, 16 Go RAM
# ==============================================================================

MAX_THREADS=10
MAX_MEMORY="14G"
BWA_THREADS=10
SAMTOOLS_THREADS=6
BCFTOOLS_THREADS=6
FASTP_THREADS=6


# ==============================================================================
# SECTION 6 — PARAMÈTRES QC (fastp)
# ==============================================================================

FASTP_MIN_QUALITY=20
FASTP_MIN_LENGTH=50
FASTP_MEAN_QUALITY=20


# ==============================================================================
# SECTION 7 — PARAMÈTRES VARIANT CALLING (bcftools)
# ==============================================================================

VCF_QUAL_THRESHOLD=30
VCF_DP_MIN=5
VCF_DP_MAX=100
MAF_THRESHOLD=0.05


# ==============================================================================
# SECTION 8 — ANNOTATION (SnpEff)
# ==============================================================================

SNPEFF_DB="Oryza_sativa"
SNPEFF_CFG=$(find ~/miniforge3/envs/${MAMBA_ENV}/ -name "snpEff.config" 2>/dev/null | head -1)


# ==============================================================================
# SECTION 9 — LOGGING
# ==============================================================================

mkdir -p "$LOGS_DIR"
LOG_FILE="$LOGS_DIR/pipeline_$(date +%Y%m%d_%H%M%S).log"
touch "$LOG_FILE"


# ==============================================================================
# SECTION 10 — FONCTIONS UTILITAIRES
# Disponibles dans tous les scripts via : source config.sh
# ==============================================================================

# log_message LEVEL "message" — affiche + écrit dans le log
log_message() {
    local level="$1"; shift
    local message="$*"
    local timestamp; timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] [$level] $message" | tee -a "$LOG_FILE"
}

# error_exit "message" — log erreur et arrête le script
error_exit() {
    log_message "ERROR" "$1"
    exit 1
}

# check_tool "outil" — vérifie qu'un outil est dans le PATH
check_tool() {
    local tool="$1"
    if ! command -v "$tool" &> /dev/null; then
        log_message "WARN" "Outil introuvable : $tool"
        return 1
    fi
    return 0
}

# create_dir "/chemin" — crée un dossier s'il n'existe pas
create_dir() {
    local dir="$1"
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        log_message "INFO" "Dossier créé : $dir"
    fi
}

# check_disk_space MIN_GB — vérifie l'espace disque disponible
check_disk_space() {
    local min_gb="$1"
    local available_gb
    available_gb=$(df "$PROJECT_DIR" | awk 'NR==2 {print int($4/1024/1024)}')
    if [ "$available_gb" -lt "$min_gb" ]; then
        error_exit "Espace insuffisant. Requis : ${min_gb}Go, Disponible : ${available_gb}Go"
    fi
    log_message "INFO" "Espace disque OK : ${available_gb}Go disponibles"
}

# init_project — crée toute l'arborescence du projet
init_project() {
    log_message "INFO" "=== Initialisation de l'arborescence ==="
    for dir in "$DATA_DIR" "$RESULTS_DIR" "$LOGS_DIR" "$SCRIPTS_DIR" \
               "$REF_DIR" "$RAW_FASTQ_DIR" "$QC_DIR" "$TRIM_DIR" \
               "$BAM_DIR" "$VCF_DIR" "$ANNO_DIR" "$ANALYSIS_DIR"; do
        create_dir "$dir"
    done
    log_message "INFO" "Arborescence initialisée"
}

# print_config — affiche la configuration courante
print_config() {
    cat << EOF
================================================================================
              RICE GENOMICS WGS PIPELINE — CONFIGURATION
================================================================================
Projet          : $PROJECT_DIR
Environnement   : $MAMBA_ENV
Référence       : $REF_FASTA
CPUs max        : $MAX_THREADS  |  RAM max : $MAX_MEMORY
Espace disque   : $(df -h "$HOME" | awk 'NR==2 {print $4}') disponibles
================================================================================
EOF
}


# ==============================================================================
# EXPORTS
# ==============================================================================

export PROJECT_DIR DATA_DIR RESULTS_DIR LOGS_DIR SCRIPTS_DIR
export REF_DIR RAW_FASTQ_DIR QC_DIR TRIM_DIR BAM_DIR VCF_DIR ANNO_DIR ANALYSIS_DIR
export SAMPLES_FILE SAMPLE_INFO SAMPLE_INFO_TRIM ENA_TSV
export MAMBA_ENV CONDA_ACTIVATE
export MAX_THREADS MAX_MEMORY BWA_THREADS SAMTOOLS_THREADS BCFTOOLS_THREADS FASTP_THREADS
export REF_FASTA REF_GFF REF_INDEX REF_DICT
export VCF_QUAL_THRESHOLD VCF_DP_MIN VCF_DP_MAX MAF_THRESHOLD
export SNPEFF_DB SNPEFF_CFG LOG_FILE
