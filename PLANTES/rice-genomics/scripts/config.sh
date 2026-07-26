#!/bin/bash

################################################################################
# RICE GENOMICS WGS ANALYSIS PIPELINE - config.sh
# Variables et fonctions partagées par tous les scripts du pipeline
# Project : 3000 Rice Genomes (3K-RGP) - Population Genetics
# Author  : Karl Mounguele | GenoLabGab
# Date    : 2026-04-18
#
# COMMENT LIRE CE FICHIER
# Ce fichier est "sourcé" par chaque script avec : source config.sh
# Cela signifie que toutes les variables ($PROJECT_DIR, $MAX_THREADS...)
# et toutes les fonctions (log_message, check_tool...) deviennent
# automatiquement disponibles dans le script qui l'appelle.
# Tu n'as normalement pas à modifier les autres scripts — seulement ce fichier.
################################################################################


# ==============================================================================
# SECTION 1 : ARBORESCENCE DU PROJET
# Toutes les variables de chemins sont définies ici une seule fois.
# ==============================================================================

PROJECT_DIR="$HOME/Documents/github_projet/PLANTES/rice-genomics"          # Racine du projet
DATA_DIR="$PROJECT_DIR/data"               # Données brutes + référence
RESULTS_DIR="$PROJECT_DIR/results"         # Tous les résultats
LOGS_DIR="$PROJECT_DIR/logs"               # Fichiers de log
SCRIPTS_DIR="$PROJECT_DIR/scripts"         # Scripts auxiliaires

# Sous-dossiers données
REF_DIR="$DATA_DIR/reference"              # Génome de référence
RAW_FASTQ_DIR="$DATA_DIR/raw_fastq"        # FASTQ bruts téléchargés

# Sous-dossiers résultats (numérotés = ordre du pipeline)
QC_DIR="$RESULTS_DIR/01_qc"               # Rapports fastp
TRIM_DIR="$RESULTS_DIR/02_trimmed"        # FASTQ nettoyés
BAM_DIR="$RESULTS_DIR/03_bam"             # Fichiers d'alignement
VCF_DIR="$RESULTS_DIR/04_vcf"             # Variants bruts et filtrés
ANNO_DIR="$RESULTS_DIR/05_annotation"     # Variants annotés
ANALYSIS_DIR="$RESULTS_DIR/06_analysis"   # Analyses population


# ==============================================================================
# SECTION 2 : ENVIRONNEMENT CONDA/MAMBA
# ==============================================================================

MAMBA_ENV="genomics_env"
# Cette commande active l'environnement dans un sous-shell bash
CONDA_ACTIVATE="source ~/miniforge3/etc/profile.d/mamba.sh && mamba activate $MAMBA_ENV"


# ==============================================================================
# SECTION 3 : GÉNOME DE RÉFÉRENCE
# Oryza sativa japonica — IRGSP v1.0 (standard international)
# C'est la référence japonica mais elle sert pour toutes les sous-espèces
# car c'est le génome le mieux annoté du riz cultivé.
# ==============================================================================

REF_GENOME_URL="https://rapdb.dna.affrc.go.jp/download/irgsp1/IRGSP-1.0_genome.fasta.gz"
REF_ANNOTATION_URL="https://rapdb.dna.affrc.go.jp/download/irgsp1/IRGSP-1.0_representative_2024-09.gff3.gz"

REF_FASTA="$REF_DIR/IRGSP-1.0_genome.fasta"
REF_GFF="$REF_DIR/IRGSP-1.0.gff3"
REF_INDEX="$REF_DIR/IRGSP-1.0_genome.bwt"  # Extension BWA-MEM2
REF_DICT="$REF_DIR/IRGSP-1.0_genome.dict"           # Dictionnaire GATK/samtools


# ==============================================================================
# SECTION 4 : FICHIER DES ÉCHANTILLONS
# Format TSV : sample_accession <TAB> run_accession <TAB> country <TAB> variety_group
# Les 8 échantillons sélectionnés manuellement depuis le projet PRJEB6180
# ==============================================================================

SAMPLES_FILE="$PROJECT_DIR/samples.tsv"

# Répertoire où les FASTQ sont déjà téléchargés (structure : fastq/SAMEA.../run_R1.fastq.gz)
FASTQ_DIR="$HOME/Documents/github_projet/PLANTES/Rice_analysis/fastq"


# ==============================================================================
# SECTION 5 : RESSOURCES COMPUTATIONNELLES
# Optimisé pour AMD Ryzen 5 PRO — 12 CPUs, 16 Go RAM
# Règle : on laisse toujours 2 CPUs et 2 Go pour le système
# ==============================================================================

MAX_THREADS=10        # Utilisé pour les outils qui acceptent -t ou --threads
MAX_MEMORY="14G"      # Utilisé pour Java (SnpEff) et samtools sort
BWA_THREADS=10        # BWA-MEM2 est très bien parallélisé → max threads
SAMTOOLS_THREADS=6    # Moins gourmand, 6 suffisent
BCFTOOLS_THREADS=6    # Idem
FASTP_THREADS=6       # fastp est rapide, 6 threads est largement suffisant


# ==============================================================================
# SECTION 6 : PARAMÈTRES QC ET TRIMMING (fastp)
# fastp remplace FastQC + Trimmomatic en une seule commande
# Il génère lui-même un rapport HTML + JSON
# ==============================================================================

FASTP_MIN_QUALITY=20       # Score Phred minimum par base
FASTP_MIN_LENGTH=50        # Longueur minimale après trimming (WGS riz : reads ~150bp)
FASTP_WINDOW_SIZE=4        # Fenêtre pour le sliding window quality
FASTP_MEAN_QUALITY=20      # Qualité moyenne minimale sur la fenêtre


# ==============================================================================
# SECTION 7 : PARAMÈTRES VARIANT CALLING (bcftools)
# ==============================================================================

VCF_QUAL_THRESHOLD=30    # Score qualité variant minimum (Phred)
VCF_DP_MIN=5             # Profondeur minimale (évite les faux positifs)
VCF_DP_MAX=100           # Profondeur maximale (évite les régions dupliquées)
MAF_THRESHOLD=0.05       # Minor Allele Frequency ≥ 5% (variants communs)


# ==============================================================================
# SECTION 8 : ANNOTATION (SnpEff)
# Base de données osativa = Oryza sativa dans SnpEff
# ==============================================================================

SNPEFF_DB="osativa"


# ==============================================================================
# SECTION 9 : LOGGING
# Chaque exécution du pipeline crée un fichier log horodaté
# tee -a "$LOG_FILE" = affiche ET écrit dans le fichier simultanément
# ==============================================================================

# Note : LOGS_DIR doit exister avant d'appeler touch
# C'est pour ça que init_project() est appelé en premier dans 01_setup.sh
LOG_FILE="$LOGS_DIR/pipeline_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$LOGS_DIR"
touch "$LOG_FILE"


# ==============================================================================
# SECTION 10 : FONCTIONS UTILITAIRES
# Ces fonctions sont disponibles dans tous les scripts via source config.sh
# ==============================================================================

# ------------------------------------------------------------------------------
# log_message LEVEL "message"
# Affiche un message horodaté dans le terminal ET dans le fichier log
# Exemple : log_message "INFO" "Alignement terminé"
# ------------------------------------------------------------------------------
log_message() {
    local level="$1"
    shift
    local message="$*"
    local timestamp
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$timestamp] [$level] $message" | tee -a "$LOG_FILE"
}

# ------------------------------------------------------------------------------
# error_exit "message"
# Affiche une erreur et arrête immédiatement le script (exit 1)
# ------------------------------------------------------------------------------
error_exit() {
    log_message "ERROR" "$1"
    exit 1
}

# ------------------------------------------------------------------------------
# check_tool "nom_outil"
# Vérifie qu'un outil est accessible dans le PATH
# Retourne 0 (succès) si trouvé, 1 (échec) si absent
# ------------------------------------------------------------------------------
check_tool() {
    local tool="$1"
    if ! command -v "$tool" &> /dev/null; then
        log_message "WARN" "Outil introuvable : $tool"
        return 1
    fi
    return 0
}

# ------------------------------------------------------------------------------
# create_dir "/chemin/vers/dossier"
# Crée un dossier s'il n'existe pas encore (mkdir -p = crée les parents aussi)
# ------------------------------------------------------------------------------
create_dir() {
    local dir="$1"
    if [ ! -d "$dir" ]; then
        mkdir -p "$dir"
        log_message "INFO" "Dossier créé : $dir"
    fi
}

# ------------------------------------------------------------------------------
# check_disk_space MIN_GB
# Vérifie qu'il reste au moins MIN_GB Go disponibles sur la partition du projet
# ------------------------------------------------------------------------------
check_disk_space() {
    local min_gb="$1"
    local available_gb
    available_gb=$(df "$PROJECT_DIR" | awk 'NR==2 {print int($4/1024/1024)}')
    if [ "$available_gb" -lt "$min_gb" ]; then
        error_exit "Espace disque insuffisant. Requis : ${min_gb}Go, Disponible : ${available_gb}Go"
    fi
    log_message "INFO" "Espace disque OK : ${available_gb}Go disponibles"
}

# ------------------------------------------------------------------------------
# init_project
# Crée toute l'arborescence du projet d'un seul coup
# Appelé une seule fois dans 01_setup.sh
# ------------------------------------------------------------------------------
init_project() {
    log_message "INFO" "=== Initialisation de l'arborescence du projet ==="
    create_dir "$DATA_DIR"
    create_dir "$RESULTS_DIR"
    create_dir "$LOGS_DIR"
    create_dir "$SCRIPTS_DIR"
    create_dir "$REF_DIR"
    create_dir "$RAW_FASTQ_DIR"
    create_dir "$QC_DIR"
    create_dir "$TRIM_DIR"
    create_dir "$BAM_DIR"
    create_dir "$VCF_DIR"
    create_dir "$ANNO_DIR"
    create_dir "$ANALYSIS_DIR"
    # Créer le log maintenant que le dossier existe
    touch "$LOG_FILE"
    log_message "INFO" "Arborescence initialisée"
}

# ------------------------------------------------------------------------------
# print_config
# Affiche un résumé de la configuration courante (utile pour déboguer)
# ------------------------------------------------------------------------------
print_config() {
    cat << EOF
================================================================================
              RICE GENOMICS WGS PIPELINE - CONFIGURATION
================================================================================
Projet          : $PROJECT_DIR
Environnement   : $MAMBA_ENV
Référence       : $REF_FASTA
CPUs max        : $MAX_THREADS
RAM max         : $MAX_MEMORY
Espace disque   : $(df -h "$HOME" | awk 'NR==2 {print $4}') disponibles
================================================================================
EOF
}


# ==============================================================================
# EXPORTS
# Ces variables sont exportées pour être accessibles dans les sous-shells
# (utile si un script appelle un autre script ou une fonction externe)
# ==============================================================================

export PROJECT_DIR DATA_DIR RESULTS_DIR LOGS_DIR SCRIPTS_DIR
export REF_DIR RAW_FASTQ_DIR QC_DIR TRIM_DIR BAM_DIR VCF_DIR ANNO_DIR ANALYSIS_DIR
export FASTQ_DIR SAMPLES_FILE
export MAMBA_ENV CONDA_ACTIVATE
export MAX_THREADS MAX_MEMORY BWA_THREADS SAMTOOLS_THREADS BCFTOOLS_THREADS FASTP_THREADS
export REF_FASTA REF_GFF REF_INDEX REF_DICT
export VCF_QUAL_THRESHOLD VCF_DP_MIN VCF_DP_MAX MAF_THRESHOLD
export SNPEFF_DB LOG_FILE
