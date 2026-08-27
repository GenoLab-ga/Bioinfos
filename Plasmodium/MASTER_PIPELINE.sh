#!/bin/bash

# ====================================================================
# MASTER SCRIPT: Pipeline Single-Cell RNA-seq COMPLET
# Plasmodium falciparum Gamétocytes - PRJEB55754
# ====================================================================
# CORRECTION: SCRIPTS_DIR pointe vers le répertoire du script lui-même
#             (suppression du sous-dossier /scripts/ inexistant)
# ====================================================================

PROJECT_ROOT=~/Documents/github_projet/Plasmodium
SCRIPTS_DIR="$(cd "$(dirname "$0")" && pwd)"   # CORRIGÉ: était $(dirname "$0")/scripts
log_dir="$PROJECT_ROOT/logs"

mkdir -p "$log_dir"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_header()  { echo ""; echo -e "${BLUE}======================================================${NC}"; echo -e "${BLUE}$1${NC}"; echo -e "${BLUE}======================================================${NC}"; echo ""; }
print_success() { echo -e "${GREEN}✓ $1${NC}"; }
print_error()   { echo -e "${RED}❌ $1${NC}"; }
print_info()    { echo -e "${YELLOW}ℹ $1${NC}"; }

check_step_completed() { [ -f "$log_dir/step_${1}_completed" ]; }
mark_step_completed()  { touch "$log_dir/step_${1}_completed"; echo "[$(date)] Étape $1 complétée" >> "$log_dir/analysis.log"; }

# ====================================================================
# VÉRIFICATION DE L'ENVIRONNEMENT
# ====================================================================

check_environment() {
    print_header "Vérification de l'environnement"
    local ok=true
    for tool in fastqc multiqc fastp salmon kallisto python3 Rscript; do
        if command -v "$tool" &> /dev/null; then
            echo "  ✓ $tool"
        else
            echo "  ⚠ $tool non trouvé (optionnel selon étape)"
        fi
    done
    echo ""
}

# ====================================================================
# ÉTAPES
# ====================================================================

run_step_1() {
    print_header "ÉTAPE 1: FastQC"
    check_step_completed 1 && { print_info "Déjà complétée"; return 0; }
    local script="$SCRIPTS_DIR/01_fastqc_style.sh"
    [ ! -f "$script" ] && { print_error "Script introuvable: $script"; return 1; }
    bash "$script" && { print_success "FastQC OK"; mark_step_completed 1; } || { print_error "Erreur FastQC"; return 1; }
}

run_step_2() {
    print_header "ÉTAPE 2: Indexation Salmon"
    check_step_completed 2 && { print_info "Déjà complétée"; return 0; }
    local script="$SCRIPTS_DIR/02_salmon_indexing.sh"
    [ ! -f "$script" ] && { print_error "Script introuvable: $script"; return 1; }
    bash "$script" && { print_success "Index Salmon OK"; mark_step_completed 2; } || { print_error "Erreur indexation"; return 1; }
}

run_step_2k() {
    print_header "ÉTAPE 2k: Indexation Kallisto (optionnel)"
    check_step_completed 2k && { print_info "Déjà complétée"; return 0; }
    local script="$SCRIPTS_DIR/02_kallisto_indexing.sh"
    [ ! -f "$script" ] && { print_error "Script introuvable: $script"; return 1; }
    bash "$script" && { print_success "Index Kallisto OK"; mark_step_completed 2k; } || { print_error "Erreur indexation"; return 1; }
}

run_step_2b() {
    print_header "ÉTAPE 2b: Fastp Trimming"
    check_step_completed 2b && { print_info "Déjà complétée"; return 0; }
    command -v fastp &> /dev/null || { print_error "fastp non trouvé — mamba install fastp"; return 1; }
    local script="$SCRIPTS_DIR/02b_fastp_trimming.sh"
    [ ! -f "$script" ] && { print_error "Script introuvable: $script"; return 1; }
    bash "$script" && { print_success "Fastp trimming OK"; mark_step_completed 2b; } || { print_error "Erreur trimming"; return 1; }
}

run_step_3() {
    print_header "ÉTAPE 3: Quantification Salmon"
    check_step_completed 3 && { print_info "Déjà complétée"; return 0; }
    local script="$SCRIPTS_DIR/03_salmon_quantification.sh"
    [ ! -f "$script" ] && { print_error "Script introuvable: $script"; return 1; }
    bash "$script" && { print_success "Quantification Salmon OK"; mark_step_completed 3; } || { print_error "Erreur quantification"; return 1; }
}

run_step_4() {
    print_header "ÉTAPE 4: Conversion Transcripts → Gènes"
    check_step_completed 4 && { print_info "Déjà complétée"; return 0; }
    local script="$SCRIPTS_DIR/04_salmon_to_genecounts.py"
    [ ! -f "$script" ] && { print_error "Script introuvable: $script"; return 1; }
    python3 "$script" && { print_success "Conversion OK"; mark_step_completed 4; } || { print_error "Erreur conversion"; return 1; }
}

run_step_5() {
    print_header "ÉTAPE 5: Seurat QC"
    check_step_completed 5 && { print_info "Déjà complétée"; return 0; }
    local script="$SCRIPTS_DIR/05_seurat_qc.R"
    [ ! -f "$script" ] && { print_error "Script introuvable: $script"; return 1; }
    Rscript "$script" && { print_success "Seurat QC OK"; mark_step_completed 5; } || { print_error "Erreur Seurat QC"; return 1; }
}

# ====================================================================
# STATUS
# ====================================================================

show_status() {
    print_header "STATUS DU PIPELINE"
    for step in 1 2 2k 2b 3 4 5; do
        if check_step_completed $step; then
            print_success "Étape $step complétée"
        else
            echo "  ○ Étape $step: à faire"
        fi
    done
    echo ""
    echo "Fichiers importants:"
    echo "  MultiQC  : $PROJECT_ROOT/qc_reports/multiqc_report.html"
    echo "  Fastp    : $PROJECT_ROOT/qc_reports/fastp/"
    echo "  Matrice  : $PROJECT_ROOT/quantification/matrices/combined_counts_matrix.csv"
    echo "  Seurat   : $PROJECT_ROOT/analysis/01_qc/seurat_qc_filtered.rds"
    echo "  Logs     : $log_dir/analysis.log"

    # Vérifier que les scripts sont bien trouvés
    echo ""
    echo "Scripts dans: $SCRIPTS_DIR"
    for s in 01_fastqc_style.sh 02_salmon_indexing.sh 02_kallisto_indexing.sh \
             02b_fastp_trimming.sh 03_salmon_quantification.sh 03_kallisto_quantification.sh \
             04_salmon_to_genecounts.py 05_seurat_qc.R; do
        if [ -f "$SCRIPTS_DIR/$s" ]; then
            echo "  ✓ $s"
        else
            echo "  ❌ $s INTROUVABLE"
        fi
    done
}

# Réinitialiser une étape spécifique
reset_step() {
    local step=$1
    rm -f "$log_dir/step_${step}_completed"
    print_info "Étape $step réinitialisée"
}

# ====================================================================
# MAIN
# ====================================================================

print_header "PIPELINE SINGLE-CELL RNA-SEQ — Pf Gamétocytes"

STEP="${1:-all}"

case "$STEP" in
    env)    check_environment ;;
    status) show_status ;;
    reset)  reset_step "${2:-1}" ;;
    1)      run_step_1 ;;
    2)      run_step_1 && run_step_2 ;;
    2k)     run_step_1 && run_step_2k ;;
    2b)     run_step_1 && run_step_2 && run_step_2b ;;
    3)      run_step_1 && run_step_2 && run_step_2b && run_step_3 ;;
    4)      run_step_1 && run_step_2 && run_step_2b && run_step_3 && run_step_4 ;;
    5|all)
        check_environment
        run_step_1 && \
        run_step_2 && \
        run_step_2b && \
        run_step_3 && \
        run_step_4 && \
        run_step_5 && \
        {
            print_header "✓ PIPELINE TERMINÉ"
            echo "  Matrice : $PROJECT_ROOT/quantification/matrices/combined_counts_matrix.csv"
            echo "  Seurat  : $PROJECT_ROOT/analysis/01_qc/seurat_qc_filtered.rds"
            echo ""
        } || print_error "Pipeline interrompu"
        ;;
    *)
        echo "Usage: bash MASTER_PIPELINE.sh [commande]"
        echo ""
        echo "Commandes:"
        echo "  env     → Vérifier les outils installés"
        echo "  status  → Voir l'avancement + vérifier scripts"
        echo "  reset N → Réinitialiser l'étape N"
        echo "  1..5    → Exécuter jusqu'à l'étape N"
        echo "  all     → Tout exécuter"
        exit 1
        ;;
esac
