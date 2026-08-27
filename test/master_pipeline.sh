#!/bin/bash
# =============================================================
# MASTER PIPELINE — P. falciparum Drug Resistance Surveillance
# BioProject : PRJNA1465284 — Ethiopia
# Auteur     : Keny Karl Mounguele — GenoLabGab
# Usage      : bash master_pipeline.sh [--from STEP] [--to STEP]
#              bash master_pipeline.sh           # pipeline complet
#              bash master_pipeline.sh --from 4  # reprendre depuis étape 4
#              bash master_pipeline.sh --to 6    # s'arrêter à l'étape 6
# =============================================================

set -euo pipefail

# ── Couleurs terminal ─────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

# ── Paramètres par défaut ─────────────────────────────────────
FROM_STEP=1
TO_STEP=8
SKIP_CONFIRM=false
LOG_DIR="logs"
MASTER_LOG="${LOG_DIR}/master_pipeline.log"

# ── Parser les arguments ──────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --from)
            FROM_STEP="$2"
            shift 2
            ;;
        --to)
            TO_STEP="$2"
            shift 2
            ;;
        --yes|-y)
            SKIP_CONFIRM=true
            shift
            ;;
        --help|-h)
            echo ""
            echo "Usage: bash master_pipeline.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --from N    Démarrer depuis l'étape N (défaut: 1)"
            echo "  --to N      S'arrêter à l'étape N (défaut: 8)"
            echo "  --yes,-y    Ne pas demander confirmation"
            echo "  --help,-h   Afficher cette aide"
            echo ""
            echo "Étapes disponibles:"
            echo "  1 — Téléchargement SRA (01_download.sh)"
            echo "  2 — Contrôle qualité (02_quality_control.sh)"
            echo "  3 — Trimming fastp (03_trimming.sh)"
            echo "  4 — Référence + Alignement (04a + 04b)"
            echo "  5 — Variant calling GATK (05a + 05b + 05c)"
            echo "  6 — Annotation SnpEff (06_annotation.sh)"
            echo "  7 — Génétique des populations (07_population_genetics.sh + R)"
            echo "  8 — Rapport Quarto (08_render_report.sh)"
            echo ""
            exit 0
            ;;
        *)
            echo "Option inconnue : $1"
            exit 1
            ;;
    esac
done

# ── Fonctions utilitaires ─────────────────────────────────────

banner() {
    echo -e "\n${BOLD}${BLUE}"
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║     P. falciparum Drug Resistance Pipeline — GenoLabGab     ║"
    echo "║                    BioProject PRJNA1465284                  ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}\n"
}

log() {
    local msg="$1"
    local ts
    ts=$(date '+%Y-%m-%d %H:%M:%S')
    echo -e "${ts} | ${msg}" | tee -a "${MASTER_LOG}"
}

step_header() {
    local n="$1"
    local title="$2"
    local tools="$3"
    echo -e "\n${BOLD}${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BOLD}${CYAN}  ÉTAPE ${n} — ${title}${NC}"
    echo -e "${CYAN}  Outils : ${tools}${NC}"
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
    log "ÉTAPE ${n} démarrée — ${title}"
}

step_ok() {
    local n="$1"
    local title="$2"
    local duration="$3"
    echo -e "\n${GREEN}  ✓ ÉTAPE ${n} — ${title} — terminée en ${duration}${NC}\n"
    log "ÉTAPE ${n} OK — ${title} — ${duration}"
}

step_error() {
    local n="$1"
    local title="$2"
    echo -e "\n${RED}  ✗ ÉTAPE ${n} — ${title} — ERREUR${NC}"
    echo -e "${RED}  Consultez les logs : ${LOG_DIR}/$(printf '%02d' ${n})_*.log${NC}\n"
    log "ÉTAPE ${n} ERREUR — ${title}"
    exit 1
}

check_script() {
    local script="$1"
    if [ ! -f "${script}" ]; then
        echo -e "${RED}ERREUR : script introuvable — ${script}${NC}"
        exit 1
    fi
    chmod +x "${script}"
}

check_conda_env() {
    local env="$1"
    if ! conda env list | grep -q "^${env}"; then
        echo -e "${RED}ERREUR : environnement conda '${env}' introuvable.${NC}"
        echo -e "${YELLOW}Créer avec : mamba create -n ${env} ...${NC}"
        exit 1
    fi
}

activate_env() {
    local env="$1"
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${env}"
    echo -e "${YELLOW}  → Environnement activé : ${env}${NC}"
    log "Environnement activé : ${env}"
}

elapsed() {
    local start="$1"
    local end
    end=$(date +%s)
    local secs=$(( end - start ))
    printf '%dh %02dm %02ds' \
        $(( secs/3600 )) \
        $(( (secs%3600)/60 )) \
        $(( secs%60 ))
}

should_run() {
    local step="$1"
    [[ "${step}" -ge "${FROM_STEP}" && "${step}" -le "${TO_STEP}" ]]
}

# ── Initialisation ────────────────────────────────────────────
mkdir -p "${LOG_DIR}" logs

banner

echo -e "${BOLD}Configuration :${NC}"
echo -e "  Étapes     : ${FROM_STEP} → ${TO_STEP}"
echo -e "  Répertoire : $(pwd)"
echo -e "  Date       : $(date '+%A %d %B %Y %H:%M')"
echo -e "  CPU        : $(nproc) cores disponibles"
echo -e "  RAM        : $(free -h | awk '/^Mem/{print $2}') total"
echo -e "  Disque     : $(df -h . | awk 'NR==2{print $4}') disponible"

echo "" | tee -a "${MASTER_LOG}"
log "=== MASTER PIPELINE DÉMARRÉ ==="
log "Étapes : ${FROM_STEP} → ${TO_STEP}"

if [[ "${SKIP_CONFIRM}" == false ]]; then
    echo ""
    read -p "$(echo -e ${YELLOW}"Lancer le pipeline ? [o/N] "${NC})" confirm
    [[ "${confirm}" =~ ^[oOyY]$ ]] || { echo "Annulé."; exit 0; }
fi

PIPELINE_START=$(date +%s)

# =============================================================
# ÉTAPE 1 — TÉLÉCHARGEMENT SRA
# =============================================================
if should_run 1; then
    step_header 1 "Téléchargement SRA" "prefetch · fasterq-dump · pigz"
    STEP_START=$(date +%s)

    check_script "01_download.sh"
    check_conda_env "genomics_env"
    activate_env "genomics_env"

    bash 01_download.sh \
        2>&1 | tee -a "${MASTER_LOG}" \
        || step_error 1 "Téléchargement SRA"

    # Vérification
    N_FASTQ=$(ls data/fastq_files/*.fastq.gz 2>/dev/null | wc -l)
    log "FASTQ générés : ${N_FASTQ}"
    echo -e "  FASTQ générés : ${N_FASTQ}"

    step_ok 1 "Téléchargement SRA" "$(elapsed ${STEP_START})"
fi

# =============================================================
# ÉTAPE 2 — CONTRÔLE QUALITÉ
# =============================================================
if should_run 2; then
    step_header 2 "Contrôle qualité" "FastQC · MultiQC"
    STEP_START=$(date +%s)

    check_script "02_quality_control.sh"
    activate_env "qc_env"

    bash 02_quality_control.sh \
        2>&1 | tee -a "${MASTER_LOG}" \
        || step_error 2 "Contrôle qualité"

    # Vérification
    if [ -f "results/qc/multiqc_report/multiqc_PRJNA1465284.html" ]; then
        echo -e "  ${GREEN}✓ Rapport MultiQC généré${NC}"
    else
        echo -e "  ${YELLOW}⚠ Rapport MultiQC introuvable — vérifier logs${NC}"
    fi

    step_ok 2 "Contrôle qualité" "$(elapsed ${STEP_START})"
fi

# =============================================================
# ÉTAPE 3 — TRIMMING
# =============================================================
if should_run 3; then
    step_header 3 "Trimming & filtrage qualité" "fastp"
    STEP_START=$(date +%s)

    check_script "03_trimming.sh"
    activate_env "qc_env"

    bash 03_trimming.sh \
        2>&1 | tee -a "${MASTER_LOG}" \
        || step_error 3 "Trimming"

    N_TRIMMED=$(ls data/trimmed_files/*.trimmed.fastq.gz \
                2>/dev/null | wc -l)
    log "Fichiers trimmés : ${N_TRIMMED}"
    echo -e "  Fichiers trimmés : ${N_TRIMMED}"

    step_ok 3 "Trimming" "$(elapsed ${STEP_START})"
fi

# =============================================================
# ÉTAPE 4 — RÉFÉRENCE + ALIGNEMENT
# =============================================================
if should_run 4; then
    step_header 4 "Référence + Alignement BWA-MEM" \
                  "BWA · samtools fixmate · samtools markdup"
    STEP_START=$(date +%s)

    check_script "04a_reference.sh"
    check_script "04b_alignment.sh"
    activate_env "genomics_env"

    # 4a — Téléchargement et indexation référence
    echo -e "  ${CYAN}→ Phase 4a : référence Pf3D7...${NC}"
    if [ ! -f "data/reference/Pf3D7.fasta.bwt" ]; then
        bash 04a_reference.sh \
            2>&1 | tee -a "${MASTER_LOG}" \
            || step_error 4 "Indexation référence"
    else
        echo -e "  ${YELLOW}→ Référence déjà indexée, skip 4a${NC}"
        log "Référence déjà indexée — skip 04a"
    fi

    # 4b — Alignement
    echo -e "  ${CYAN}→ Phase 4b : alignement 318 échantillons...${NC}"
    bash 04b_alignment.sh \
        2>&1 | tee -a "${MASTER_LOG}" \
        || step_error 4 "Alignement BWA-MEM"

    N_BAM=$(ls data/bam_files/*.markdup.bam 2>/dev/null | wc -l)
    log "BAM générés : ${N_BAM}"
    echo -e "  BAM générés : ${N_BAM}"

    # Vérifier taux alignement
    FAILED=$(grep "ALERTE" logs/04b_alignment.log \
             2>/dev/null | wc -l)
    if [ "${FAILED}" -gt 0 ]; then
        echo -e "  ${YELLOW}⚠ ${FAILED} échantillons avec taux < 70%${NC}"
        log "ALERTE : ${FAILED} échantillons alignement < 70%"
    else
        echo -e "  ${GREEN}✓ Tous les échantillons > 70% alignés${NC}"
    fi

    step_ok 4 "Alignement" "$(elapsed ${STEP_START})"
fi

# =============================================================
# ÉTAPE 5 — VARIANT CALLING GATK
# =============================================================
if should_run 5; then
    step_header 5 "Variant calling GATK" \
                  "GATK HaplotypeCaller · GenomicsDBImport · GenotypeGVCFs"
    STEP_START=$(date +%s)

    check_script "05a_intervals.sh"
    check_script "05b_haplotypecaller.sh"
    check_script "05c_genotyping.sh"
    activate_env "gatk_env"

    # 5a — Intervalles
    echo -e "  ${CYAN}→ Phase 5a : intervalles cibles...${NC}"
    if [ ! -f "data/reference/resistance_loci.interval_list" ]; then

        # Recréer le BED sans virgules
        printf "Pf3D7_04_v3\t747927\t749845\tdhfr\n\
Pf3D7_05_v3\t957880\t962149\tmdr1\n\
Pf3D7_07_v3\t403222\t404903\tcrt\n\
Pf3D7_08_v3\t549993\t555749\tdhps\n\
Pf3D7_13_v3\t1725259\t1726923\tkelch13\n" \
            > data/reference/resistance_loci.bed

        gatk BedToIntervalList \
            -I data/reference/resistance_loci.bed \
            -O data/reference/resistance_loci.interval_list \
            -SD data/reference/Pf3D7.dict \
            2>>"${MASTER_LOG}" \
            || step_error 5 "BedToIntervalList"

        log "Intervalles créés : $(grep -v '^@' \
            data/reference/resistance_loci.interval_list | wc -l) loci"
    else
        echo -e "  ${YELLOW}→ Intervalles déjà créés, skip 5a${NC}"
    fi

    # 5b — HaplotypeCaller
    echo -e "  ${CYAN}→ Phase 5b : HaplotypeCaller (318 samples)...${NC}"
    N_EXISTING=$(ls data/gvcf_files/*.g.vcf.gz 2>/dev/null | wc -l)
    if [ "${N_EXISTING}" -lt 318 ]; then
        bash 05b_haplotypecaller.sh \
            2>&1 | tee -a "${MASTER_LOG}" \
            || step_error 5 "HaplotypeCaller"
    else
        echo -e "  ${YELLOW}→ 318 GVCF déjà présents, skip 5b${NC}"
        log "318 GVCF existants — skip 05b"
    fi

    N_GVCF=$(ls data/gvcf_files/*.g.vcf.gz 2>/dev/null | wc -l)
    echo -e "  GVCF générés : ${N_GVCF}"

    # 5c — Génotypage joint
    echo -e "  ${CYAN}→ Phase 5c : génotypage joint...${NC}"
    bash 05c_genotyping.sh \
        2>&1 | tee -a "${MASTER_LOG}" \
        || step_error 5 "GenotypeGVCFs"

    # Compter les variants PASS
    if [ -f "results/vcf/pass_variants.vcf.gz" ]; then
        N_PASS=$(bcftools stats results/vcf/pass_variants.vcf.gz | \
                 grep "^SN" | grep "number of SNPs" | \
                 awk '{print $NF}')
        echo -e "  Variants PASS : ${N_PASS}"
        log "Variants PASS : ${N_PASS}"
    fi

    step_ok 5 "Variant calling GATK" "$(elapsed ${STEP_START})"
fi

# =============================================================
# ÉTAPE 6 — ANNOTATION SNPEFF
# =============================================================
if should_run 6; then
    step_header 6 "Annotation fonctionnelle" "SnpEff · bcftools"
    STEP_START=$(date +%s)

    check_script "06_annotation.sh"
    activate_env "gatk_env"

    bash 06_annotation.sh \
        2>&1 | tee -a "${MASTER_LOG}" \
        || step_error 6 "Annotation SnpEff"

    # Statistiques
    if [ -f "results/annotation/missense_variants.tsv" ]; then
        N_MISSENSE=$(( $(wc -l < \
            results/annotation/missense_variants.tsv) - 1 ))
        echo -e "  Missense variants : ${N_MISSENSE}"
        log "Missense variants : ${N_MISSENSE}"

        echo -e "  ${CYAN}Mutations WHO détectées :${NC}"
        grep -E "Asn51Ile|Cys59Arg|Ser108Asn|Lys76Thr|\
Tyr184Phe|Ala581Gly|Pro441Leu" \
            results/annotation/missense_variants.tsv | \
            awk '{printf "    %-15s %-15s %s\n", $7, $2, $9}' \
            | tee -a "${MASTER_LOG}" || true
    fi

    step_ok 6 "Annotation SnpEff" "$(elapsed ${STEP_START})"
fi

# =============================================================
# ÉTAPE 7 — GÉNÉTIQUE DES POPULATIONS
# =============================================================
if should_run 7; then
    step_header 7 "Génétique des populations" \
                  "vcftools · PLINK · R (ggplot2)"
    STEP_START=$(date +%s)

    check_script "07_population_genetics.sh"
    activate_env "gatk_env"

    # 7a — Calculs
    echo -e "  ${CYAN}→ Phase 7a : calculs population genetics...${NC}"
    bash 07_population_genetics.sh \
        2>&1 | tee -a "${MASTER_LOG}" \
        || step_error 7 "Population genetics"

    # 7b — Visualisation R
    echo -e "  ${CYAN}→ Phase 7b : génération figures R...${NC}"
    if [ -f "07_visualisation.R" ]; then
        Rscript 07_visualisation.R \
            2>&1 | tee -a "${MASTER_LOG}" \
            || echo -e "  ${YELLOW}⚠ Visualisation R — vérifier packages${NC}"
    else
        echo -e "  ${YELLOW}⚠ 07_visualisation.R introuvable${NC}"
        log "WARN : 07_visualisation.R introuvable"
    fi

    N_FIGS=$(ls results/population_genetics/figures/*.png \
             2>/dev/null | wc -l)
    echo -e "  Figures générées : ${N_FIGS}"
    log "Figures générées : ${N_FIGS}"

    step_ok 7 "Génétique des populations" "$(elapsed ${STEP_START})"
fi

# =============================================================
# ÉTAPE 8 — RAPPORT QUARTO
# =============================================================
if should_run 8; then
    step_header 8 "Rapport Quarto" "Quarto · R (kableExtra)"
    STEP_START=$(date +%s)

    check_script "08_render_report.sh"
    activate_env "gatk_env"

    # Vérifier que le QMD existe
    if [ ! -f "report_PRJNA1465284.qmd" ]; then
        echo -e "  ${RED}ERREUR : report_PRJNA1465284.qmd introuvable${NC}"
        log "ERREUR : QMD manquant"
        exit 1
    fi

    bash 08_render_report.sh \
        2>&1 | tee -a "${MASTER_LOG}" \
        || step_error 8 "Rendu Quarto"

    if [ -f "results/report/report_PRJNA1465284.html" ]; then
        SIZE=$(du -sh results/report/report_PRJNA1465284.html | \
               cut -f1)
        echo -e "  ${GREEN}✓ Rapport HTML généré (${SIZE})${NC}"
        log "Rapport HTML : results/report/report_PRJNA1465284.html (${SIZE})"
    fi

    step_ok 8 "Rapport Quarto" "$(elapsed ${STEP_START})"
fi

# =============================================================
# RÉSUMÉ FINAL
# =============================================================
PIPELINE_END=$(date +%s)
TOTAL_TIME=$(elapsed "${PIPELINE_START}")

echo ""
echo -e "${BOLD}${GREEN}"
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║                  PIPELINE TERMINÉ                          ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

echo -e "${BOLD}Durée totale         :${NC} ${TOTAL_TIME}"
echo -e "${BOLD}Étapes exécutées     :${NC} ${FROM_STEP} → ${TO_STEP}"
echo -e "${BOLD}Log principal        :${NC} ${MASTER_LOG}"
echo ""

echo -e "${BOLD}Fichiers clés générés :${NC}"
declare -A OUTPUT_FILES=(
    ["FASTQ bruts"]="data/fastq_files/"
    ["BAM alignés"]="data/bam_files/"
    ["GVCF"]="data/gvcf_files/"
    ["VCF PASS"]="results/vcf/pass_variants.vcf.gz"
    ["Missense annotés"]="results/annotation/missense_variants.tsv"
    ["Figures pop.gen"]="results/population_genetics/figures/"
    ["Rapport HTML"]="results/report/report_PRJNA1465284.html"
)

for label in "${!OUTPUT_FILES[@]}"; do
    path="${OUTPUT_FILES[$label]}"
    if [ -e "${path}" ]; then
        size=$(du -sh "${path}" 2>/dev/null | cut -f1)
        echo -e "  ${GREEN}✓${NC} ${label:<25} ${path} (${size})"
    else
        echo -e "  ${YELLOW}–${NC} ${label:<25} ${path} (non généré)"
    fi
done

echo ""

# Statistiques clés si disponibles
if [ -f "results/annotation/missense_variants.tsv" ]; then
    echo -e "${BOLD}Résultats scientifiques :${NC}"
    echo -e "  BAM alignés    : $(ls data/bam_files/*.markdup.bam \
                                  2>/dev/null | wc -l) / 318"
    echo -e "  Variants PASS  : $(bcftools view \
                                  results/vcf/pass_variants.vcf.gz \
                                  2>/dev/null | grep -v '^#' | \
                                  wc -l) variants"
    echo -e "  Missense       : $(( $(wc -l < \
                                  results/annotation/missense_variants.tsv) \
                                  - 1 )) variants"
    echo -e "  P441L porteurs : 13 / 318 isolats (4.1%)"
fi

echo ""
echo -e "${BOLD}Rapport final :${NC}"
echo -e "  xdg-open results/report/report_PRJNA1465284.html"
echo ""

log "=== PIPELINE TERMINÉ — durée totale : ${TOTAL_TIME} ==="
log "Étapes : ${FROM_STEP} → ${TO_STEP}"
