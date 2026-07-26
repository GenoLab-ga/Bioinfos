#!/bin/bash
################################################################################
# RICE GENOMICS WGS PIPELINE — 05_variant_calling.sh
# Détection des variants avec bcftools mpileup/call
#
# EXÉCUTION : bash 05_variant_calling.sh
# DURÉE     : ~30-40 min pour 8 samples
#
# PIPELINE :
#   1. bcftools mpileup  — calcul du pileup multi-samples
#   2. bcftools call -m  — variant calling (modèle multi-allélique)
#   3. Filtre QUAL/DP    — qualité et profondeur
#   4. Filtre MAF ≥ 5%   — variants communs seulement
#
# SORTIES :
#   results/04_vcf/all_samples.raw.vcf.gz      — variants bruts
#   results/04_vcf/all_samples.filtered.vcf.gz — après QUAL/DP
#   results/04_vcf/all_samples.maf0.05.vcf.gz  — final (MAF ≥ 5%)
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

log_message "INFO" "=== ÉTAPE 5 : VARIANT CALLING ==="
eval "$CONDA_ACTIVATE"

BAM_LIST="$BAM_DIR/bam_list.txt"
[ -f "$BAM_LIST" ] && [ -s "$BAM_LIST" ] || error_exit "bam_list.txt vide — lancer 04_mapping.sh"
BAM_COUNT=$(wc -l < "$BAM_LIST")
log_message "INFO" "$BAM_COUNT fichiers BAM"

mkdir -p "$VCF_DIR"

# --- Pileup multi-samples ---
PILEUP="$VCF_DIR/all_samples.pileup.bcf"
if [ ! -f "$PILEUP" ]; then
    log_message "INFO" "Calcul du pileup (étape la plus longue ~25 min)..."
    bcftools mpileup \
        -f "$REF_FASTA" \
        -b "$BAM_LIST" \
        -q 30 -Q 20 -d 500 \
        --threads "$BCFTOOLS_THREADS" \
        -o "$PILEUP" -O b \
        2>> "$LOG_FILE"
    log_message "INFO" "✓ Pileup : $PILEUP"
else
    log_message "INFO" "✓ Pileup déjà présent"
fi

# --- Variant calling ---
VCF_RAW="$VCF_DIR/all_samples.raw.vcf.gz"
if [ ! -f "$VCF_RAW" ]; then
    log_message "INFO" "Appel des variants..."
    bcftools call -m -v -o "$VCF_RAW" -O z \
        --threads "$BCFTOOLS_THREADS" "$PILEUP" 2>> "$LOG_FILE"
    bcftools index "$VCF_RAW"
    N_RAW=$(bcftools view -H "$VCF_RAW" 2>/dev/null | wc -l)
    log_message "INFO" "✓ VCF brut : $N_RAW variants"
else
    log_message "INFO" "✓ VCF brut déjà présent"
fi

# --- Filtre QUAL et profondeur ---
VCF_FILTERED="$VCF_DIR/all_samples.filtered.vcf.gz"
if [ ! -f "$VCF_FILTERED" ]; then
    log_message "INFO" "Filtrage QUAL≥$VCF_QUAL_THRESHOLD | DP $VCF_DP_MIN-$VCF_DP_MAX..."
    bcftools filter \
        -i "QUAL>=$VCF_QUAL_THRESHOLD & INFO/DP>=$VCF_DP_MIN & INFO/DP<=$VCF_DP_MAX" \
        -o "$VCF_FILTERED" -O z "$VCF_RAW" 2>> "$LOG_FILE"
    bcftools index "$VCF_FILTERED"
    N_FILT=$(bcftools view -H "$VCF_FILTERED" 2>/dev/null | wc -l)
    log_message "INFO" "✓ VCF filtré : $N_FILT variants"
else
    log_message "INFO" "✓ VCF filtré déjà présent"
fi

# --- Filtre MAF ---
VCF_MAF="$VCF_DIR/all_samples.maf${MAF_THRESHOLD}.vcf.gz"
if [ ! -f "$VCF_MAF" ] || [ ! -s "$VCF_MAF" ] || [ "$(wc -c < "$VCF_MAF")" -lt 1000 ]; then
    log_message "INFO" "Filtrage MAF ≥ $MAF_THRESHOLD..."
    bcftools view --threads "$BCFTOOLS_THREADS" \
        -q "${MAF_THRESHOLD}:minor" \
        -o "$VCF_MAF" -O z "$VCF_FILTERED" 2>> "$LOG_FILE"
    bcftools index "$VCF_MAF"
    N_MAF=$(bcftools view -H "$VCF_MAF" 2>/dev/null | wc -l)
    log_message "INFO" "✓ VCF final (MAF≥$MAF_THRESHOLD) : $N_MAF variants"
else
    log_message "INFO" "✓ VCF MAF déjà présent"
fi

# --- Statistiques ---
cat > "$VCF_DIR/variant_stats.txt" << STATS
=== STATISTIQUES VARIANT CALLING ===
Référence    : Oryza sativa IRGSP v1.0
Appeleur     : bcftools mpileup/call
Paramètres   : QUAL≥${VCF_QUAL_THRESHOLD}, DP ${VCF_DP_MIN}-${VCF_DP_MAX}, MAF≥${MAF_THRESHOLD}

VCF brut     : $(bcftools view -H "$VCF_RAW" 2>/dev/null | wc -l) variants
VCF filtré   : $(bcftools view -H "$VCF_FILTERED" 2>/dev/null | wc -l) variants
VCF final    : $(bcftools view -H "$VCF_MAF" 2>/dev/null | wc -l) variants
Samples      : $(bcftools query -l "$VCF_MAF" 2>/dev/null | tr '\n' ' ')
STATS

log_message "INFO" "=== VARIANT CALLING TERMINÉ ==="
log_message "INFO" "VCF final      : $VCF_MAF"
log_message "INFO" "Étape suivante : bash 06_annotation.sh"
