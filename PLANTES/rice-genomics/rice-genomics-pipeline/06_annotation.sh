#!/bin/bash
################################################################################
# RICE GENOMICS WGS PIPELINE — 06_annotation.sh
# Annotation fonctionnelle des variants avec SnpEff
#
# EXÉCUTION : bash 06_annotation.sh
# DURÉE     : ~5 min
#
# PRÉREQUIS :
#   - snpeff 5.1 (compatible Java 17) : mamba install -c bioconda snpeff=5.1
#   - Base Oryza_sativa téléchargée lors du setup
#
# SORTIES :
#   results/05_annotation/all_samples.annotated.vcf.gz — VCF annoté
#   results/05_annotation/high_impact_variants.tsv     — variants HIGH impact
#   results/05_annotation/snpEff_summary.html          — rapport SnpEff
#   results/05_annotation/snpeff_stats.csv             — stats CSV
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

log_message "INFO" "=== ÉTAPE 6 : ANNOTATION SNPEFF ==="
eval "$CONDA_ACTIVATE"

check_tool "snpEff" || error_exit "snpEff introuvable — mamba install -c bioconda snpeff=5.1"

VCF_MAF="$VCF_DIR/all_samples.maf${MAF_THRESHOLD}.vcf.gz"
[ -f "$VCF_MAF" ] && [ -s "$VCF_MAF" ] || error_exit "VCF introuvable — lancer 05_variant_calling.sh"

mkdir -p "$ANNO_DIR"

VCF_ANNOTATED="$ANNO_DIR/all_samples.annotated.vcf.gz"

if [ ! -f "$VCF_ANNOTATED" ] || [ "$(wc -c < "$VCF_ANNOTATED")" -lt 10000 ]; then

    # Décompresser le VCF avec bcftools (le fichier est en BGZF, pas gzip standard)
    log_message "INFO" "Décompression du VCF..."
    bcftools view "$VCF_MAF" -O v -o "$ANNO_DIR/input.vcf" \
        --threads "$BCFTOOLS_THREADS" 2>> "$LOG_FILE"
    N_INPUT=$(grep -v "^#" "$ANNO_DIR/input.vcf" | wc -l)
    log_message "INFO" "✓ VCF décompressé : $N_INPUT variants"

    # Annoter avec SnpEff
    log_message "INFO" "Annotation SnpEff ($SNPEFF_DB)..."
    snpEff ann \
        -config "$SNPEFF_CFG" \
        -stats "$ANNO_DIR/snpEff_summary.html" \
        -csvStats "$ANNO_DIR/snpeff_stats.csv" \
        -v "$SNPEFF_DB" \
        "$ANNO_DIR/input.vcf" \
        > "$ANNO_DIR/all_samples.annotated.vcf" \
        2> "$ANNO_DIR/snpeff.log"

    N_ANNO=$(grep -v "^#" "$ANNO_DIR/all_samples.annotated.vcf" | wc -l)
    log_message "INFO" "✓ Variants annotés : $N_ANNO"

    # Compresser et indexer
    log_message "INFO" "Compression et indexation..."
    bgzip "$ANNO_DIR/all_samples.annotated.vcf"
    bcftools index "$VCF_ANNOTATED"

    # Nettoyer le fichier temporaire
    rm -f "$ANNO_DIR/input.vcf"
    log_message "INFO" "✓ VCF annoté : $VCF_ANNOTATED"
else
    log_message "INFO" "✓ VCF annoté déjà présent"
fi

# --- Distribution par impact ---
log_message "INFO" "Calcul de la distribution par impact..."
echo "=== DISTRIBUTION PAR IMPACT ===" > "$ANNO_DIR/impact_summary.txt"
bcftools query -f '%INFO/ANN\n' "$VCF_ANNOTATED" 2>/dev/null | \
    grep -oP '\|HIGH\||\|MODERATE\||\|LOW\||\|MODIFIER\|' | \
    sort | uniq -c | sort -rn >> "$ANNO_DIR/impact_summary.txt"
cat "$ANNO_DIR/impact_summary.txt"

# --- Extraction des variants HIGH impact ---
log_message "INFO" "Extraction des variants HIGH impact..."
bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANN\n' \
    "$VCF_ANNOTATED" 2>/dev/null | \
    grep "|HIGH|" | \
    awk -F'\t' '{
        split($5, ann, ",");
        for(i in ann) {
            split(ann[i], f, "|");
            if(f[3]=="HIGH")
                print $1"\t"$2"\t"$3"\t"$4"\t"f[2]"\t"f[4]"\t"f[10]
        }
    }' > "$ANNO_DIR/high_impact_variants.tsv"

N_HIGH=$(wc -l < "$ANNO_DIR/high_impact_variants.tsv")
log_message "INFO" "✓ Variants HIGH impact : $N_HIGH → $ANNO_DIR/high_impact_variants.tsv"

log_message "INFO" "=== ANNOTATION TERMINÉE ==="
log_message "INFO" "Rapport SnpEff : $ANNO_DIR/snpEff_summary.html"
log_message "INFO" "Étape suivante : bash 07_analysis.sh"
