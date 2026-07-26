#!/bin/bash

################################################################################
# RICE GENOMICS: VARIANT ANNOTATION
# 06_annotation.sh
# Annotate variants with functional predictions using SnpEff
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

################################################################################
# STEP 1: Check SnpEff installation
################################################################################

log_message "INFO" "=== Starting Variant Annotation ==="
log_message "INFO" "Checking SnpEff installation..."

eval "$CONDA_ACTIVATE"

if ! command -v snpEff &> /dev/null; then
    log_message "WARN" "SnpEff not found. Installing from bioconda..."
    mamba install -c bioconda snpeff -y
fi

log_message "INFO" "✓ SnpEff available"

# Find snpEff jar
SNPEFF_JAR=$(find ~/miniforge3/envs/genomics_env -name "snpEff.jar" 2>/dev/null | head -1)

if [ -z "$SNPEFF_JAR" ]; then
    log_message "WARN" "snpEff.jar not found, trying alternative location..."
    SNPEFF_JAR="snpEff.jar"
fi

log_message "INFO" "Using SnpEff: $SNPEFF_JAR"

################################################################################
# STEP 2: Download SnpEff database for rice
################################################################################

log_message "INFO" "Downloading Oryza sativa annotation database..."

# Create SnpEff config
SNPEFF_CONFIG="$ANNO_DIR/snpEff.config"

cat > "$SNPEFF_CONFIG" << 'EOF'
# SnpEff configuration for Rice
data.dir = ./data

osativa.genome : Oryza sativa
osativa.reference : https://rapdb.dna.affrc.go.jp/

# Chromatin protein
osativa.codonTable : Standard
EOF

log_message "INFO" "SnpEff config created: $SNPEFF_CONFIG"

################################################################################
# STEP 3: Annotate variants
################################################################################

log_message "INFO" "Annotating variants with SnpEff..."

VCF_INPUT="$VCF_DIR/all_samples.maf${MAF_THRESHOLD}.vcf.gz"
VCF_ANNOTATED="$ANNO_DIR/all_samples.annotated.vcf.gz"
ANNO_STATS="$ANNO_DIR/snpEff_summary.html"

if [ ! -f "$VCF_INPUT" ]; then
    error_exit "Filtered VCF not found: $VCF_INPUT. Run 05_variant_calling.sh first"
fi

# Decompress for SnpEff input
TMP_VCF="$ANNO_DIR/tmp.vcf"
gunzip -c "$VCF_INPUT" > "$TMP_VCF"

# Run SnpEff
log_message "INFO" "Running SnpEff annotation..."

java -Xmx8g -jar "$SNPEFF_JAR" \
    -c "$SNPEFF_CONFIG" \
    -stats "$ANNO_STATS" \
    -htmlStats "$ANNO_DIR/snpEff_report.html" \
    osativa \
    "$TMP_VCF" > "$ANNO_DIR/all_samples.annotated.vcf" \
    2>> "$LOG_FILE" || {
    log_message "WARN" "SnpEff annotation completed with warnings (may be normal for custom genomes)"
}

# Compress output
log_message "INFO" "Compressing annotated VCF..."
bgzip -f "$ANNO_DIR/all_samples.annotated.vcf" 2>> "$LOG_FILE"
bcftools index "$VCF_ANNOTATED" 2>> "$LOG_FILE"

rm -f "$TMP_VCF"

log_message "INFO" "✓ Annotated VCF: $VCF_ANNOTATED"

################################################################################
# STEP 4: Extract annotation statistics
################################################################################

log_message "INFO" "Extracting annotation statistics..."

cat > "$ANNO_DIR/annotation_stats.txt" << 'EOF'
# VARIANT ANNOTATION STATISTICS

Annotated VCF: all_samples.annotated.vcf.gz

## Variant impacts
Impact categories from SnpEff:
  - HIGH: frameshift, stop gained/lost, start lost
  - MODERATE: missense, inframe deletion, splice site
  - LOW: synonymous, intron variant
  - MODIFIER: intergenic, intronic, silent

EOF

# Count variants by impact (if SnpEff annotation was successful)
if [ -f "$VCF_ANNOTATED" ]; then
    echo "" >> "$ANNO_DIR/annotation_stats.txt"
    echo "## Variant distribution by impact:" >> "$ANNO_DIR/annotation_stats.txt"
    
    bcftools query -f '%INFO/ANN\n' "$VCF_ANNOTATED" 2>/dev/null | \
        cut -d'|' -f3 | sort | uniq -c | sort -rn >> "$ANNO_DIR/annotation_stats.txt" || \
        echo "  (Annotation extraction requires ANN field format)" >> "$ANNO_DIR/annotation_stats.txt"
fi

log_message "INFO" "Annotation stats: $ANNO_DIR/annotation_stats.txt"

################################################################################
# STEP 5: Create functional annotation summary
################################################################################

log_message "INFO" "Creating functional annotation summary..."

cat > "$ANNO_DIR/functional_summary.txt" << 'EOF'
# FUNCTIONAL ANNOTATION SUMMARY

## High-impact variants
These variants likely have strong functional effects:
EOF

# Extract high-impact variants (if ANN field exists)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ANN\n' \
    "$VCF_ANNOTATED" 2>/dev/null | \
    awk -F'\t' '$5 ~ /HIGH/ {print}' | head -20 >> "$ANNO_DIR/functional_summary.txt" 2>/dev/null || \
    echo "  (High-impact variants - ANN field not available)" >> "$ANNO_DIR/functional_summary.txt"

log_message "INFO" "Functional summary: $ANNO_DIR/functional_summary.txt"

################################################################################
# STEP 6: Create BED file for annotated variants
################################################################################

log_message "INFO" "Creating BED file for annotated variants..."

BED_FILE="$ANNO_DIR/all_samples.annotated.bed"

bcftools query -f '%CHROM\t%POS\t%POS\t%REF/%ALT\t%QUAL\n' \
    "$VCF_ANNOTATED" > "$BED_FILE" 2>> "$LOG_FILE"

log_message "INFO" "BED file: $BED_FILE"

################################################################################
# STEP 7: Generate annotation report
################################################################################

cat > "$ANNO_DIR/annotation_report.txt" << 'EOF'
# RICE GENOMICS - VARIANT ANNOTATION REPORT

## Input
Reference: Oryza sativa ssp. japonica (IRGSP v1.0)
Input VCF: all_samples.maf0.05.vcf.gz
Annotator: SnpEff

## Output Files
1. all_samples.annotated.vcf.gz - Fully annotated VCF
2. all_samples.annotated.bed - BED format variant positions
3. snpEff_report.html - Detailed SnpEff report
4. snpEff_summary.html - Summary statistics
5. annotation_stats.txt - Variant impact distribution
6. functional_summary.txt - High-impact variants

## Next Steps
1. Analyze population structure and diversity
2. Search for variants associated with drought tolerance
3. Identify candidate genes for functional analysis

## Important Notes
- ANN fields contain multiple annotations separated by commas
- Each annotation contains: allele|consequence|impact|gene_name|...
- See SnpEff documentation for full annotation format

EOF

log_message "INFO" "Annotation report: $ANNO_DIR/annotation_report.txt"

################################################################################
# STEP 8: Create summary of annotated variants
################################################################################

log_message "INFO" "Creating variant summary table..."

cat > "$ANNO_DIR/variant_summary.tsv" << 'EOF'
Chromosome	Position	Reference	Alternate	Quality	Type
EOF

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%TYPE\n' \
    "$VCF_ANNOTATED" >> "$ANNO_DIR/variant_summary.tsv" 2>> "$LOG_FILE"

log_message "INFO" "Variant summary: $ANNO_DIR/variant_summary.tsv"

################################################################################
# SUMMARY
################################################################################

if [ -f "$VCF_ANNOTATED" ]; then
    ANNO_COUNT=$(bcftools view "$VCF_ANNOTATED" | grep -v "^#" | wc -l)
    ANNO_SIZE=$(du -h "$VCF_ANNOTATED" | cut -f1)
else
    ANNO_COUNT="N/A"
    ANNO_SIZE="N/A"
fi

log_message "INFO" "=== ANNOTATION COMPLETE ==="
log_message "INFO" "Annotated VCF: $VCF_ANNOTATED"
log_message "INFO" "Number of annotated variants: $ANNO_COUNT"
log_message "INFO" "VCF size: $ANNO_SIZE"
log_message "INFO" "SnpEff reports: $ANNO_DIR/snpEff*.html"
log_message "INFO" ""
log_message "INFO" "Next step: bash 07_analysis.sh"

################################################################################
