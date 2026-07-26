#!/bin/bash

################################################################################
# RICE GENOMICS: VARIANT CALLING
# 05_variant_calling.sh
# Call SNPs and indels using bcftools mpileup/call
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

################################################################################
# STEP 1: Verify BAM files exist
################################################################################

log_message "INFO" "=== Starting Variant Calling ==="
log_message "INFO" "Checking BAM files..."

eval "$CONDA_ACTIVATE"

if [ ! -f "$BAM_DIR/bam_list.txt" ] || [ ! -s "$BAM_DIR/bam_list.txt" ]; then
    error_exit "BAM list not found. Run 04_mapping.sh first"
fi

BAM_COUNT=$(wc -l < "$BAM_DIR/bam_list.txt")
log_message "INFO" "Found $BAM_COUNT BAM files for variant calling"

################################################################################
# STEP 2: Call variants with bcftools
################################################################################

log_message "INFO" "Calling variants with bcftools mpileup/call..."
log_message "INFO" "Quality threshold: $VCF_QUAL_THRESHOLD"
log_message "INFO" "Depth thresholds: ${VCF_DP_MIN}-${VCF_DP_MAX}"

mkdir -p "$VCF_DIR"

# Create pileup from all samples
log_message "INFO" "Step 1: Computing pileup from all samples..."

PILEUP="$VCF_DIR/all_samples.pileup.bcf"

bcftools mpileup \
    -f "$REF_FASTA" \
    -b "$BAM_DIR/bam_list.txt" \
    -q 30 \
    -Q 20 \
    -d 500 \
    --threads "$BCFTOOLS_THREADS" \
    -o "$PILEUP" \
    -O b \
    2>> "$LOG_FILE"

log_message "INFO" "✓ Pileup created: $PILEUP"

################################################################################
# STEP 3: Call variants
################################################################################

log_message "INFO" "Step 2: Calling variants..."

VCF_RAW="$VCF_DIR/all_samples.raw.vcf.gz"

bcftools call \
    -m \
    -v \
    -o "$VCF_RAW" \
    -O z \
    --threads "$BCFTOOLS_THREADS" \
    "$PILEUP" \
    2>> "$LOG_FILE"

bcftools index "$VCF_RAW"

log_message "INFO" "✓ Raw VCF created: $VCF_RAW"

################################################################################
# STEP 4: Filter variants
################################################################################

log_message "INFO" "Step 3: Filtering variants..."

VCF_FILTERED="$VCF_DIR/all_samples.filtered.vcf.gz"

# Apply filters:
# - QUAL >= 30 (high confidence)
# - DP between 5 and 100 (depth filter)
# - MAF >= 0.05 (minor allele frequency, remove rare variants)

bcftools filter \
    -i "QUAL>=$VCF_QUAL_THRESHOLD & INFO/DP>=$VCF_DP_MIN & INFO/DP<=$VCF_DP_MAX" \
    -o "$VCF_FILTERED" \
    -O z \
    "$VCF_RAW" \
    2>> "$LOG_FILE"

bcftools index "$VCF_FILTERED"

log_message "INFO" "✓ Filtered VCF created: $VCF_FILTERED"

################################################################################
# STEP 5: Apply MAF filter
################################################################################

log_message "INFO" "Step 4: Applying MAF filter (threshold: $MAF_THRESHOLD)..."

VCF_MAF="$VCF_DIR/all_samples.maf${MAF_THRESHOLD}.vcf.gz"

# Count alleles and filter
bcftools view \
    --threads "$BCFTOOLS_THREADS" \
    -q "$MAF_THRESHOLD":minor \
    -o "$VCF_MAF" \
    -O z \
    "$VCF_FILTERED" \
    2>> "$LOG_FILE"

bcftools index "$VCF_MAF"

log_message "INFO" "✓ MAF-filtered VCF created: $VCF_MAF"

################################################################################
# STEP 6: Generate VCF statistics
################################################################################

log_message "INFO" "Computing VCF statistics..."

cat > "$VCF_DIR/variant_stats.txt" << 'EOF'
# VARIANT CALLING STATISTICS

## Raw VCF
EOF

bcftools view "$VCF_RAW" | grep -v "^#" | wc -l >> "$VCF_DIR/variant_stats.txt"

cat >> "$VCF_DIR/variant_stats.txt" << 'EOF'
 variants in raw VCF

## After QUAL/DP filtering
EOF

bcftools view "$VCF_FILTERED" | grep -v "^#" | wc -l >> "$VCF_DIR/variant_stats.txt"

cat >> "$VCF_DIR/variant_stats.txt" << 'EOF'
 variants

## After MAF filtering
EOF

bcftools view "$VCF_MAF" | grep -v "^#" | wc -l >> "$VCF_DIR/variant_stats.txt"

cat >> "$VCF_DIR/variant_stats.txt" << 'EOF'
 variants

## Variant types (final VCF)
EOF

# Count SNPs vs indels
SNP_COUNT=$(bcftools view "$VCF_MAF" | grep -v "^#" | awk '{print length($4) + length($5)}' | awk '$1==2 {count++} END {print count}')
INDEL_COUNT=$(bcftools view "$VCF_MAF" | grep -v "^#" | awk '{print length($4) + length($5)}' | awk '$1>2 {count++} END {print count}')

echo "SNPs: $SNP_COUNT" >> "$VCF_DIR/variant_stats.txt"
echo "Indels: $INDEL_COUNT" >> "$VCF_DIR/variant_stats.txt"

log_message "INFO" "Variant stats: $VCF_DIR/variant_stats.txt"

################################################################################
# STEP 7: Create per-sample VCFs (optional but useful)
################################################################################

log_message "INFO" "Extracting per-sample genotypes..."

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%GT\t]\n' \
    "$VCF_MAF" > "$VCF_DIR/genotypes.tsv" \
    2>> "$LOG_FILE"

log_message "INFO" "Genotype matrix: $VCF_DIR/genotypes.tsv"

################################################################################
# STEP 8: Create VCF summary report
################################################################################

cat > "$VCF_DIR/vcf_summary.txt" << 'EOF'
# VARIANT CALLING SUMMARY

Reference Genome: Oryza sativa ssp. japonica (IRGSP v1.0)
Variant Caller: bcftools mpileup/call

## Calling Parameters
- Minimum base quality: 20
- Minimum mapping quality: 30
- Maximum depth: 500
- VCF quality threshold: >= 30
- Depth filter: 5-100
- MAF threshold: >= 0.05

## Output Files
- raw VCF: all_samples.raw.vcf.gz
- filtered VCF: all_samples.filtered.vcf.gz
- final VCF: all_samples.maf0.05.vcf.gz
- genotypes matrix: genotypes.tsv
- statistics: variant_stats.txt

## Next Steps
1. Annotate variants (SnpEff)
2. Perform population genetics analyses
3. Generate visualizations
EOF

log_message "INFO" "VCF summary: $VCF_DIR/vcf_summary.txt"

################################################################################
# SUMMARY
################################################################################

FINAL_VARIANTS=$(bcftools view "$VCF_MAF" | grep -v "^#" | wc -l)
VCF_SIZE=$(du -h "$VCF_MAF" | cut -f1)

log_message "INFO" "=== VARIANT CALLING COMPLETE ==="
log_message "INFO" "Final VCF: $VCF_MAF"
log_message "INFO" "Number of variants: $FINAL_VARIANTS"
log_message "INFO" "VCF size: $VCF_SIZE"
log_message "INFO" ""
log_message "INFO" "Next step: bash 06_annotation.sh"

################################################################################
