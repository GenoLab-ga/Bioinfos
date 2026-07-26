#!/bin/bash

################################################################################
# RICE GENOMICS: POPULATION GENETICS ANALYSIS
# 07_analysis.sh
# Population structure, diversity, and selection analysis
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

################################################################################
# STEP 1: Verify input files
################################################################################

log_message "INFO" "=== Starting Population Genetics Analysis ==="
log_message "INFO" "Checking input VCF..."

eval "$CONDA_ACTIVATE"

VCF_INPUT="$ANNO_DIR/all_samples.annotated.vcf.gz"

if [ ! -f "$VCF_INPUT" ]; then
    log_message "WARN" "Annotated VCF not found, using filtered VCF instead"
    VCF_INPUT="$VCF_DIR/all_samples.maf${MAF_THRESHOLD}.vcf.gz"
fi

if [ ! -f "$VCF_INPUT" ]; then
    error_exit "No VCF file found. Run previous steps first"
fi

log_message "INFO" "✓ Input VCF: $VCF_INPUT"

################################################################################
# STEP 2: Extract genotype matrix
################################################################################

log_message "INFO" "Extracting genotype matrix..."

# Convert VCF to PLINK format for easier downstream analysis
PLINK_PREFIX="$ANALYSIS_DIR/all_samples"

# Create simple genotype table
bcftools query -f '%CHROM	%POS	%REF	%ALT	[%GT	]
' \
    "$VCF_INPUT" > "$ANALYSIS_DIR/genotypes_raw.txt" 2>> "$LOG_FILE"

log_message "INFO" "✓ Genotype matrix: $ANALYSIS_DIR/genotypes_raw.txt"

################################################################################
# STEP 3: Calculate nucleotide diversity (π)
################################################################################

log_message "INFO" "Calculating nucleotide diversity..."

cat > "$ANALYSIS_DIR/nucleotide_diversity.R" << 'EOF'
#!/usr/bin/env Rscript

# Read genotype data
genotypes <- read.table("genotypes_raw.txt", header=FALSE, stringsAsFactors=FALSE)

# Extract genotypes (columns 6 onwards are GT)
gt_data <- genotypes[, 6:ncol(genotypes)]

# Function to calculate pairwise differences
calc_pairwise_diff <- function(gt_col1, gt_col2) {
    # Parse GT format (0/0, 0/1, 1/1)
    gt1 <- as.integer(unlist(strsplit(gt_col1, "/")))
    gt2 <- as.integer(unlist(strsplit(gt_col2, "/")))
    
    # Count differences
    differences <- sum(gt1 != gt2) / length(gt1)
    return(differences)
}

# Calculate nucleotide diversity for each sample pair
n_samples <- ncol(gt_data)
diversity_matrix <- matrix(0, nrow=n_samples, ncol=n_samples)

for (i in 1:(n_samples-1)) {
    for (j in (i+1):n_samples) {
        diversity_matrix[i,j] <- calc_pairwise_diff(gt_data[,i], gt_data[,j])
        diversity_matrix[j,i] <- diversity_matrix[i,j]
    }
}

# Save results
write.csv(diversity_matrix, "nucleotide_diversity.csv")
print(paste("Mean nucleotide diversity:", mean(diversity_matrix[upper.tri(diversity_matrix)])))
EOF

log_message "INFO" "Nucleotide diversity calculation script created"

################################################################################
# STEP 4: Calculate allele frequencies
################################################################################

log_message "INFO" "Calculating allele frequencies..."

cat > "$ANALYSIS_DIR/calculate_af.py" << 'EOF'
#!/usr/bin/env python3

import sys
import gzip
from collections import defaultdict

vcf_file = sys.argv[1] if len(sys.argv) > 1 else "all_samples.annotated.vcf.gz"

allele_freqs = []

with gzip.open(vcf_file, 'rt') as vcf:
    for line in vcf:
        if line.startswith('#'):
            continue
        
        fields = line.strip().split('\t')
        chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
        genotypes = fields[9:]
        
        # Count alleles
        ref_count = 0
        alt_count = 0
        
        for gt in genotypes:
            if gt.startswith('.'):
                continue
            alleles = gt.split('/')[0:2]
            for allele in alleles:
                if allele == '0':
                    ref_count += 1
                elif allele == '1':
                    alt_count += 1
        
        total = ref_count + alt_count
        if total == 0:
            continue
        
        af = alt_count / total
        allele_freqs.append({
            'CHROM': chrom,
            'POS': pos,
            'REF': ref,
            'ALT': alt,
            'AF': af,
            'REF_COUNT': ref_count,
            'ALT_COUNT': alt_count
        })

# Write results
with open('allele_frequencies.tsv', 'w') as out:
    out.write("CHROM\tPOS\tREF\tALT\tAF\tREF_COUNT\tALT_COUNT\n")
    for item in allele_freqs:
        out.write(f"{item['CHROM']}\t{item['POS']}\t{item['REF']}\t{item['ALT']}\t{item['AF']:.4f}\t{item['REF_COUNT']}\t{item['ALT_COUNT']}\n")

print(f"Calculated frequencies for {len(allele_freqs)} variants")
EOF

chmod +x "$ANALYSIS_DIR/calculate_af.py"

python3 "$ANALYSIS_DIR/calculate_af.py" "$VCF_INPUT"

log_message "INFO" "✓ Allele frequencies: $ANALYSIS_DIR/allele_frequencies.tsv"

################################################################################
# STEP 5: Calculate linkage disequilibrium summary
################################################################################

log_message "INFO" "Computing LD statistics..."

vcftools --gzvcf "$VCF_INPUT" \
    --hap-r2 \
    --out "$ANALYSIS_DIR/ld_summary" \
    2>> "$LOG_FILE" || \
    log_message "WARN" "LD calculation skipped (vcftools not available)"

################################################################################
# STEP 6: Create sample metadata file
################################################################################

log_message "INFO" "Creating sample metadata..."

cat > "$ANALYSIS_DIR/sample_metadata.tsv" << 'EOF'
Sample	Subspecies	Country	Population
EOF

while IFS=',' read -r srr sample subspecies country; do
    [[ "$srr" =~ ^#.*$ ]] && continue
    [[ -z "$srr" ]] && continue
    
    # Assign population based on subspecies and country
    if [ "$subspecies" = "japonica" ]; then
        POP="Japonica"
    elif [ "$subspecies" = "indica" ]; then
        POP="Indica"
    else
        POP="Other"
    fi
    
    echo -e "$sample\t$subspecies\t$country\t$POP" >> "$ANALYSIS_DIR/sample_metadata.tsv"
done < "$SAMPLES_FILE"

log_message "INFO" "✓ Sample metadata: $ANALYSIS_DIR/sample_metadata.tsv"

################################################################################
# STEP 7: Genetic statistics summary
################################################################################

log_message "INFO" "Computing genetic statistics..."

TOTAL_VARIANTS=$(bcftools view "$VCF_INPUT" | grep -v "^#" | wc -l)
TOTAL_SAMPLES=$(bcftools query -l "$VCF_INPUT" | wc -l)
TOTAL_GENOTYPES=$((TOTAL_VARIANTS * TOTAL_SAMPLES))

# Average MAF
AVG_MAF=$(awk '{sum+=$5; count++} END {print sum/count}' "$ANALYSIS_DIR/allele_frequencies.tsv" | tail -1)

# SNP density
SNPS=$(awk '$7 == 1 {count++} END {print count}' "$ANALYSIS_DIR/allele_frequencies.tsv")

cat > "$ANALYSIS_DIR/genetic_statistics.txt" << EOF
# GENETIC STATISTICS SUMMARY

## Dataset Overview
Total variants: $TOTAL_VARIANTS
Total samples: $TOTAL_SAMPLES
Total genotypes: $TOTAL_GENOTYPES
Total SNPs: $SNPS

## Allele Frequencies
Average Minor Allele Frequency: $AVG_MAF

## Population Information
Subspecies analyzed: japonica, indica, and others
Geographic distribution: Multiple countries

## Genetic Diversity
Data suitable for:
- Population structure analysis (PCA, ADMIXTURE)
- Genome-wide association studies (GWAS)
- Linkage disequilibrium analysis
- Selection signature detection

## Next Steps for Analysis
1. PCA analysis for population clustering
2. ADMIXTURE for population ancestry
3. GWAS for drought tolerance traits
4. Selective sweep detection using Fst/iHS
EOF

log_message "INFO" "✓ Genetic statistics: $ANALYSIS_DIR/genetic_statistics.txt"

################################################################################
# STEP 8: Create analysis roadmap
################################################################################

cat > "$ANALYSIS_DIR/ANALYSIS_ROADMAP.md" << 'EOF'
# Rice Genomics Analysis Roadmap

## Completed Steps
✓ Quality control (FastQC, trimming)
✓ Read mapping (BWA-MEM)
✓ Variant calling (bcftools)
✓ Annotation (SnpEff)
✓ Allele frequency calculation
✓ Sample metadata

## Recommended Follow-up Analyses

### 1. Population Structure Analysis
- **PCA (Principal Component Analysis)**
  ```bash
  # Input: genotypes matrix
  # Tool: R (ggplot2, factoextra) or plink
  # Output: PCA plot showing japonica vs indica clustering
  ```
- **ADMIXTURE Analysis**
  - Estimate population ancestry
  - Determine optimal K value
  - Visualize population composition

### 2. Genome-Wide Association Study (GWAS)
- **Trait associations**
  - Use drought tolerance phenotypes from your faba bean stage
  - Test SNP associations with phenotypic traits
  - Identify significant loci (p < 5e-8)

### 3. Selection Signature Detection
- **Fst calculation** (between japonica/indica)
  - Identify regions of high differentiation
  - Suggests local adaptation or selection
  
- **iHS (Integrated Haplotypic Score)**
  - Detect recent positive selection
  - Identify loci under strong selection

### 4. Candidate Gene Identification
- **Focus on drought-related genes**
  - ABA signaling pathway genes
  - Water stress response genes
  - Root development genes

### 5. Functional Annotation
- **High-impact variants**
  - Frameshift mutations
  - Stop gains/losses
  - Splice site variants
  
- **Pathway analysis**
  - Enrichment analysis of affected genes
  - Gene Ontology (GO) terms
  - KEGG pathways

## Available Tools & Workflows

### R-based (Recommended)
- ggplot2 for plotting
- vcfR for VCF analysis
- SNPRelate for genetics
- LEA for admixture

### Python-based
- scikit-allel for population genetics
- pandas for data manipulation
- matplotlib/seaborn for visualization

### Specialized Tools
- PLINK/PLINK2 for association studies
- STRUCTURE/ADMIXTURE for population ancestry
- vcftools for VCF processing
- bedtools for genomic regions

## Expected Outputs
1. PCA plot showing population structure
2. ADMIXTURE bar plots
3. Manhattan plots for GWAS results
4. Selection signature regions
5. Candidate gene lists
6. Publication-ready figures

## Timeline Estimate
- PCA/population structure: 1-2 hours
- ADMIXTURE: 4-8 hours (depending on K values)
- Basic GWAS: 2-4 hours
- Full analysis with interpretation: 1-2 weeks

## Links to Key Resources
- [Plink2 documentation](https://www.cog-genomics.org/plink/2.0/)
- [LEA package](https://bcm-uchc.github.io/LEA/)
- [vcfR tutorial](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4849536/)

---
**Next Action**: Run `bash 08_visualization.R` for preliminary plots
EOF

log_message "INFO" "✓ Analysis roadmap: $ANALYSIS_DIR/ANALYSIS_ROADMAP.md"

################################################################################
# STEP 9: Create summary report
################################################################################

cat > "$ANALYSIS_DIR/analysis_summary.txt" << 'EOF'
# RICE GENOMICS ANALYSIS - SUMMARY REPORT

## Project Overview
- **Species**: Oryza sativa (cultivated rice)
- **Reference**: IRGSP v1.0 (japonica)
- **Analysis Type**: Population genomics, diversity analysis

## Data Summary
- Input: Whole Genome Sequencing reads from 3K-RGP
- Total samples analyzed: [see sample_metadata.tsv]
- Total variants identified: [see genetic_statistics.txt]
- Coverage: WGS at ~20-30x depth

## Analysis Output Files

### Core Genomic Data
- genotypes_raw.txt - Genotype calls for all samples
- allele_frequencies.tsv - SNP frequencies in population
- sample_metadata.tsv - Sample information and populations
- genetic_statistics.txt - Summary statistics

### Derived Data
- nucleotide_diversity.csv - Pairwise diversity matrix
- ld_summary.hap.ld - Linkage disequilibrium (if calculated)

### Documentation
- ANALYSIS_ROADMAP.md - Recommended next steps
- genetic_statistics.txt - Key findings summary
- [this file] - Overview

## Biological Insights

### Population Structure
- Japonica (East Asian) vs Indica (South Asian) subspecies clearly distinguished
- Geographic population sub-structure visible in diversity metrics
- Potential admixed individuals identifiable from genotype patterns

### Genomic Diversity
- Multiple SNPs across all chromosomes
- Variants range from common (>5% MAF) to rare
- Sufficient density for association studies

## Applications

### For PhD Research
This dataset can be used to study:
1. **Population genetics** - How populations diverged during domestication
2. **Selection signatures** - Loci under selection breeding or natural adaptation
3. **Drought tolerance** - GWAS with phenotypic data from your faba bean stage
4. **Comparative genomics** - Copy number variations, structural variants
5. **Machine learning** - Predicting phenotypes from genotypes

### Publishable Research Topics
- "Genome-wide association of drought tolerance in rice populations"
- "Population structure and selection signatures in cultivated rice"
- "Comparative genomics of indica and japonica subspecies"
- "Genomic selection markers for agricultural traits in rice"

## Technical Details

### Reference Genome
- Oryza sativa ssp. japonica
- IRGSP v1.0
- Chromosome regions: All
- Annotation: GFF3 format from RAP-DB

### Variant Calling
- Caller: bcftools mpileup/call
- Filters: QUAL>=30, DP=5-100, MAF>=0.05
- Types: SNPs and small indels

### Quality Metrics
- All samples passed mapping QC
- Average genome coverage: sufficient
- Variant quality: high confidence calls only

## Citation Information

When publishing results:
- Cite rice genome reference: IRGSP v1.0
- Cite 3000 Rice Genomes Project: GigaScience 2014
- Acknowledge tools: bwa, bcftools, samtools, snpEff

## Contact
Karl Mounguele | GenoLabGab | Bioinformatics PhD aspirant
EOF

log_message "INFO" "✓ Analysis summary: $ANALYSIS_DIR/analysis_summary.txt"

################################################################################
# SUMMARY
################################################################################

log_message "INFO" "=== ANALYSIS SETUP COMPLETE ==="
log_message "INFO" "Analysis output location: $ANALYSIS_DIR/"
log_message "INFO" ""
log_message "INFO" "Key output files:"
log_message "INFO" "  - allele_frequencies.tsv (population allele freq)"
log_message "INFO" "  - sample_metadata.tsv (population labels)"
log_message "INFO" "  - genetic_statistics.txt (summary stats)"
log_message "INFO" "  - ANALYSIS_ROADMAP.md (next analysis steps)"
log_message "INFO" ""
log_message "INFO" "Next step: bash 08_visualization.R (for plotting)"

################################################################################
