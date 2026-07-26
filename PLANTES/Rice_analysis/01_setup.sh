#!/bin/bash

################################################################################
# RICE GENOMICS: SETUP & INITIALIZATION
# 01_setup.sh
# Initialise l'environnement, vérifie les outils, prépare les répertoires
################################################################################

set -euo pipefail

# Load configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

################################################################################
# STEP 1: Initialize directories and logging
################################################################################

log_message "INFO" "Starting Rice Genomics Pipeline Setup"
print_config

init_project

################################################################################
# STEP 2: Check disk space (require 500GB minimum)
################################################################################

log_message "INFO" "Checking disk space requirements..."
check_disk_space 500

################################################################################
# STEP 3: Verify mamba environment exists
################################################################################

log_message "INFO" "Verifying mamba environment: $MAMBA_ENV"

if eval "$CONDA_ACTIVATE" &> /dev/null; then
    log_message "INFO" "Mamba environment activated successfully"
else
    error_exit "Failed to activate mamba environment: $MAMBA_ENV"
fi

################################################################################
# STEP 4: Check essential tools
################################################################################

log_message "INFO" "Checking bioinformatics tools..."

eval "$CONDA_ACTIVATE"

tools=(
    "bwa"
    "samtools"
    "bcftools"
    "bedtools"
    "fastqc"
    "trimmomatic"
    "vcftools"
)

missing_tools=()

for tool in "${tools[@]}"; do
    if check_tool "$tool"; then
        log_message "INFO" "✓ Found: $tool ($(command -v $tool))"
    else
        log_message "WARN" "✗ Missing: $tool"
        missing_tools+=("$tool")
    fi
done

if [ ${#missing_tools[@]} -gt 0 ]; then
    log_message "WARN" "Missing tools: ${missing_tools[*]}"
    log_message "WARN" "Install them with: mamba install -n genomics_env ${missing_tools[*]}"
fi

################################################################################
# STEP 5: Download and index reference genome
################################################################################

log_message "INFO" "Checking reference genome..."

if [ ! -f "$REF_FASTA" ]; then
    log_message "INFO" "Reference genome not found. Downloading IRGSP v1.0..."
    
    cd "$REF_DIR"
    
    # Download genome
    log_message "INFO" "Downloading reference FASTA..."
    wget -q --show-progress "$REF_GENOME_URL" -O IRGSP-1.0_genome.fasta.gz
    
    # Decompress
    log_message "INFO" "Decompressing..."
    gunzip -f IRGSP-1.0_genome.fasta.gz
    
    # Download annotation
    log_message "INFO" "Downloading GFF3 annotation..."
    wget -q --show-progress "$REF_ANNOTATION_URL" -O IRGSP-1.0.gff3.gz
    gunzip -f IRGSP-1.0.gff3.gz
    
    log_message "INFO" "✓ Reference genome downloaded"
else
    log_message "INFO" "✓ Reference genome exists: $REF_FASTA"
fi

################################################################################
# STEP 6: Index reference genome for BWA
################################################################################

if [ ! -f "$REF_INDEX" ]; then
    log_message "INFO" "Indexing reference genome for BWA..."
    bwa index -p "$REF_DIR/IRGSP-1.0_genome" "$REF_FASTA" 2>> "$LOG_FILE"
    log_message "INFO" "✓ BWA index created"
else
    log_message "INFO" "✓ BWA index already exists"
fi

################################################################################
# STEP 7: Create sequence dictionary for GATK
################################################################################

if [ ! -f "$REF_DICT" ]; then
    log_message "INFO" "Creating sequence dictionary..."
    samtools dict "$REF_FASTA" > "$REF_DICT"
    log_message "INFO" "✓ Sequence dictionary created"
else
    log_message "INFO" "✓ Sequence dictionary exists"
fi

################################################################################
# STEP 8: Index reference with samtools
################################################################################

if [ ! -f "$REF_FASTA.fai" ]; then
    log_message "INFO" "Creating samtools index..."
    samtools faidx "$REF_FASTA"
    log_message "INFO" "✓ Samtools index created"
else
    log_message "INFO" "✓ Samtools index exists"
fi

################################################################################
# STEP 9: Create samples manifest file
################################################################################

if [ ! -f "$SAMPLES_FILE" ]; then
    log_message "INFO" "Creating samples manifest template..."
    
    cat > "$SAMPLES_FILE" << 'EOF'
# Samples manifest for Rice Genomics Analysis
# Format: SRR_ID,SAMPLE_NAME,SUBSPECIES,COUNTRY
# Download FASTQ files using: prefetch and fasterq-dump from SRA Toolkit
# Example:
# SRR1234567,IR64,indica,India
# SRR2345678,Nipponbare,japonica,Japan

SRR8817819,IRIS_313-10118,indica,India
SRR8817821,IRIS_313-10177,indica,India
SRR8817823,IRIS_313-10272,indica,India
SRR8817825,IRIS_313-10312,indica,India
SRR8817827,IRIS_313-10390,indica,India
SRR8817829,IRIS_313-10406,japonica,Japan
SRR8817831,IRIS_313-10450,japonica,Japan
SRR8817833,IRIS_313-10474,japonica,Japan
SRR8817835,IRIS_313-10529,japonica,Japan
SRR8817837,IRIS_313-10576,japonica,Japan
EOF
    
    log_message "INFO" "✓ Sample manifest created: $SAMPLES_FILE"
    log_message "INFO" "Edit this file to add your 15 accessions"
else
    log_message "INFO" "✓ Sample manifest exists"
fi

################################################################################
# STEP 10: Create README
################################################################################

cat > "$PROJECT_DIR/README.md" << 'EOF'
# Rice Genomics WGS Analysis Pipeline

## Overview
Complete pipeline for analyzing Whole Genome Sequencing (WGS) data from the 3000 Rice Genomes Project (3K-RGP).

### Pipeline Steps:
1. **01_setup.sh** - Initialize project, download reference genome
2. **02_download_data.sh** - Download FASTQ files from SRA
3. **03_qc_preprocessing.sh** - Quality control and trimming
4. **04_mapping.sh** - Map reads to reference with BWA
5. **05_variant_calling.sh** - Call variants with bcftools
6. **06_annotation.sh** - Annotate variants with SnpEff
7. **07_analysis.sh** - Population genetics analysis
8. **08_visualization.R** - Generate plots

## Getting Started

### 1. Edit samples.txt
Add your 15 rice accessions:
```bash
nano samples.txt
```

### 2. Run setup
```bash
bash 01_setup.sh
```

### 3. Download data
```bash
bash 02_download_data.sh
```

### 4. Execute pipeline
```bash
bash 03_qc_preprocessing.sh
bash 04_mapping.sh
bash 05_variant_calling.sh
bash 06_annotation.sh
bash 07_analysis.sh
```

## Requirements
- Mamba/Conda with genomics_env activated
- 16GB RAM minimum
- 500GB disk space
- 12 CPUs recommended

## Tools Required
- bwa (alignment)
- samtools (BAM processing)
- bcftools (variant calling)
- fastqc (quality control)
- trimmomatic (read trimming)
- snpEff (annotation)
- bedtools, vcftools (utilities)

## Output Structure
```
results/
├── 01_qc/              # FastQC reports
├── 02_trimmed/         # Trimmed FASTQ
├── 03_bam/             # Aligned BAM files
├── 04_vcf/             # VCF variants
├── 05_annotation/      # SnpEff annotated
└── 06_analysis/        # Population genetics results
```

## Contact
Karl Mounguele | GenoLabGab | karl@genolabgab.com
EOF

log_message "INFO" "✓ README.md created"

################################################################################
# SUMMARY
################################################################################

log_message "INFO" "=== SETUP COMPLETE ==="
log_message "INFO" "Next steps:"
log_message "INFO" "1. Edit samples.txt to add your 15 rice accessions"
log_message "INFO" "2. Run: bash 02_download_data.sh"
log_message "INFO" "Pipeline log: $LOG_FILE"

################################################################################
