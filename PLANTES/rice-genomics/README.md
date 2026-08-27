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
