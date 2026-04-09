#!/bin/bash

# Liste des identifiants
srr_list=(
SRR13687743 SRR13687744 SRR13687745 SRR13687746 SRR13687747 SRR13687748 SRR13687749 SRR13687750 SRR13687751 # Témoins sans rifampicine
SRR13687752 SRR13687753 SRR13687754 SRR13687755 SRR13687756 SRR13687757 SRR13687758 SRR13687759 SRR13687760 # Traités à 0.008µg/ml
SRR13687761 SRR13687762 SRR13687763 SRR13687764 SRR13687765 SRR13687766 SRR13687767 SRR13687768 SRR13687769 # Traités à 0.02µg/ml
)

# Repertoires de travail
mkdir -p qc_reports nettoyage mapping reference comptage qc_clean

input_dir="raw_data"
qc_dir="qc_reports"
clean_dir="nettoyage"
aligned_dir="mapping"
reference_dir="reference"
genome_ref="$reference_dir/h37rv.fna"
gtf_ref="$reference_dir/annotation.gtf"
compte_dir="comptage"

# Traitement 
# ================================================================
# 1) Contrôle qualité des données brutes
# ================================================================

echo "=== Etape 1: Contrôle qualité avec FastQC ==="
for srr_id in "${srr_list[@]}"; do
    fastqc "$input_dir/${srr_id}_1.fastq.gz" "$input_dir/${srr_id}_2.fastq.gz" -o "$qc_dir" -t 4
done
multiqc "$qc_dir" -o "$qc_dir"


# ================================================================
# 2) Nettoyage de reads + suppression des duplications
# ================================================================
echo "=== Etape 2: Nettoyage avec Fastp ==="
for srr_id in "${srr_list[@]}"; do
   fastp -i "$input_dir/${srr_id}_1.fastq.gz" \
         -I "$input_dir/${srr_id}_2.fastq.gz" \
         -o "$clean_dir/${srr_id}_1_clean.fastq.gz" \
         -O "$clean_dir/${srr_id}_2_clean.fastq.gz" \
         --dedup \
         --dup_calc_accuracy 4 \
         --thread 4 \
         --html "$clean_dir/${srr_id}_fastp.html" \
         --json "$clean_dir/${srr_id}_fastp.json"
done
   
   # 2.2) Contrôle qualité post-nettoyage
# Cette étape de contrôle qualité n'est pas forcement necessaire, vu que les reads sont déjà de bonne qualité. Mais si vous voulez confirmer, alors faites-le. Vous allez constatez aucun changement.

echo "=== Etape 2.2 : Contrôle qualité post-nettoyage =="

for srr_id in "${srr_list[@]}"; do
     fastqc "$clean_dir/${srr_id}_1_clean.fastq.gz" "$clean_dir/${srr_id}_2_clean.fastq.gz" -o qc_clean -t 4
done
multiqc qc_clean -o qc_clean
   
# ================================================================
# 3) Indexation du génome de référence
# ================================================================
echo "=== Etape 3 : Indexation du génome de référence ==="
if [ ! -f "$reference_dir/h37rv_index.1.bt2" ]; then
     bowtie2-build "$genome_ref" "$reference_dir/h37rv_index"
fi

# ================================================================
# 4) Alignement
# ================================================================
echo "=== Etape 4 : Alignement avec Bowtie2 ==="
for srr_id in "${srr_list[@]}"; do
    bowtie2 -x "$reference_dir/h37rv_index" \
            -1 "$clean_dir/${srr_id}_1_clean.fastq.gz" \
            -2 "$clean_dir/${srr_id}_2_clean.fastq.gz" \
            -S "$aligned_dir/${srr_id}.sam" \
            --very-sensitive -p 8
            
    # Conversion SAM > BAM trié       
    samtools view -bS "$aligned_dir/${srr_id}.sam" | samtools sort -o "$aligned_dir/${srr_id}.sorted.bam"
    
    # Indexation 
    samtools index "$aligned_dir/${srr_id}.sorted.bam"
    
    # Suppression du SAM
    rm -f "$aligned_dir/${srr_id}.sam"
            
done

# ================================================================
# 5) Comptage des reads
# ================================================================
echo "=== Étape 4 : Comptage avec featureCounts ==="
featureCounts -T 8 -p -t gene -g gene_id \
    -a reference/annotation.gtf \
    -o comptage/counts.txt \
    mapping/*.sorted.bam

echo ""
echo "=== Résumé ==="
cat comptage/counts.txt.summary

echo ""
echo "=== Aperçu des comptages ==="
head -20 comptage/counts.txt

echo "=== Alignement et comptage terminés ! ==="
