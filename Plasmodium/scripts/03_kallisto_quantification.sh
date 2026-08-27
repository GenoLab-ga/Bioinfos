#!/bin/bash

# ====================================================================
# QUANTIFICATION: Kallisto (mode bulk / Smart-seq2)
# PRJEB55754 - Plasmodium falciparum Gamétocytes
# ====================================================================
# CORRECTION: input_dir pointe sur data_trimmed/ (cohérence avec Salmon)
# ====================================================================

PROJECT_ROOT=~/Documents/github_projet/Plasmodium
input_dir="$PROJECT_ROOT/data_trimmed"        # CORRIGÉ: était data/ (reads bruts)
reference_dir="$PROJECT_ROOT/reference"
kallisto_index="$reference_dir/Pfalciparum3D7_kallisto_index.idx"
kallisto_output="$PROJECT_ROOT/quantification/kallisto_output"
log_dir="$PROJECT_ROOT/logs"

sample_list=(
    ERR11471971
    ERR11471975
    ERR11471979
    ERR11471985
)

THREADS=12

mkdir -p "$kallisto_output" "$log_dir"

echo "======================================================"
echo "QUANTIFICATION: Kallisto"
echo "======================================================"
echo ""

# ====================================================================
# 1. VÉRIFICATIONS
# ====================================================================

echo "======================================================"
echo "ÉTAPE 1: Vérifications"
echo "======================================================"
echo ""

if ! command -v kallisto &> /dev/null; then
    echo "❌ Erreur: kallisto non trouvé"
    echo "   Installation: mamba install kallisto"
    exit 1
fi
echo "✓ Kallisto: $(kallisto version)"

if [ ! -f "$kallisto_index" ]; then
    echo "❌ Index Kallisto non trouvé: $kallisto_index"
    echo "   Exécutez d'abord: bash 02_kallisto_indexing.sh"
    exit 1
fi
echo "✓ Index: $(basename $kallisto_index)"

echo ""
echo "Vérification des fichiers FASTQ trimmés:"
for sample in "${sample_list[@]}"; do
    r1="$input_dir/${sample}_1.trimmed.fastq.gz"
    r2="$input_dir/${sample}_2.trimmed.fastq.gz"
    if [ -f "$r1" ] && [ -f "$r2" ]; then
        echo "  ✓ $sample"
    else
        echo "  ❌ $sample: manquant dans $input_dir"
        echo "     Lancez d'abord: bash 02b_fastp_trimming.sh"
        exit 1
    fi
done

echo ""

# ====================================================================
# 2. QUANTIFICATION PAR ÉCHANTILLON
# ====================================================================

echo "======================================================"
echo "ÉTAPE 2: Quantification Kallisto (mode paired-end)"
echo "======================================================"
echo ""

for sample in "${sample_list[@]}"; do
    R1="$input_dir/${sample}_1.trimmed.fastq.gz"
    R2="$input_dir/${sample}_2.trimmed.fastq.gz"
    sample_output="$kallisto_output/$sample"

    if [ -f "$sample_output/abundance.tsv" ]; then
        echo "⚠ $sample: déjà quantifié (ignoré)"
        continue
    fi

    mkdir -p "$sample_output"

    echo "=========================================="
    echo "Traitement: $sample"
    echo "=========================================="

    kallisto quant \
        -i "$kallisto_index" \
        -o "$sample_output" \
        -t "$THREADS" \
        "$R1" "$R2" \
        2>&1 | tee "$log_dir/kallisto_${sample}.log"

    if [ -f "$sample_output/abundance.tsv" ]; then
        echo "✓ Quantification réussie"

        # Extraire stats depuis le log
        n_processed=$(grep "processed" "$log_dir/kallisto_${sample}.log" | grep -oP '[0-9,]+ reads' | head -1)
        n_pseudoaligned=$(grep "pseudoaligned" "$log_dir/kallisto_${sample}.log" | grep -oP '[0-9,]+ reads' | head -1)
        echo "  Reads traités      : $n_processed"
        echo "  Reads pseudoalignés: $n_pseudoaligned"
    else
        echo "❌ Erreur quantification $sample"
        echo "   Consultez: $log_dir/kallisto_${sample}.log"
        exit 1
    fi

    echo ""
done

echo ""

# ====================================================================
# 3. RÉSUMÉ
# ====================================================================

echo "======================================================"
echo "RÉSUMÉ"
echo "======================================================"
echo ""

echo "Fichiers abundance.tsv générés:"
for sample in "${sample_list[@]}"; do
    f="$kallisto_output/$sample/abundance.tsv"
    if [ -f "$f" ]; then
        n=$(( $(wc -l < "$f") - 1 ))
        echo "  ✓ $sample: $n transcripts"
    else
        echo "  ❌ $sample manquant"
    fi
done

echo ""
echo "======================================================"
echo "✓ QUANTIFICATION KALLISTO COMPLÉTÉE"
echo "======================================================"
echo ""

echo "[$(date)] Kallisto quantification completed" >> "$log_dir/analysis.log"

echo "PROCHAINES ÉTAPES:"
echo "  bash 04_kallisto_to_genecounts.sh"
echo ""
