#!/bin/bash

# ====================================================================
# QUANTIFICATION: Salmon (RNA-seq bulk / Smart-seq2)
# PRJEB55754 - Plasmodium falciparum Gamétocytes
# ====================================================================
# CORRECTION: flag --meta supprimé (métagénomique uniquement)
# ====================================================================

PROJECT_ROOT=~/Documents/github_projet/Plasmodium
input_dir="$PROJECT_ROOT/data_trimmed"
reference_dir="$PROJECT_ROOT/reference"
salmon_index_dir="$reference_dir/Pfalciparum3D7_salmon_index"
quantification_dir="$PROJECT_ROOT/quantification/salmon_output"
log_dir="$PROJECT_ROOT/logs"

sample_list=(
    ERR11471971
    ERR11471975
    ERR11471979
    ERR11471985
)

SALMON_THREADS=12
SALMON_LIBTYPE="A"

mkdir -p "$quantification_dir" "$log_dir"

echo "======================================================"
echo "QUANTIFICATION: Salmon"
echo "Plasmodium falciparum gamétocytes"
echo "======================================================"
echo ""

# ====================================================================
# 1. VÉRIFICATIONS
# ====================================================================

echo "======================================================"
echo "ÉTAPE 1: Vérifications préalables"
echo "======================================================"
echo ""

if ! command -v salmon &> /dev/null; then
    echo "❌ Erreur: Salmon non trouvé"
    exit 1
fi
echo "✓ Salmon: $(salmon --version)"

if [ ! -d "$salmon_index_dir" ]; then
    echo "❌ Index Salmon non trouvé: $salmon_index_dir"
    echo "   Exécutez d'abord: bash 02_salmon_indexing.sh"
    exit 1
fi
echo "✓ Index Salmon: $salmon_index_dir"

echo ""
echo "Vérification des fichiers FASTQ trimmés:"
all_present=true
for sample in "${sample_list[@]}"; do
    r1="$input_dir/${sample}_1.trimmed.fastq.gz"
    r2="$input_dir/${sample}_2.trimmed.fastq.gz"
    if [ ! -f "$r1" ] || [ ! -f "$r2" ]; then
        echo "  ❌ $sample: fichiers manquants dans $input_dir"
        all_present=false
    else
        r1_size=$(du -h "$r1" | cut -f1)
        r2_size=$(du -h "$r2" | cut -f1)
        echo "  ✓ $sample: R1=$r1_size, R2=$r2_size"
    fi
done

if [ "$all_present" = false ]; then
    echo "❌ Certains fichiers manquent. Lancez d'abord: bash 02b_fastp_trimming.sh"
    exit 1
fi

echo ""
echo "Paramètres Salmon:"
echo "  Library type: $SALMON_LIBTYPE (auto-detect)"
echo "  Threads: $SALMON_THREADS"
echo "  --validateMappings: activé"
echo "  NOTE: --meta supprimé (flag métagénomique, invalide ici)"
echo ""

# ====================================================================
# 2. QUANTIFICATION PAR ÉCHANTILLON
# ====================================================================

echo "======================================================"
echo "ÉTAPE 2: Quantification"
echo "======================================================"
echo ""

for sample in "${sample_list[@]}"; do
    R1="$input_dir/${sample}_1.trimmed.fastq.gz"
    R2="$input_dir/${sample}_2.trimmed.fastq.gz"
    sample_output="$quantification_dir/$sample"

    if [ -f "$sample_output/quant.sf" ]; then
        echo "⚠ $sample: déjà quantifié (ignoré)"
        continue
    fi

    echo "=========================================="
    echo "Traitement: $sample"
    echo "=========================================="
    echo "  R1: $(basename $R1)"
    echo "  R2: $(basename $R2)"
    echo "  Sortie: $sample_output"
    echo ""

    mkdir -p "$sample_output"

    # NOTE: --meta retiré. Options conservées:
    #   --validateMappings : améliore la précision du mapping
    #   -l A               : auto-détection de la strandedness
    salmon quant -i "$salmon_index_dir" \
                 -l "$SALMON_LIBTYPE" \
                 -1 "$R1" \
                 -2 "$R2" \
                 --validateMappings \
                 -p "$SALMON_THREADS" \
                 -o "$sample_output" \
                 2>&1 | tee "$log_dir/salmon_${sample}.log"

    if [ $? -eq 0 ] && [ -f "$sample_output/quant.sf" ]; then
        echo ""
        echo "✓ Quantification réussie: $sample"

        # Afficher le taux de mapping depuis le log JSON
        meta_json="$sample_output/aux_info/meta_info.json"
        if [ -f "$meta_json" ]; then
            mapping_rate=$(grep "percent_mapped" "$meta_json" | grep -oP '[0-9]+\.[0-9]+')
            total_reads=$(grep "num_processed" "$meta_json" | grep -oP '[0-9]+')
            echo "  Reads traités  : $total_reads"
            echo "  Taux de mapping: ${mapping_rate}%"
            
            # Avertir si taux trop bas
            if [ ! -z "$mapping_rate" ]; then
                threshold=50
                if (( $(echo "$mapping_rate < $threshold" | bc -l) )); then
                    echo "  ⚠ Taux de mapping faible (<${threshold}%)"
                    echo "     Vérifiez: bon organisme de référence? Reads contaminés?"
                fi
            fi
        fi
    else
        echo "❌ Erreur lors de la quantification de $sample"
        echo "   Consultez: $log_dir/salmon_${sample}.log"
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

echo "Fichiers quant.sf générés:"
for sample in "${sample_list[@]}"; do
    quant_file="$quantification_dir/$sample/quant.sf"
    if [ -f "$quant_file" ]; then
        n_lines=$(( $(wc -l < "$quant_file") - 1 ))
        echo "  ✓ $sample: $n_lines transcripts"
    else
        echo "  ❌ $sample: quant.sf absent"
    fi
done

echo ""
echo "======================================================"
echo "✓ QUANTIFICATION SALMON COMPLÉTÉE"
echo "======================================================"
echo ""

echo "[$(date)] Salmon quantification completed" >> "$log_dir/analysis.log"

echo "PROCHAINES ÉTAPES:"
echo "  python3 04_salmon_to_genecounts.py"
echo ""
