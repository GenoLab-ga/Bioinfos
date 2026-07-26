#!/bin/bash

# ====================================================================
# PIPELINE SINGLE-CELL RNA-SEQ: Plasmodium falciparum Gamétocytes
# PRJEB55754 - 4 échantillons
# ====================================================================

# ====================================================================
# 0. VARIABLES ET CHEMINS
# ====================================================================

# Chemins principaux
PROJECT_ROOT=~/Documents/github_projet/Plasmodium
input_dir="$PROJECT_ROOT/data"
fastqc_dir="$PROJECT_ROOT/qc_reports/fastqc"
multiqc_dir="$PROJECT_ROOT/qc_reports"
log_dir="$PROJECT_ROOT/logs"

# Échantillons à traiter
sample_list=(
    ERR11471971  # Day_1_pf_gam_tc_v2 (CONTRÔLE)
    ERR11471975  # Unknown
    ERR11471979  # Unknown
    ERR11471985  # Unknown
)

# Créer les répertoires s'ils n'existent pas
mkdir -p "$fastqc_dir" "$multiqc_dir" "$log_dir"

echo "======================================================"
echo "PIPELINE SINGLE-CELL: Plasmodium falciparum"
echo "Projet: PRJEB55754"
echo "======================================================"
echo ""
echo "Répertoire projet: $PROJECT_ROOT"
echo "Échantillons: ${#sample_list[@]}"
echo ""

# ====================================================================
# 1. CONTRÔLE QUALITÉ DES READS BRUTES (FastQC)
# ====================================================================

echo "======================================================"
echo "ÉTAPE 1: FastQC - Contrôle qualité des reads brutes"
echo "======================================================"
echo ""

# Vérifier que FastQC est disponible
if ! command -v fastqc &> /dev/null; then
    echo "❌ Erreur: FastQC non trouvé"
    echo "Solution: mamba install fastqc"
    exit 1
fi

echo "Version FastQC: $(fastqc --version)"
echo ""

# Déterminer le nombre de threads
THREADS=$(nproc)
echo "Threads disponibles: $THREADS"
echo ""

# Vérifier les fichiers d'entrée
echo "Vérification des fichiers FASTQ..."
for sample in "${sample_list[@]}"; do
    if [ ! -f "$input_dir/${sample}_1.fastq.gz" ] || [ ! -f "$input_dir/${sample}_2.fastq.gz" ]; then
        echo "❌ Erreur: Fichiers manquants pour $sample"
        echo "   Attendu: $input_dir/${sample}_1.fastq.gz et _2.fastq.gz"
        exit 1
    fi
done
echo "✓ Tous les fichiers présents"
echo ""

# Afficher les échantillons
echo "Échantillons à traiter:"
for sample in "${sample_list[@]}"; do
    size1=$(du -h "$input_dir/${sample}_1.fastq.gz" | cut -f1)
    size2=$(du -h "$input_dir/${sample}_2.fastq.gz" | cut -f1)
    echo "  $sample: R1=$size1, R2=$size2"
done
echo ""

# Boucle FastQC sur chaque échantillon
echo "Exécution FastQC (peut prendre 10-20 minutes)..."
echo ""

for sample in "${sample_list[@]}"; do
    echo "Traitement: $sample"
    
    fastqc "$input_dir/${sample}_1.fastq.gz" \
           "$input_dir/${sample}_2.fastq.gz" \
           --outdir "$fastqc_dir" \
           --threads "$THREADS" \
           --nogroup \
           -q
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Terminé"
    else
        echo "  ❌ Erreur sur $sample"
        exit 1
    fi
done

echo ""
echo "✓ FastQC complété pour tous les échantillons"
echo ""

# ====================================================================
# 2. AGRÉGATION DES RAPPORTS (MultiQC)
# ====================================================================

echo "======================================================"
echo "ÉTAPE 2: MultiQC - Agrégation des rapports FastQC"
echo "======================================================"
echo ""

# Vérifier MultiQC
if ! command -v multiqc &> /dev/null; then
    echo "⚠ MultiQC non trouvé. Installation..."
    pip install multiqc --quiet
fi

echo "Création du rapport MultiQC..."
cd "$fastqc_dir"

multiqc . --outdir "$multiqc_dir" \
          --force \
          --title "Single-cell RNA-seq QC: Pf gamétocytes (PRJEB55754)" \
          --quiet

if [ -f "$multiqc_dir/multiqc_report.html" ]; then
    echo "✓ Rapport MultiQC généré"
else
    echo "⚠ Erreur MultiQC"
fi

cd - > /dev/null

echo ""

# ====================================================================
# 3. STATISTIQUES RÉSUMÉES
# ====================================================================

echo "======================================================"
echo "ÉTAPE 3: Résumé statistiques"
echo "======================================================"
echo ""

echo "Taille totale des données:"
du -sh "$input_dir"
echo ""

echo "Nombre de reads par fichier (approximatif):"
for sample in "${sample_list[@]}"; do
    reads1=$(($(zcat "$input_dir/${sample}_1.fastq.gz" | wc -l) / 4))
    reads2=$(($(zcat "$input_dir/${sample}_2.fastq.gz" | wc -l) / 4))
    printf "  %-20s R1: %12d reads | R2: %12d reads\n" "$sample" "$reads1" "$reads2"
done
echo ""

# ====================================================================
# 4. LOG ET FINALISATION
# ====================================================================

echo "======================================================"
echo "RAPPORT FINAL"
echo "======================================================"
echo ""
echo "Fichiers générés:"
echo "  - Rapports FastQC: $fastqc_dir/"
echo "  - Rapport MultiQC: $multiqc_dir/multiqc_report.html"
echo ""
echo "Consulter le rapport:"
echo "  firefox $multiqc_dir/multiqc_report.html"
echo ""

# Log
echo "[$(date)] FastQC completed for all samples" >> "$log_dir/analysis.log"

echo "======================================================"
echo "✓ ÉTAPE 1 COMPLÉTÉE"
echo "======================================================"
echo ""
echo "PROCHAINES ÉTAPES:"
echo "1. Consulter le rapport MultiQC"
echo "2. Vérifier la qualité des reads"
echo "3. Lancer: bash 02_salmon_indexing.sh"
echo ""

