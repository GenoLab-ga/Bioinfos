#!/bin/bash

# ====================================================================
# INDEXATION: Kallisto (meilleur pour single-cell que Salmon)
# ====================================================================

PROJECT_ROOT=~/Documents/github_projet/Plasmodium
reference_dir="$PROJECT_ROOT/reference"
transcriptome_fasta="$reference_dir/Pfalciparum3D7_transcripts.fasta"
kallisto_index="$reference_dir/Pfalciparum3D7_kallisto_index.idx"
log_dir="$PROJECT_ROOT/logs"

mkdir -p "$reference_dir" "$log_dir"

echo "======================================================"
echo "INDEXATION: Kallisto"
echo "======================================================"
echo ""

# ====================================================================
# 1. VÉRIFICATIONS
# ====================================================================

echo "======================================================"
echo "ÉTAPE 1: Vérifications préalables"
echo "======================================================"
echo ""

if ! command -v kallisto &> /dev/null; then
    echo "❌ Erreur: kallisto non trouvé"
    echo "Installation: mamba install kallisto"
    exit 1
fi

echo "Version kallisto: $(kallisto version)"
echo ""

if [ ! -f "$transcriptome_fasta" ]; then
    echo "❌ Transcriptome non trouvé: $transcriptome_fasta"
    exit 1
fi

echo "Transcriptome trouvé: $(basename $transcriptome_fasta)"
echo "  Taille: $(du -h $transcriptome_fasta | cut -f1)"
echo "  Transcripts: $(grep -c '^>' $transcriptome_fasta)"
echo ""

# ====================================================================
# 2. CRÉER L'INDEX KALLISTO
# ====================================================================

echo "======================================================"
echo "ÉTAPE 2: Construction de l'index Kallisto"
echo "======================================================"
echo ""

if [ -f "$kallisto_index" ]; then
    echo "⚠ Index déjà existe. Recréation..."
    rm "$kallisto_index"
fi

echo "Construction en cours (peut prendre 2-3 minutes)..."
echo ""

kallisto index -i "$kallisto_index" \
               -k 31 \
               "$transcriptome_fasta" \
               2>&1 | tee "$log_dir/kallisto_index.log"

if [ -f "$kallisto_index" ]; then
    echo ""
    echo "✓ Index Kallisto créé avec succès"
    
    index_size=$(du -h "$kallisto_index" | cut -f1)
    echo "  Taille: $index_size"
    echo "  Localisation: $kallisto_index"
else
    echo "❌ Erreur: Index non créé"
    exit 1
fi

echo ""
echo "======================================================"
echo "✓ INDEXATION KALLISTO COMPLÉTÉE"
echo "======================================================"
echo ""

echo "[$(date)] Kallisto index created" >> "$log_dir/analysis.log"

