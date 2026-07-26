#!/bin/bash

# ====================================================================
# PIPELINE SINGLE-CELL RNA-SEQ: Plasmodium falciparum Gamétocytes
# PRJEB55754
# ====================================================================

# ====================================================================
# 0. VARIABLES ET CHEMINS
# ====================================================================

PROJECT_ROOT=~/Documents/github_projet/Plasmodium
reference_dir="$PROJECT_ROOT/reference"
transcriptome_fasta="$reference_dir/Pfalciparum3D7_transcripts.fasta"
salmon_index_dir="$reference_dir/Pfalciparum3D7_salmon_index"
genome_fasta="$reference_dir/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"
gff_file="$reference_dir/PlasmoDB-68_Pfalciparum3D7.gff"
log_dir="$PROJECT_ROOT/logs"

# Paramètres Salmon
SALMON_KMER=31
SALMON_THREADS=12

# URLs pour téléchargement (utilisé seulement si transcriptome absent)
PF_VERSION="68"
TRANSCRIPT_URL="https://plasmodb.org/common/downloads/release-${PF_VERSION}/Pfalciparum3D7/fasta/data/PlasmoDB-${PF_VERSION}_Pfalciparum3D7_AnnotatedTranscripts.fasta"

# Créer répertoires
mkdir -p "$reference_dir" "$log_dir"

echo "======================================================"
echo "INDEXATION GÉNOME: Plasmodium falciparum"
echo "======================================================"
echo ""

# ====================================================================
# 1. TÉLÉCHARGER LE TRANSCRIPTOME
# ====================================================================

echo "======================================================"
echo "ÉTAPE 1: Vérification/Téléchargement du transcriptome Pf"
echo "======================================================"
echo ""

echo "Fichiers de référence utilisateur détectés:"
if [ -f "$genome_fasta" ]; then
    genome_size=$(du -h "$genome_fasta" | cut -f1)
    echo "  ✓ Génome: $(basename $genome_fasta) ($genome_size)"
else
    echo "  ⚠ Génome non trouvé: $genome_fasta"
fi

if [ -f "$gff_file" ]; then
    gff_size=$(du -h "$gff_file" | cut -f1)
    echo "  ✓ Annotation GFF: $(basename $gff_file) ($gff_size)"
else
    echo "  ⚠ GFF non trouvé: $gff_file"
fi

echo ""

if [ ! -f "$transcriptome_fasta" ]; then
    echo "Transcriptome non trouvé localement. Téléchargement..."
    echo "Source: PlasmoDB version $PF_VERSION"
    echo ""
    
    echo "URL: $TRANSCRIPT_URL"
    echo ""
    
    cd "$reference_dir"
    
    wget -q --show-progress "$TRANSCRIPT_URL" -O "Pfalciparum3D7_transcripts.fasta"
    
    if [ $? -eq 0 ]; then
        echo "✓ Téléchargement réussi"
    else
        echo "❌ Erreur lors du téléchargement"
        echo "Téléchargez manuellement: $TRANSCRIPT_URL"
        exit 1
    fi
    
    cd - > /dev/null
else
    echo "✓ Transcriptome trouvé localement (réutilisation)"
fi

echo ""
echo "Informations transcriptome:"
file_size=$(du -h "$transcriptome_fasta" | cut -f1)
num_transcripts=$(grep -c "^>" "$transcriptome_fasta")
echo "  Fichier: $(basename $transcriptome_fasta)"
echo "  Taille: $file_size"
echo "  Transcripts: $num_transcripts"
echo ""

# ====================================================================
# 2. VÉRIFIER QUE SALMON EST DISPONIBLE
# ====================================================================

echo "======================================================"
echo "ÉTAPE 2: Vérification de Salmon"
echo "======================================================"
echo ""

if ! command -v salmon &> /dev/null; then
    echo "❌ Erreur: Salmon non trouvé"
    echo "Solution: mamba install salmon"
    exit 1
fi

echo "Version Salmon: $(salmon --version)"
echo "Threads disponibles: $SALMON_THREADS"
echo ""

# ====================================================================
# 3. CRÉER L'INDEX SALMON
# ====================================================================

echo "======================================================"
echo "ÉTAPE 3: Construction de l'index Salmon"
echo "======================================================"
echo ""

if [ -d "$salmon_index_dir" ]; then
    echo "⚠ Index déjà existe: $salmon_index_dir"
    echo "Suppression et recréation..."
    rm -rf "$salmon_index_dir"
fi

echo "Paramètres d'indexation:"
echo "  k-mer size: $SALMON_KMER"
echo "  Threads: $SALMON_THREADS"
echo "  Décoys: Inclus (améliore mapping)"
echo ""

echo "Construction en cours (peut prendre 5-10 minutes)..."
echo ""

salmon index -t "$transcriptome_fasta" \
             -i "$salmon_index_dir" \
             -k "$SALMON_KMER" \
             -p "$SALMON_THREADS" \
             2>&1 | tee "$log_dir/salmon_index.log"

if [ -d "$salmon_index_dir" ]; then
    echo ""
    echo "✓ Index créé avec succès"
else
    echo ""
    echo "❌ Erreur: Index non créé"
    exit 1
fi

echo ""

# ====================================================================
# 4. VÉRIFICATION INDEX
# ====================================================================

echo "======================================================"
echo "ÉTAPE 4: Vérification de l'index"
echo "======================================================"
echo ""

echo "Contenu de l'index:"
ls -lh "$salmon_index_dir"/ | grep -E "^d|json"
echo ""

index_size=$(du -sh "$salmon_index_dir" | cut -f1)
echo "Taille de l'index: $index_size"
echo ""

# ====================================================================
# 5. LOG ET FINALISATION
# ====================================================================

echo "======================================================"
echo "RAPPORT FINAL"
echo "======================================================"
echo ""
echo "✓ Index Salmon prêt pour quantification"
echo ""
echo "Localisation: $salmon_index_dir"
echo ""
echo "Fichiers importants:"
echo "  - Transcriptome: $transcriptome_fasta"
echo "  - Index: $salmon_index_dir/"
echo "  - Log: $log_dir/salmon_index.log"
echo ""

# Log
echo "[$(date)] Salmon index created successfully" >> "$log_dir/analysis.log"

echo "======================================================"
echo "✓ ÉTAPE 2 COMPLÉTÉE"
echo "======================================================"
echo ""
echo "PROCHAINES ÉTAPES:"
echo "1. Vérifier rapport FastQC (étape 1)"
echo "2. Si OK → Lancer: bash 03_salmon_quantification.sh"
echo ""

