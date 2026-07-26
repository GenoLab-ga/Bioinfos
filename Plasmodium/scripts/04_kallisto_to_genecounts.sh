#!/bin/bash

# ====================================================================
# CRÉER LES MATRICES: À partir de Kallisto outputs
# Kallisto crée matrix.ec + matrix.cells + matrix.genes
# On crée une matrice combinée (gènes × échantillons)
# ====================================================================

PROJECT_ROOT=~/Documents/github_projet/Plasmodium
kallisto_output="$PROJECT_ROOT/quantification/kallisto_output"
reference_dir="$PROJECT_ROOT/reference"
transcriptome_fasta="$reference_dir/Pfalciparum3D7_transcripts.fasta"
matrices_dir="$PROJECT_ROOT/quantification/matrices"
log_dir="$PROJECT_ROOT/logs"

sample_list=(
    ERR11471971
    ERR11471975
    ERR11471979
    ERR11471985
)

mkdir -p "$matrices_dir" "$log_dir"

echo "======================================================"
echo "CRÉATION DES MATRICES: Kallisto → Gènes"
echo "======================================================"
echo ""

# ====================================================================
# 1. VÉRIFICATIONS
# ====================================================================

echo "======================================================"
echo "ÉTAPE 1: Vérifications"
echo "======================================================"
echo ""

# Vérifier les outputs Kallisto
echo "Vérification outputs Kallisto:"
for sample in "${sample_list[@]}"; do
    sample_output="$kallisto_output/$sample"
    if [ -f "$sample_output/abundance.tsv" ]; then
        echo "  ✓ $sample: abundance.tsv trouvé"
    else
        echo "  ❌ $sample: fichier manquant"
        exit 1
    fi
done

echo ""

# ====================================================================
# 2. TRAITER CHAQUE ÉCHANTILLON
# ====================================================================

echo "======================================================"
echo "ÉTAPE 2: Créer les matrices par échantillon"
echo "======================================================"
echo ""

# Python inline pour traiter les fichiers Kallisto
python3 << 'PYTHON_EOF'
import pandas as pd
import numpy as np
from pathlib import Path

project_root = Path.home() / "Documents/github_projet/Plasmodium"
kallisto_output = project_root / "quantification/kallisto_output"
matrices_dir = project_root / "quantification/matrices"
transcriptome_file = project_root / "reference/Pfalciparum3D7_transcripts.fasta"

samples = ["ERR11471971", "ERR11471975", "ERR11471979", "ERR11471985"]

print("Traitement des outputs Kallisto...")
print("")

# Créer mapping transcript → gène
transcript_to_gene = {}
if transcriptome_file.exists():
    with open(transcriptome_file) as f:
        for line in f:
            if line.startswith('>'):
                header = line.strip('>\n')
                parts = header.split('|')
                if len(parts) >= 2:
                    gene_id = parts[0]
                    transcript_id = parts[1].split()[0]
                    transcript_to_gene[transcript_id] = gene_id

print(f"✓ Mapping: {len(transcript_to_gene)} transcripts")
print("")

# Charger et agréger par échantillon
all_matrices = {}

for sample in samples:
    sample_output = kallisto_output / sample
    abundance_file = sample_output / "abundance.tsv"
    
    if not abundance_file.exists():
        print(f"⚠ {sample}: fichier manquant")
        continue
    
    print(f"Traitement: {sample}")
    
    # Charger abundance.tsv (kallisto output)
    # Colonnes: target_id, length, eff_length, est_counts, tpm
    df = pd.read_csv(abundance_file, sep='\t')
    
    # Mapper transcripts → gènes
    df['gene_id'] = df['target_id'].map(transcript_to_gene)
    
    # Agréger par gène (sum est_counts)
    gene_counts = df.groupby('gene_id')['est_counts'].sum().reset_index()
    gene_counts.columns = ['gene_id', 'counts']
    gene_counts['counts'] = gene_counts['counts'].astype(int)
    
    print(f"  Gènes: {len(gene_counts)}")
    print(f"  Counts total: {gene_counts['counts'].sum():.0f}")
    print("")
    
    # Sauvegarder
    output_file = matrices_dir / f"{sample}_gene_counts.csv"
    gene_counts.to_csv(output_file, index=False)
    
    all_matrices[sample] = gene_counts

print("")

# ====================================================================
# 3. CRÉER MATRICE COMBINÉE
# ====================================================================

print("====================================================")
print("ÉTAPE 3: Matrice combinée (gènes × échantillons)")
print("====================================================")
print("")

if not all_matrices:
    print("❌ Aucune matrice créée")
    exit(1)

# Fusionner
combined = None
for sample, gene_counts in all_matrices.items():
    gc = gene_counts.set_index('gene_id')
    gc.columns = [sample]
    
    if combined is None:
        combined = gc
    else:
        combined = combined.join(gc, how='outer')

# Remplir NaN par 0
combined = combined.fillna(0).astype(int)

print(f"Dimensions: {combined.shape[0]} gènes × {combined.shape[1]} échantillons")
print(f"Sparsité: {(combined == 0).sum().sum() / (combined.shape[0] * combined.shape[1]) * 100:.1f}%")
print("")

# Sauvegarder
output_file = matrices_dir / "combined_counts_matrix.csv"
combined.to_csv(output_file)
print(f"✓ Matrice sauvegardée: {output_file}")

print("")
print("====================================================")
print("✓ MATRICES CRÉÉES AVEC SUCCÈS")
print("====================================================")

PYTHON_EOF

echo ""

# ====================================================================
# 4. FINALISATION
# ====================================================================

echo "======================================================"
echo "RÉSUMÉ FINAL"
echo "======================================================"
echo ""

if [ -f "$matrices_dir/combined_counts_matrix.csv" ]; then
    echo "✓ Matrice combinée créée"
    
    # Afficher dimensions
    lines=$(wc -l < "$matrices_dir/combined_counts_matrix.csv")
    cols=$(head -1 "$matrices_dir/combined_counts_matrix.csv" | tr ',' '\n' | wc -l)
    
    echo "  Fichier: combined_counts_matrix.csv"
    echo "  Dimensions: $((lines-1)) gènes × $((cols-1)) échantillons"
else
    echo "❌ Erreur: matrice non créée"
    exit 1
fi

echo ""
echo "[$(date)] Kallisto to gene matrix conversion completed" >> "$log_dir/analysis.log"

echo "======================================================"
echo "✓ ÉTAPE CONVERSION COMPLÉTÉE"
echo "======================================================"
echo ""
echo "PROCHAINES ÉTAPES:"
echo "1. Lancer Seurat QC: bash 05_seurat_qc.sh"
echo ""

