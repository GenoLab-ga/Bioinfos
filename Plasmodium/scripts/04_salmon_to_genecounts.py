#!/usr/bin/env python3

"""
====================================================================
CONVERSION TRANSCRIPTS → GÈNES
De: quant.sf (transcripts Salmon) → Matrice gènes × échantillons
====================================================================
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import json

# ====================================================================
# 0. VARIABLES ET CHEMINS
# ====================================================================

PROJECT_ROOT = Path.home() / "Documents/github_projet/Plasmodium"
salmon_dir = PROJECT_ROOT / "quantification" / "salmon_output"
output_dir = PROJECT_ROOT / "quantification" / "matrices"
reference_dir = PROJECT_ROOT / "reference"
transcriptome_file = reference_dir / "Pfalciparum3D7_transcripts.fasta"
log_dir = PROJECT_ROOT / "logs"

# Échantillons
samples = [
    "ERR11471971",  # Day_1_pf_gam_tc_v2 (CONTRÔLE)
    "ERR11471975",  # Unknown
    "ERR11471979",  # Unknown
    "ERR11471985"   # Unknown
]

# Créer répertoires
output_dir.mkdir(parents=True, exist_ok=True)
log_dir.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("CONVERSION TRANSCRIPTS → GÈNES")
print("=" * 60)
print("")

# ====================================================================
# 1. VÉRIFICATIONS PRÉALABLES
# ====================================================================

print("=" * 60)
print("ÉTAPE 1: Vérifications préalables")
print("=" * 60)
print("")

# Vérifier que Salmon outputs existent
print("Vérification des fichiers Salmon:")
for sample in samples:
    quant_file = salmon_dir / sample / "quant.sf"
    if quant_file.exists():
        num_lines = sum(1 for _ in open(quant_file)) - 1  # -1 pour header
        print(f"  ✓ {sample}: {num_lines} transcripts")
    else:
        print(f"  ❌ {sample}: fichier quant.sf manquant")
        sys.exit(1)

print("")
print("Chemins:")
print(f"  Salmon outputs: {salmon_dir}")
print(f"  Output matrices: {output_dir}")
print("")

# ====================================================================
# 2. CHARGER MAPPING TRANSCRIPT → GÈNE
# ====================================================================

print("=" * 60)
print("ÉTAPE 2: Création mapping transcript → gène")
print("=" * 60)
print("")

# Parser le FASTA pour créer mapping
transcript_to_gene = {}

if transcriptome_file.exists():
    print(f"Parsing: {transcriptome_file}")
    with open(transcriptome_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Format Pf: >PF3D7_0100100|PF3D7_0100100.1 ...
                header = line.strip('>\n')
                parts = header.split('|')
                
                if len(parts) >= 2:
                    gene_id = parts[0]
                    transcript_id = parts[1].split()[0]
                else:
                    # Fallback
                    transcript_id = header.split()[0]
                    gene_id = transcript_id
                
                transcript_to_gene[transcript_id] = gene_id
    
    print(f"  ✓ {len(transcript_to_gene)} transcripts mappés")
else:
    print(f"  ⚠ Fichier FASTA non trouvé: {transcriptome_file}")
    print("  Utilisation IDs transcripts comme gènes")

print("")

# ====================================================================
# 3. CHARGER ET AGRÉGER PAR ÉCHANTILLON
# ====================================================================

print("=" * 60)
print("ÉTAPE 3: Agrégation transcripts → gènes")
print("=" * 60)
print("")

# Dictionnaire pour stocker les données
all_matrices = {}

for sample in samples:
    print(f"Traitement: {sample}")
    
    # Charger quant.sf
    quant_file = salmon_dir / sample / "quant.sf"
    quant_df = pd.read_csv(quant_file, sep='\t')
    
    print(f"  Reads total: {quant_df['NumReads'].sum():.0f}")
    
    # Ajouter mapping gène
    quant_df['gene_id'] = quant_df['Name'].map(transcript_to_gene)
    
    # Compter transcripts non-mappés
    unmapped = quant_df['gene_id'].isna().sum()
    if unmapped > 0:
        print(f"  ⚠ {unmapped} transcripts sans mapping (supprimés)")
        quant_df = quant_df.dropna(subset=['gene_id'])
    
    # Agréger par gène (somme des counts)
    gene_counts = quant_df.groupby('gene_id')['NumReads'].sum().reset_index()
    gene_counts.columns = ['gene_id', 'counts']
    
    print(f"  Gènes détectés: {len(gene_counts)}")
    print(f"  Counts total: {gene_counts['counts'].sum():.0f}")
    print("")
    
    # Sauvegarder
    output_file = output_dir / f"{sample}_gene_counts.csv"
    gene_counts.to_csv(output_file, index=False)
    
    all_matrices[sample] = gene_counts

print("")

# ====================================================================
# 4. CRÉER MATRICE COMBINÉE
# ====================================================================

print("=" * 60)
print("ÉTAPE 4: Création matrice combinée (gènes × échantillons)")
print("=" * 60)
print("")

# Fusionner toutes les matrices
combined = None
for sample, gene_counts in all_matrices.items():
    gene_counts = gene_counts.set_index('gene_id')
    gene_counts.columns = [sample]
    
    if combined is None:
        combined = gene_counts.copy()
    else:
        combined = combined.join(gene_counts, how='outer')

# Remplir NaN par 0 (gène non-détecté)
combined = combined.fillna(0).astype(int)

# Statistiques
num_genes = combined.shape[0]
num_samples = combined.shape[1]
sparsity = (combined == 0).sum().sum() / (num_genes * num_samples) * 100
total_counts = combined.sum().sum()

print(f"Matrice combinée:")
print(f"  Dimensions: {num_genes} gènes × {num_samples} échantillons")
print(f"  Sparsité: {sparsity:.1f}%")
print(f"  Counts total: {total_counts:.0f}")
print("")

# Sauvegarder
matrix_file = output_dir / "combined_counts_matrix.csv"
combined.to_csv(matrix_file)
print(f"✓ Matrice sauvegardée: {matrix_file}")

# Créer preview (premiers 500 gènes, tous les échantillons)
try:
    import openpyxl
    preview_file = output_dir / "combined_counts_matrix_preview.xlsx"
    combined.iloc[:500, :].to_excel(preview_file)
    print(f"✓ Preview Excel: {preview_file}")
except:
    print("  (openpyxl non disponible, skip preview Excel)")

print("")

# ====================================================================
# 5. STATISTIQUES PAR ÉCHANTILLON
# ====================================================================

print("=" * 60)
print("ÉTAPE 5: Statistiques par échantillon")
print("=" * 60)
print("")

for sample in samples:
    total_counts = combined[sample].sum()
    num_detected = (combined[sample] > 0).sum()
    mean_counts = combined[combined[sample] > 0][sample].mean()
    
    print(f"{sample}:")
    print(f"  Counts total: {total_counts:.0f}")
    print(f"  Gènes détectés: {num_detected}")
    print(f"  Mean counts/gène: {mean_counts:.1f}")

print("")

# ====================================================================
# 6. LOG ET FINALISATION
# ====================================================================

print("=" * 60)
print("RAPPORT FINAL")
print("=" * 60)
print("")

print("✓ Conversion complétée")
print("")
print(f"Fichiers générés dans: {output_dir}")
print(f"  - combined_counts_matrix.csv (matrice principale)")
print(f"  - *_gene_counts.csv (par échantillon)")
print("")
print(f"Dimensions finales: {num_genes} gènes × {num_samples} échantillons")
print("")

# Log
with open(log_dir / "analysis.log", "a") as f:
    f.write(f"[{pd.Timestamp.now()}] Transcript-to-gene conversion completed\n")

print("=" * 60)
print("✓ ÉTAPE 4 COMPLÉTÉE")
print("=" * 60)
print("")
print("PROCHAINES ÉTAPES:")
print("1. Vérifier: head combined_counts_matrix.csv")
print("2. Lancer R: bash 05_seurat_qc.sh")
print("")

