#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="vina_sts_project"
mkdir -p "$PROJECT_DIR"/{data,results,logs}
cd "$PROJECT_DIR"

# 1. Téléchargement de la structure 4MVF
echo "[INFO] Téléchargement de la structure 4MVF..."
wget -q -O data/4mvf.pdb "https://files.rcsb.org/download/4MVF.pdb"

# 2. Nettoyage du PDB (suppression des HETATM : eaux, ligand natif, ions)
echo "[INFO] Nettoyage du PDB (suppression des hétéroatomes)..."
pdb_delhetatm data/4mvf.pdb > data/4mvf_clean.pdb

# 3. Préparation du récepteur avec meeko (gestion automatique des erreurs et altlocs)
echo "[INFO] Préparation du récepteur..."
mk_prepare_receptor.py -i data/4mvf_clean.pdb -p data/receptor.pdbqt -a --default_altloc A

# 4. Préparation du ligand (Staurosporine) via RDKit et Meeko
echo "[INFO] Préparation du ligand..."
python3 << 'EOF'
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess

# SMILES de la staurosporine
smiles = "CN1C(=O)C=C(C2=C1C(=C(C3=C2NC4=C3C=C(C=C4)OC)C)C)C5=CC=CC=C5"
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol)

# Export au format MOL requis par meeko pour la perception des liaisons
Chem.MolToMolFile(mol, "data/ligand_raw.mol")

# Conversion en PDBQT avec meeko
subprocess.run(["mk_prepare_ligand.py", "-i", "data/ligand_raw.mol", "-o", "data/ligand.pdbqt"], check=True)
print("[INFO] Ligand préparé avec succès.")
EOF

# 5. Configuration de la boîte de docking
cat > data/config.txt << 'EOF'
receptor = data/receptor.pdbqt
ligand = data/ligand.pdbqt
center_x = 32.055
center_y = 85.045
center_z = 22.125
size_x = 20.0
size_y = 20.0
size_z = 20.0
exhaustiveness = 16
cpu = 6
num_modes = 9
energy_range = 3
EOF

# 6. Exécution de Vina
echo "[INFO] Lancement du docking avec AutoDock Vina (6 cœurs)..."
vina --config data/config.txt --out results/poses.pdbqt | tee logs/vina.log

echo "[INFO] Docking terminé. Résultats dans $PROJECT_DIR/results/"
