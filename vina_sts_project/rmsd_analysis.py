#!/usr/bin/env python3
"""
Validation du redocking 4MVF / staurosporine.
RMSD (atomes lourds) de chaque pose Vina vs ligand cristallin.
Perception des liaisons natives par distances covalentes (indépendante de CONECT).
"""
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Geometry import Point3D

STU_SMILES = "CN1C(=O)C=C(C2=C1C(=C(C3=C2NC4=C3C=C(C=C4)OC)C)C)C5=CC=CC=C5"
COV_RADII = {"C": 0.76, "N": 0.71, "O": 0.66}
TOL = 0.45  # tolérance (Å) au-delà de la somme des rayons covalents

subprocess.run(["mk_export.py", "results/poses.pdbqt", "-s", "results/poses.sdf"], check=True)

# 1. Lecture directe des atomes lourds du STU natif
atoms = []
with open("data/4mvf.pdb") as fin:
    for line in fin:
        if line.startswith("HETATM") and line[17:20].strip() == "STU":
            sym = line[76:78].strip() or line[12:13]
            x, y, z = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            atoms.append((sym, x, y, z))

if not atoms:
    raise SystemExit("Ligand natif STU introuvable dans data/4mvf.pdb")

# 2. RWMol avec liaisons perçues par distances covalentes (méthode déterministe)
raw = Chem.RWMol()
coords = []
for sym, x, y, z in atoms:
    raw.AddAtom(Chem.Atom(sym))
    coords.append((x, y, z))
pos = np.array(coords)
for i in range(len(atoms)):
    for j in range(i + 1, len(atoms)):
        if np.linalg.norm(pos[i] - pos[j]) <= COV_RADII[atoms[i][0]] + COV_RADII[atoms[j][0]] + TOL:
            raw.AddBond(i, j, Chem.BondType.SINGLE)

conf = Chem.Conformer(raw.GetNumAtoms())
for i, (x, y, z) in enumerate(coords):
    conf.SetAtomPosition(i, Point3D(x, y, z))
raw.AddConformer(conf)

# 3. Appariement topologique insensible aux H / ordres de liaison
ref = Chem.MolFromSmiles(STU_SMILES)
query = Chem.RWMol()
qidx = {}
for a in ref.GetAtoms():
    qidx[a.GetIdx()] = query.AddAtom(Chem.AtomFromSmarts(f"[#{a.GetAtomicNum()}]"))
for b in ref.GetBonds():
    query.AddBond(qidx[b.GetBeginAtomIdx()], qidx[b.GetEndAtomIdx()], Chem.BondType.UNSPECIFIED)

match = raw.GetSubstructMatch(query)
if len(match) != ref.GetNumAtoms():
    raise SystemExit(
        f"Échec d'appariement ({len(match)}/{ref.GetNumAtoms()} atomes). "
        f"raw: {raw.GetNumAtoms()} atomes, {raw.GetNumBonds()} liaisons."
    )

# 4. Injection des coordonnées natives dans le graphe de référence
ref3d = Chem.Mol(ref)
conf_ref = Chem.Conformer(ref.GetNumAtoms())
for r in range(ref.GetNumAtoms()):
    x, y, z = coords[match[r]]
    conf_ref.SetAtomPosition(r, Point3D(x, y, z))
ref3d.AddConformer(conf_ref)

# 5. RMSD par pose (superposition optimale sur atomes lourds appariés)
suppl = Chem.SDMolSupplier("results/poses.sdf", removeHs=True)
results = []
for i, pose in enumerate(suppl, start=1):
    if pose is None:
        continue
    ref2pose = pose.GetSubstructMatch(ref)
    if not ref2pose:
        continue
    atom_map = [(ref2pose[r], r) for r in range(ref.GetNumAtoms())]
    rmsd = rdMolAlign.AlignMol(pose, ref3d, atomMap=atom_map)
    results.append((i, rmsd))

# 6. Rapport
print(f"{'Pose':>4} | {'RMSD (Å)':>8}")
print("-" * 28)
for i, rmsd in results:
    flag = "  <-- succès (< 2.0 Å)" if rmsd < 2.0 else ""
    print(f"{i:>4} | {rmsd:>8.2f}{flag}")

best = min(results, key=lambda x: x[1])
print(f"\nMeilleure reproduction : pose {best[0]} (RMSD = {best[1]:.2f} Å)")
