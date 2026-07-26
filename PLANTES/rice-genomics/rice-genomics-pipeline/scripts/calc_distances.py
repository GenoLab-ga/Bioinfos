#!/usr/bin/env python3
"""Calcul de la matrice de distances IBS — Rice Genomics Pipeline"""
import sys, numpy as np

ANALYSIS_DIR = sys.argv[1] if len(sys.argv) > 1 else "."

SAMPLE_NAMES = [
    "SAMEA2567482","SAMEA2567483","SAMEA2567484","SAMEA2567493",
    "SAMEA2567500","SAMEA2567501","SAMEA2567504","SAMEA2567509"
]
LABELS = {
    "SAMEA2567482": "Brazil_Indica",
    "SAMEA2567483": "SouthKorea_TempJaponica",
    "SAMEA2567484": "IvoryCoast_TropJaponica",
    "SAMEA2567493": "India_Basmati",
    "SAMEA2567500": "Bangladesh_Aus1",
    "SAMEA2567501": "India_Aus1",
    "SAMEA2567504": "India_Aus2",
    "SAMEA2567509": "Bangladesh_Aus2",
}

print("Chargement des génotypes...")
genotypes = [[] for _ in range(8)]

with open(f"{ANALYSIS_DIR}/genotypes_raw.txt") as f:
    for i, line in enumerate(f):
        if i % 500000 == 0 and i > 0:
            print(f"  {i} variants lus...")
        fields = line.strip().split('\t')
        gts = fields[4:]
        for j, gt in enumerate(gts):
            if gt in ('./.', '.|.', '.'):
                genotypes[j].append(-1)
            elif gt in ('0/0', '0|0'):
                genotypes[j].append(0)
            elif gt in ('0/1', '1/0', '0|1', '1|0'):
                genotypes[j].append(1)
            elif gt in ('1/1', '1|1'):
                genotypes[j].append(2)
            else:
                genotypes[j].append(-1)

print(f"Variants chargés : {len(genotypes[0])}")

n = len(SAMPLE_NAMES)
G = np.array(genotypes, dtype=np.int8)
dist_matrix = np.zeros((n, n))

print("Calcul des distances IBS...")
for i in range(n):
    for j in range(i+1, n):
        gi, gj = G[i].astype(np.float32), G[j].astype(np.float32)
        mask = (gi >= 0) & (gj >= 0)
        dist = np.sum(gi[mask] != gj[mask]) / mask.sum() if mask.sum() > 0 else 1.0
        dist_matrix[i,j] = dist_matrix[j,i] = dist
        print(f"  {SAMPLE_NAMES[i]} vs {SAMPLE_NAMES[j]} : {dist:.4f}")

label_list = [LABELS[s] for s in SAMPLE_NAMES]
with open(f"{ANALYSIS_DIR}/distance_matrix.tsv", "w") as out:
    out.write("\t" + "\t".join(label_list) + "\n")
    for i, s in enumerate(SAMPLE_NAMES):
        row = "\t".join([f"{dist_matrix[i,j]:.4f}" for j in range(n)])
        out.write(f"{LABELS[s]}\t{row}\n")

print(f"Matrice sauvegardée : {ANALYSIS_DIR}/distance_matrix.tsv")
