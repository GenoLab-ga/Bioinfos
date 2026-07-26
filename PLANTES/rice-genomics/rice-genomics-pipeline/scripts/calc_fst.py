#!/usr/bin/env python3
"""Calcul du Fst de Weir & Cockerham entre groupes variétaux — Rice Genomics Pipeline"""
import sys, numpy as np
from itertools import combinations

ANALYSIS_DIR = sys.argv[1] if len(sys.argv) > 1 else "."

GROUPS = {
    "Japonica": [1, 2],
    "Basmati":  [3],
    "Aus/boro": [4, 5, 6, 7],
    "Indica":   [0],
}

print("Chargement des génotypes...")
genotypes = [[] for _ in range(8)]

with open(f"{ANALYSIS_DIR}/genotypes_raw.txt") as f:
    for i, line in enumerate(f):
        if i % 500000 == 0 and i > 0:
            print(f"  {i} variants lus...")
        fields = line.strip().split('\t')
        for j, gt in enumerate(fields[4:]):
            if gt in ('./.', '.|.', '.'):   genotypes[j].append(-1)
            elif gt in ('0/0', '0|0'):      genotypes[j].append(0)
            elif gt in ('0/1','1/0','0|1','1|0'): genotypes[j].append(1)
            elif gt in ('1/1', '1|1'):      genotypes[j].append(2)
            else:                           genotypes[j].append(-1)

G = np.array(genotypes, dtype=np.int8)
print(f"Variants chargés : {G.shape[1]}")

def calc_fst(idx1, idx2, G):
    """Fst simplifié : (Ht - Hs) / Ht"""
    all_idx = idx1 + idx2
    results = []
    chunk = 100000
    for start in range(0, G.shape[1], chunk):
        g = G[:, start:start+chunk].astype(np.float32)
        mask1 = np.all(g[idx1, :] >= 0, axis=0)
        mask2 = np.all(g[idx2, :] >= 0, axis=0)
        mask = mask1 & mask2
        if mask.sum() == 0: continue
        g = g[:, mask]
        p1 = g[idx1, :].mean(axis=0) / 2.0
        p2 = g[idx2, :].mean(axis=0) / 2.0
        pt = g[all_idx, :].mean(axis=0) / 2.0
        n1, n2 = len(idx1), len(idx2)
        Ht = 2 * pt * (1 - pt)
        Hs = (n1 * 2*p1*(1-p1) + n2 * 2*p2*(1-p2)) / (n1+n2)
        valid = Ht > 0.01
        if valid.sum() == 0: continue
        fst_vals = np.clip((Ht[valid] - Hs[valid]) / Ht[valid], 0, 1)
        results.extend(fst_vals.tolist())
    return (float(np.mean(results)), len(results)) if results else (0.0, 0)

print("\nCalcul du Fst...")
print("="*55)
with open(f"{ANALYSIS_DIR}/fst_results.tsv", "w") as out:
    out.write("Group1\tGroup2\tFst\n")
    for g1, g2 in combinations(GROUPS.keys(), 2):
        fst, n = calc_fst(GROUPS[g1], GROUPS[g2], G)
        print(f"Fst {g1:12s} vs {g2:12s} = {fst:.4f}  (n={n})")
        out.write(f"{g1}\t{g2}\t{fst:.4f}\n")

print(f"\nFst sauvegardé : {ANALYSIS_DIR}/fst_results.tsv")
