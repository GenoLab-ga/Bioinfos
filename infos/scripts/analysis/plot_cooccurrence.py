import json
import os
import numpy as np
import matplotlib.pyplot as plt

results_dir = 'data/results/tbprofiler_v6/results'

top_genes = ['rpoB','katG','gyrA','embB','pncA','rpsL',
             'inhA','rrs','ethA','rplC','eis','gyrB']

# Matrice de co-occurrence
matrix = np.zeros((len(top_genes), len(top_genes)), dtype=int)
gene_counts = {g: 0 for g in top_genes}

for fname in os.listdir(results_dir):
    if not fname.endswith('.results.json'):
        continue
    with open(f'{results_dir}/{fname}') as f:
        data = json.load(f)
    genes_in_sample = set()
    for v in data.get('dr_variants', []):
        gene = v.get('gene_name')
        if gene in top_genes:
            genes_in_sample.add(gene)
    for g in genes_in_sample:
        gene_counts[g] += 1
    genes_list = list(genes_in_sample)
    for i, g1 in enumerate(genes_list):
        for g2 in genes_list[i:]:
            idx1 = top_genes.index(g1)
            idx2 = top_genes.index(g2)
            matrix[idx1][idx2] += 1
            if idx1 != idx2:
                matrix[idx2][idx1] += 1

# Masquer diagonale pour lisibilité
np.fill_diagonal(matrix, 0)

fig, ax = plt.subplots(figsize=(11, 9))
im = ax.imshow(matrix, cmap='Blues', aspect='auto')
ax.set_xticks(range(len(top_genes)))
ax.set_yticks(range(len(top_genes)))
ax.set_xticklabels(top_genes, rotation=45, ha='right', fontsize=11)
ax.set_yticklabels(top_genes, fontsize=11)
plt.colorbar(im, ax=ax, label='Samples co-porteurs')

for i in range(len(top_genes)):
    for j in range(len(top_genes)):
        val = matrix[i, j]
        if val > 0:
            color = 'white' if val > matrix.max()*0.6 else 'black'
            ax.text(j, i, str(val), ha='center', va='center',
                   fontsize=9, color=color, fontweight='bold')

ax.set_title('Co-occurrence des mutations de résistance par gène\n'
             'PRJNA1177198 — MTB Nanopore Targeted Sequencing (n=277)',
             fontsize=13, fontweight='bold', pad=15)

plt.tight_layout()
plt.savefig('data/results/tbprofiler_v6/summary/cooccurrence_matrix.png',
            dpi=150, bbox_inches='tight')
print("Matrice sauvegardée : cooccurrence_matrix.png")
