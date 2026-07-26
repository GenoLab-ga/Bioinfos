import json
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from collections import Counter

results_dir = 'data/results/tbprofiler_v6/results'

# Top mutations à afficher
top_mutations = [
    'katG:p.Ser315Thr', 'rpoB:p.Ser450Leu', 'rpsL:p.Lys43Arg',
    'inhA:c.-777C>T', 'gyrA:p.Asp94Gly', 'embB:p.Met306Val',
    'gyrA:p.Ala90Val', 'embB:p.Met306Ile', 'rrs:n.1401A>G',
    'embB:p.Gln497Arg', 'rplC:p.Cys154Arg', 'rpoB:p.Leu452Pro',
    'gyrA:p.Asp94Ala', 'rpoB:p.Leu430Pro', 'gyrA:p.Asp94Asn',
    'rpoB:p.His445Asn', 'pncA:p.Gly132Ala', 'gyrA:p.Asp94Tyr',
    'eis:c.-10G>A', 'embB:p.Gly406Ala'
]

# Compter par drtype
drtype_mutations = {}
drtype_order = ['Susceptible','HR-TB','RR-TB','MDR-TB','Pre-XDR-TB','XDR-TB']

for drtype in drtype_order:
    drtype_mutations[drtype] = Counter()

for fname in os.listdir(results_dir):
    if not fname.endswith('.results.json'):
        continue
    with open(f'{results_dir}/{fname}') as f:
        data = json.load(f)
    drtype = data.get('drtype','Other')
    if drtype not in drtype_mutations:
        continue
    for v in data.get('dr_variants', []):
        gene = v.get('gene_name') or 'N/A'
        change = v.get('change') or 'N/A'
        mut = f"{gene}:{change}"
        if mut in top_mutations:
            drtype_mutations[drtype][mut] += 1

# Matrice pour heatmap
matrix = np.zeros((len(top_mutations), len(drtype_order)))
for j, drtype in enumerate(drtype_order):
    for i, mut in enumerate(top_mutations):
        matrix[i, j] = drtype_mutations[drtype].get(mut, 0)

fig, ax = plt.subplots(figsize=(12, 10))
im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto')
ax.set_xticks(range(len(drtype_order)))
ax.set_xticklabels(drtype_order, rotation=30, ha='right', fontsize=11)
ax.set_yticks(range(len(top_mutations)))
ax.set_yticklabels(top_mutations, fontsize=9)
plt.colorbar(im, ax=ax, label='Nombre de samples')

for i in range(len(top_mutations)):
    for j in range(len(drtype_order)):
        val = int(matrix[i, j])
        if val > 0:
            color = 'white' if val > matrix.max()*0.6 else 'black'
            ax.text(j, i, str(val), ha='center', va='center',
                   fontsize=9, color=color, fontweight='bold')

ax.set_title('Distribution des mutations de résistance par profil clinique\n'
             'PRJNA1177198 — MTB Nanopore Targeted Sequencing (n=277)',
             fontsize=13, fontweight='bold', pad=15)
ax.set_xlabel('Classification de résistance (WHO 2021)', fontsize=11)
ax.set_ylabel('Mutation', fontsize=11)

plt.tight_layout()
plt.savefig('data/results/tbprofiler_v6/summary/mutation_heatmap.png',
            dpi=150, bbox_inches='tight')
print("Heatmap sauvegardée : mutation_heatmap.png")
