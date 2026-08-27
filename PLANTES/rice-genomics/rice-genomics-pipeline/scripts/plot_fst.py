#!/usr/bin/env python3
"""Visualisation Fst — Rice Genomics Pipeline"""
import sys, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

ANALYSIS_DIR = sys.argv[1] if len(sys.argv) > 1 else "."

groups = ["Japonica","Basmati","Aus/boro","Indica"]
n = len(groups)
fst_data = {}
with open(f"{ANALYSIS_DIR}/fst_results.tsv") as f:
    f.readline()
    for line in f:
        g1, g2, fst = line.strip().split('\t')
        fst_data[(g1, g2)] = float(fst)

matrix = np.zeros((n,n))
for (g1,g2), fst in fst_data.items():
    i, j = groups.index(g1), groups.index(g2)
    matrix[i,j] = matrix[j,i] = fst

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

ax1 = axes[0]
display = np.where(np.eye(n, dtype=bool), np.nan, matrix)
im = ax1.imshow(display, cmap='YlOrRd', vmin=0, vmax=0.65)
plt.colorbar(im, ax=ax1, label='Fst')
ax1.set_xticks(range(n)); ax1.set_yticks(range(n))
ax1.set_xticklabels(groups, fontsize=11, rotation=20)
ax1.set_yticklabels(groups, fontsize=11)
for i in range(n):
    for j in range(n):
        if i != j:
            v = matrix[i,j]
            ax1.text(j, i, f'{v:.3f}', ha='center', va='center',
                     fontsize=12, fontweight='bold',
                     color='white' if v > 0.45 else 'black')
ax1.set_title("Matrice Fst entre groupes variétaux\nOryza sativa — 3.3M SNPs", fontsize=12)

ax2 = axes[1]
pairs = [f"{g1}\nvs\n{g2}" for (g1,g2) in fst_data.keys()]
values = list(fst_data.values())
bar_colors = ['#4DAF4A' if v<0.30 else '#FF7F00' if v<0.45 else '#E41A1C' for v in values]
bars = ax2.barh(pairs, values, color=bar_colors, height=0.6)
for bar, val in zip(bars, values):
    ax2.text(val+0.005, bar.get_y()+bar.get_height()/2,
             f'{val:.4f}', va='center', fontsize=10, fontweight='bold')
ax2.set_xlabel("Fst", fontsize=12)
ax2.set_title("Différenciation génétique par paire\n(Wright 1978)", fontsize=12)
ax2.set_xlim(0, 0.70)
ax2.axvline(0.25, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax2.axvline(0.45, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax2.grid(axis='x', alpha=0.3)
legend = [
    mpatches.Patch(color='#4DAF4A', label='Fst < 0.30'),
    mpatches.Patch(color='#FF7F00', label='Fst 0.30–0.45'),
    mpatches.Patch(color='#E41A1C', label='Fst > 0.45'),
]
ax2.legend(handles=legend, fontsize=9, loc='lower right')

plt.suptitle("Différenciation génétique entre sous-populations du riz cultivé\n"
             "3000 Rice Genomes Project — 8 accessions",
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f"{ANALYSIS_DIR}/rice_fst_plot.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{ANALYSIS_DIR}/rice_fst_plot.pdf", bbox_inches='tight')
print(f"Sauvegardé : {ANALYSIS_DIR}/rice_fst_plot.png")
