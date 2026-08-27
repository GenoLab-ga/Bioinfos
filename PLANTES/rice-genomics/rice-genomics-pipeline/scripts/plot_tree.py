#!/usr/bin/env python3
"""Arbre phylogénétique UPGMA — Rice Genomics Pipeline"""
import sys, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

ANALYSIS_DIR = sys.argv[1] if len(sys.argv) > 1 else "."

labels = [
    "Brazil (Indica)", "S.Korea (Temp.japonica)", "IvoryCoast (Trop.japonica)",
    "India (Basmati)", "Bangladesh (Aus/boro)", "India Aus/boro-1",
    "India Aus/boro-2", "Bangladesh Aus/boro-2"
]
color_map = {
    "Brazil (Indica)":            "#E41A1C",
    "S.Korea (Temp.japonica)":    "#377EB8",
    "IvoryCoast (Trop.japonica)": "#4DAF4A",
    "India (Basmati)":            "#FF7F00",
    "Bangladesh (Aus/boro)":      "#984EA3",
    "India Aus/boro-1":           "#984EA3",
    "India Aus/boro-2":           "#984EA3",
    "Bangladesh Aus/boro-2":      "#984EA3",
}

# Charger la matrice de distances
dist = np.zeros((8,8))
with open(f"{ANALYSIS_DIR}/distance_matrix.tsv") as f:
    f.readline()
    for i, line in enumerate(f):
        vals = line.strip().split('\t')[1:]
        for j, v in enumerate(vals):
            dist[i,j] = float(v)

Z = linkage(squareform(dist), method='average')

fig, ax = plt.subplots(figsize=(11, 7))
dend = dendrogram(Z, labels=labels, orientation='right', ax=ax,
                  color_threshold=0, above_threshold_color='#aaaaaa',
                  leaf_font_size=12, leaf_rotation=0)

for lbl in ax.get_yticklabels():
    txt = lbl.get_text()
    if txt in color_map:
        lbl.set_color(color_map[txt])
        lbl.set_fontweight('bold')

ax.set_xlabel("Distance IBS", fontsize=12)
ax.set_title("Arbre phylogénétique UPGMA — 3000 Rice Genomes Project\n"
             "8 accessions · 3.3M SNPs · distances IBS", fontsize=13)

legend = [
    Patch(color="#E41A1C", label="Indica"),
    Patch(color="#377EB8", label="Temp. japonica"),
    Patch(color="#4DAF4A", label="Trop. japonica"),
    Patch(color="#FF7F00", label="Basmati/sadri"),
    Patch(color="#984EA3", label="Aus/boro"),
]
ax.legend(handles=legend, loc='lower right', fontsize=11, framealpha=0.9)
ax.grid(axis='x', alpha=0.3)
plt.tight_layout()
plt.savefig(f"{ANALYSIS_DIR}/rice_phylo_tree.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{ANALYSIS_DIR}/rice_phylo_tree.pdf", bbox_inches='tight')
print(f"Sauvegardé : {ANALYSIS_DIR}/rice_phylo_tree.png")
