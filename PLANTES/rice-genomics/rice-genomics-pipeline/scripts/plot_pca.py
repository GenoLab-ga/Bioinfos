#!/usr/bin/env python3
"""Visualisation PCA — Rice Genomics Pipeline"""
import sys, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

ANALYSIS_DIR = sys.argv[1] if len(sys.argv) > 1 else "."

labels = {
    "SAMEA2567482": ("Brazil",       "Indica",              "#E41A1C"),
    "SAMEA2567483": ("South Korea",  "Temp. japonica",      "#377EB8"),
    "SAMEA2567484": ("Ivory Coast",  "Trop. japonica",      "#4DAF4A"),
    "SAMEA2567493": ("India",        "Basmati/sadri",       "#FF7F00"),
    "SAMEA2567500": ("Bangladesh",   "Aus/boro",            "#984EA3"),
    "SAMEA2567501": ("India",        "Aus/boro",            "#984EA3"),
    "SAMEA2567504": ("India",        "Aus/boro",            "#984EA3"),
    "SAMEA2567509": ("Bangladesh",   "Aus/boro",            "#984EA3"),
}

# Charger eigenvec et eigenval
eigenvec_file = f"{ANALYSIS_DIR}/rice_pca.eigenvec"
eigenval_file = f"{ANALYSIS_DIR}/rice_pca.eigenval"

pca_data = {}
with open(eigenvec_file) as f:
    header = f.readline().strip().split('\t')
    for line in f:
        fields = line.strip().split('\t')
        sample = fields[0]
        pca_data[sample] = [float(x) for x in fields[1:]]

eigenvals = []
with open(eigenval_file) as f:
    eigenvals = [float(l.strip()) for l in f]

total = sum(eigenvals)
pct = [round(v/total*100, 1) for v in eigenvals]

fig, ax = plt.subplots(figsize=(10, 8))

seen_groups = set()
for sample, (country, group, color) in labels.items():
    if sample not in pca_data:
        continue
    pc1, pc2 = pca_data[sample][0], pca_data[sample][1]
    label = group if group not in seen_groups else ""
    ax.scatter(pc1, pc2, c=color, s=150, zorder=5, label=label)
    ax.annotate(f"{country}\n({group})", (pc1, pc2),
                textcoords="offset points", xytext=(8, 5), fontsize=8.5)
    seen_groups.add(group)

ax.set_xlabel(f"PC1 ({pct[0]}%)", fontsize=12)
ax.set_ylabel(f"PC2 ({pct[1]}%)", fontsize=12)
ax.set_title("PCA — 3000 Rice Genomes Project\n8 accessions · 3.3M SNPs", fontsize=13)
ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax.axvline(0, color='gray', linewidth=0.5, linestyle='--')
ax.grid(True, alpha=0.3)
ax.legend(title="Variety group", loc='upper left', framealpha=0.9)

plt.tight_layout()
plt.savefig(f"{ANALYSIS_DIR}/rice_pca_plot.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{ANALYSIS_DIR}/rice_pca_plot.pdf", bbox_inches='tight')
print(f"Sauvegardé : {ANALYSIS_DIR}/rice_pca_plot.png")
