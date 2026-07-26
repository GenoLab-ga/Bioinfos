#!/usr/bin/env python3
"""Diversité nucléotidique π par chromosome — Rice Genomics Pipeline"""
import sys, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd, numpy as np

ANALYSIS_DIR = sys.argv[1] if len(sys.argv) > 1 else "."

df = pd.read_csv(f"{ANALYSIS_DIR}/rice_diversity.windowed.pi", sep='\t')
df.columns = ['CHROM','BIN_START','BIN_END','N_VARIANTS','N_MONOMORPHIC','PI']
df = df[df['CHROM'].astype(str).isin([str(i) for i in range(1,13)])].copy()
df['CHROM'] = df['CHROM'].astype(int)
df = df.sort_values(['CHROM','BIN_START'])

colors = ['#2166AC','#4DAC26'] * 6
fig, axes = plt.subplots(12, 1, figsize=(16, 20), sharex=False)

for idx, chrom in enumerate(range(1,13)):
    ax = axes[idx]
    data = df[df['CHROM']==chrom]
    if len(data) == 0: continue
    x = data['BIN_START'] / 1e6
    y = data['PI']
    ax.fill_between(x, y, alpha=0.6, color=colors[idx])
    ax.plot(x, y, color=colors[idx], linewidth=0.8)
    pi_mean = y.mean()
    ax.axhline(pi_mean, color='red', linewidth=0.8, linestyle='--', alpha=0.7)
    ax.set_ylabel(f'Chr{chrom}', fontsize=9, rotation=0, labelpad=35, va='center')
    ax.set_xlim(0, x.max())
    ax.set_ylim(0, max(y.max()*1.1, 0.001))
    ax.tick_params(labelsize=7)
    ax.grid(axis='x', alpha=0.2)
    ax.text(0.99, 0.85, f'π={pi_mean:.4f}', transform=ax.transAxes,
            fontsize=7.5, ha='right', color='red',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

axes[-1].set_xlabel("Position génomique (Mb)", fontsize=11)
fig.suptitle("Diversité nucléotidique (π) — 12 chromosomes du riz\n"
             "Fenêtres 100 kb · 8 accessions · 3.3M SNPs",
             fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(f"{ANALYSIS_DIR}/rice_nucleotide_diversity.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{ANALYSIS_DIR}/rice_nucleotide_diversity.pdf", bbox_inches='tight')
print(f"Sauvegardé : {ANALYSIS_DIR}/rice_nucleotide_diversity.png")
