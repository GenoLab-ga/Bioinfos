#!/usr/bin/env python3
"""Manhattan plot de diversité π — Rice Genomics Pipeline"""
import sys, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd, numpy as np

ANALYSIS_DIR = sys.argv[1] if len(sys.argv) > 1 else "."

df = pd.read_csv(f"{ANALYSIS_DIR}/rice_diversity.windowed.pi", sep='\t')
df.columns = ['CHROM','BIN_START','BIN_END','N_VARIANTS','N_MONOMORPHIC','PI']
df = df[df['CHROM'].astype(str).isin([str(i) for i in range(1,13)])].copy()
df['CHROM'] = df['CHROM'].astype(int)
df = df.sort_values(['CHROM','BIN_START']).reset_index(drop=True)

chrom_sizes = df.groupby('CHROM')['BIN_END'].max()
offsets = {}
offset = 0
for c in range(1,13):
    offsets[c] = offset
    offset += chrom_sizes.get(c, 0) + 2000000
df['POS_CUM'] = df.apply(lambda r: r['BIN_START'] + offsets[r['CHROM']], axis=1)
centers = {c: df[df['CHROM']==c]['POS_CUM'].mean() for c in range(1,13)}

colors_alt = ['#2166AC','#4DAC26']
df['COLOR'] = df['CHROM'].apply(lambda c: colors_alt[(c-1)%2])

pi_mean  = df['PI'].mean()
pi_top5  = df['PI'].quantile(0.95)
pi_top1  = df['PI'].quantile(0.99)

fig, axes = plt.subplots(2, 1, figsize=(18, 12),
                          gridspec_kw={'height_ratios': [3,1]})
ax1, ax2 = axes

for c in range(1,13):
    data = df[df['CHROM']==c]
    ax1.scatter(data['POS_CUM']/1e6, data['PI'],
                c=data['COLOR'], s=4, alpha=0.6, linewidths=0)

ax1.axhline(pi_mean, color='gray', linewidth=1, linestyle='--',
            label=f'Moyenne (π={pi_mean:.4f})')
ax1.axhline(pi_top5, color='#FF7F00', linewidth=1.2, linestyle='--',
            label=f'95e pct (π={pi_top5:.4f})')
ax1.axhline(pi_top1, color='#E41A1C', linewidth=1.5, linestyle='--',
            label=f'99e pct (π={pi_top1:.4f})')
hotspots = df[df['PI'] >= pi_top1]
ax1.scatter(hotspots['POS_CUM']/1e6, hotspots['PI'],
            c='#E41A1C', s=12, alpha=0.9, zorder=5)

ax1.set_xticks([centers[c]/1e6 for c in range(1,13)])
ax1.set_xticklabels([f'Chr{c}' for c in range(1,13)], fontsize=10)
ax1.set_ylabel("Diversité nucléotidique (π)", fontsize=12)
ax1.set_title("Manhattan plot — Diversité nucléotidique (π)\n"
              "8 accessions · 3.3M SNPs · fenêtres 100 kb", fontsize=13)
ax1.legend(fontsize=9, loc='upper right')
ax1.set_xlim(0, df['POS_CUM'].max()/1e6)
ax1.grid(axis='y', alpha=0.2)

for c in range(1,13):
    x_end = df[df['CHROM']==c]['POS_CUM'].max()/1e6
    ax1.axvline(x_end, color='lightgray', linewidth=0.5)
    ax2.axvline(x_end, color='lightgray', linewidth=0.5)

for c in range(1,13):
    data = df[df['CHROM']==c]
    ax2.bar(data['POS_CUM']/1e6, data['N_VARIANTS'],
            width=0.04, color=colors_alt[(c-1)%2], alpha=0.7, linewidth=0)

ax2.set_xticks([centers[c]/1e6 for c in range(1,13)])
ax2.set_xticklabels([f'Chr{c}' for c in range(1,13)], fontsize=10)
ax2.set_ylabel("Variants / 100kb", fontsize=11)
ax2.set_xlabel("Position génomique", fontsize=11)
ax2.set_xlim(0, df['POS_CUM'].max()/1e6)
ax2.grid(axis='y', alpha=0.2)
ax2.set_title("Densité de variants par fenêtre", fontsize=11)

plt.tight_layout()
plt.savefig(f"{ANALYSIS_DIR}/rice_manhattan_plot.png", dpi=300, bbox_inches='tight')
plt.savefig(f"{ANALYSIS_DIR}/rice_manhattan_plot.pdf", bbox_inches='tight')
print(f"Sauvegardé : {ANALYSIS_DIR}/rice_manhattan_plot.png")
print(f"Hotspots (top 1%) : {len(hotspots)} fenêtres")
