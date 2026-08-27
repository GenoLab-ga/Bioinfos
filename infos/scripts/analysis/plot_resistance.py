import csv
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

with open('data/results/tbprofiler_v6/summary/all_samples.csv') as f:
    reader = csv.DictReader(f)
    rows = list(reader)

fig, axes = plt.subplots(1, 2, figsize=(16, 7))
fig.suptitle('PRJNA1177198 — MTB Drug Resistance Profile\n'
             'Nanopore Targeted Sequencing (n=277)', 
             fontsize=14, fontweight='bold')

# ── Figure 1 : Distribution des profils de résistance ──
from collections import Counter
drtype_order = ['Susceptible','HR-TB','RR-TB','MDR-TB','Pre-XDR-TB','XDR-TB','Other']
colors = ['#2ecc71','#f39c12','#e67e22','#e74c3c','#c0392b','#8e44ad','#95a5a6']
drtype_counts = Counter(r.get('drtype','') for r in rows)
counts = [drtype_counts.get(d, 0) for d in drtype_order]
bars = axes[0].bar(drtype_order, counts, color=colors, edgecolor='white', linewidth=1.5)
axes[0].set_title('Classification de résistance', fontweight='bold', pad=15)
axes[0].set_ylabel('Nombre de samples')
axes[0].set_xlabel('')
axes[0].tick_params(axis='x', rotation=35)
for bar, count in zip(bars, counts):
    if count > 0:
        pct = count/len(rows)*100
        axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    f'{count}\n({pct:.1f}%)', ha='center', va='bottom', fontsize=9)
axes[0].set_ylim(0, max(counts)*1.2)
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

# ── Figure 2 : Résistance par drogue ──
drug_exclude = ['sample','main_lineage','sub_lineage','spoligotype','drtype',
                'target_median_depth','pct_reads_mapped','num_reads_mapped',
                'num_dr_variants','num_other_variants']
drug_cols = [k for k in rows[0].keys() if k not in drug_exclude]
drug_resist = {}
for drug in drug_cols:
    resistant = sum(1 for r in rows
                   if r.get(drug,'') not in ['','Sensitive','S','-'])
    if resistant > 0:
        drug_resist[drug] = resistant

drugs_sorted = sorted(drug_resist.items(), key=lambda x: x[1], reverse=True)
drug_names = [d[0] for d in drugs_sorted]
drug_counts = [d[1] for d in drugs_sorted]

# Colorer selon la ligne thérapeutique
def drug_color(name):
    fld = ['rifampicin','rifapentine','isoniazid','ethambutol','pyrazinamide','streptomycin']
    fq  = ['moxifloxacin','levofloxacin']
    last = ['bedaquiline','clofazimine','linezolid']
    if name in fld: return '#e74c3c'
    if name in fq:  return '#e67e22'
    if name in last: return '#8e44ad'
    return '#3498db'

bar_colors = [drug_color(d) for d in drug_names]
bars2 = axes[1].barh(drug_names, drug_counts, color=bar_colors, edgecolor='white')
axes[1].set_title('Résistance par drogue', fontweight='bold', pad=15)
axes[1].set_xlabel('Nombre de samples résistants')
axes[1].invert_yaxis()
for bar, count in zip(bars2, drug_counts):
    pct = count/len(rows)*100
    axes[1].text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                f'{count} ({pct:.1f}%)', va='center', fontsize=9)
axes[1].set_xlim(0, max(drug_counts)*1.25)
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

legend_elements = [
    mpatches.Patch(color='#e74c3c', label='Première ligne (FLD)'),
    mpatches.Patch(color='#e67e22', label='Fluoroquinolones'),
    mpatches.Patch(color='#3498db', label='Injectables / Groupe B-C'),
    mpatches.Patch(color='#8e44ad', label='Dernier recours'),
]
axes[1].legend(handles=legend_elements, loc='lower right', fontsize=9)

plt.tight_layout()
plt.savefig('data/results/tbprofiler_v6/summary/resistance_profile.png',
            dpi=150, bbox_inches='tight')
print("Figure sauvegardée : data/results/tbprofiler_v6/summary/resistance_profile.png")
