import csv
from collections import Counter

with open('data/results/tbprofiler_v6/summary/all_samples.csv') as f:
    reader = csv.DictReader(f)
    rows = list(reader)

# Colonnes des drogues
drug_cols = [k for k in rows[0].keys()
             if k not in ['sample','main_lineage','sub_lineage',
                          'spoligotype','drtype','target_median_depth',
                          'pct_reads_mapped','num_reads_mapped',
                          'num_dr_variants','num_other_variants']]

print(f"Drogues testées : {len(drug_cols)}")
print("\n=== RÉSISTANCE PAR DROGUE ===")
print(f"  {'Drogue':<30} {'Résistants':>10} {'%':>8}")
print("  " + "-" * 52)

drug_resist = {}
for drug in drug_cols:
    resistant = sum(1 for r in rows
                   if r.get(drug,'') not in ['', 'Sensitive', 'S', '-'])
    drug_resist[drug] = resistant

for drug, count in sorted(drug_resist.items(),
                           key=lambda x: x[1], reverse=True):
    if count > 0:
        pct = count/len(rows)*100
        print(f"  {drug:<30} {count:>10} ({pct:.1f}%)")
