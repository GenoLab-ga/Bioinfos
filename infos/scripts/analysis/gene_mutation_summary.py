import json
import os
from collections import Counter

results_dir = 'data/results/tbprofiler_v6/results'
gene_counter = Counter()
mutation_counter = Counter()

for fname in os.listdir(results_dir):
    if not fname.endswith('.results.json'):
        continue
    with open(f'{results_dir}/{fname}') as f:
        data = json.load(f)
    for v in data.get('dr_variants', []):
        gene = v.get('gene_name') or 'N/A'
        change = v.get('change') or 'N/A'
        gene_counter[gene] += 1
        mutation_counter[f"{gene}:{change}"] += 1

print("=== GÈNES PORTEURS DE MUTATIONS DE RÉSISTANCE ===")
print(f"  {'Gène':<15} {'Samples':>10} {'%':>8}")
print("  " + "-" * 38)
for gene, count in gene_counter.most_common():
    print(f"  {gene:<15} {count:>10} ({count/277*100:.1f}%)")

print("\n=== TOP 20 MUTATIONS LES PLUS FRÉQUENTES ===")
print(f"  {'Mutation':<35} {'Samples':>10} {'%':>8}")
print("  " + "-" * 58)
for mut, count in mutation_counter.most_common(20):
    print(f"  {mut:<35} {count:>10} ({count/277*100:.1f}%)")
