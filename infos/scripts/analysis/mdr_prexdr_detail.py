import csv

with open('data/results/tbprofiler_v6/summary/all_samples.csv') as f:
    reader = csv.DictReader(f)
    rows = list(reader)

targets = ['MDR-TB', 'Pre-XDR-TB', 'XDR-TB']
print("=== SAMPLES MDR / Pre-XDR / XDR ===")
print(f"  {'Sample':<15} {'Classification':<15} {'Depth':>8} {'Reads%':>8} {'DR variants':>12}")
print("  " + "-" * 65)
for r in rows:
    if r.get('drtype') in targets:
        print(f"  {r['sample']:<15} {r['drtype']:<15} "
              f"{float(r.get('target_median_depth',0)):>8.1f} "
              f"{float(r.get('pct_reads_mapped',0)):>8.1f}% "
              f"{r.get('num_dr_variants','?'):>12}")
