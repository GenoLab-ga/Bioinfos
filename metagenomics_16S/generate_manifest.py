#!/usr/bin/env python3
# ============================================================
# generate_manifest.py — Génère le manifest QIIME2
# Format : PairedEndFastqManifestPhred33V2
# ============================================================

import os
import pandas as pd

FASTQ_DIR = os.path.expanduser(
    "~/Documents/github_projet/metagenomics_16S/fastq"
)
META_DIR = os.path.expanduser(
    "~/Documents/github_projet/metagenomics_16S/metadata"
)

# Lire le SraRunTable pour avoir sample-id propre
sra = pd.read_csv(f"{META_DIR}/SraRunTable.csv")

rows = []
for _, row in sra.iterrows():
    srr = row['Run']
    sample_id = row['Sample Name']  # ex: MS_687
    r1 = os.path.join(FASTQ_DIR, f"{srr}_1.fastq.gz")
    r2 = os.path.join(FASTQ_DIR, f"{srr}_2.fastq.gz")

    # Ne garder que les fichiers présents sur disque
    if os.path.exists(r1) and os.path.exists(r2):
        rows.append({
            'sample-id': sample_id,
            'forward-absolute-filepath': r1,
            'reverse-absolute-filepath': r2
        })

manifest = pd.DataFrame(rows)
out_path = f"{META_DIR}/manifest.tsv"
manifest.to_csv(out_path, sep='\t', index=False)

print(f"Manifest généré : {len(manifest)} samples → {out_path}")
print(manifest.head(3).to_string())
