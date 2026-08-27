#!/usr/bin/env python3
import csv
import os
import glob

# Chemins des fichiers
sra_csv = "SraRunTable.csv"
manifest_out = "data/manifest/paired_end_manifest.csv"
raw_dir = "data/raw"

# Créer le dossier s'il n'existe pas
os.makedirs("data/manifest", exist_ok=True)

# Dictionnaire pour mapper rapidement Run -> Sample Name
run_to_sample = {}
with open(sra_csv, 'r', encoding='utf-8-sig') as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        run = row['Run'].strip()
        sampname = row['Sample Name'].strip()
        run_to_sample[run] = sampname

with open(manifest_out, 'w', encoding='utf-8', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath'])

    added_count = 0
    for fwd_file in glob.glob(os.path.join(raw_dir, "*_1.fastq")):
        filename = os.path.basename(fwd_file)
        run = filename.replace("_1.fastq", "")
        rev_file = fwd_file.replace("_1.fastq", "_2.fastq")

        if os.path.exists(rev_file) and run in run_to_sample:
            sampname = run_to_sample[run]
            abs_fwd = os.path.abspath(fwd_file)
            abs_rev = os.path.abspath(rev_file)

            writer.writerow([sampname, abs_fwd, abs_rev])
            added_count += 1

print(f"\n🎉 Manifest généré avec succès ! {added_count} échantillons ajoutés.")
print(f"📁 Chemin du manifest : {os.path.abspath(manifest_out)}")
