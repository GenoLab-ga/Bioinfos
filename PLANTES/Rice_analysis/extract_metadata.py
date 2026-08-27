import urllib.request
import xml.etree.ElementTree as ET
import time
import csv

with open('sample_accessions.txt') as f:
    samples = [l.strip() for l in f if l.strip()]

print(f"Total samples: {len(samples)}")

results = []
for i, acc in enumerate(samples[:10]):  # test sur 10 d'abord
    url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{acc}"
    try:
        with urllib.request.urlopen(url, timeout=10) as r:
            xml_data = r.read()
        root = ET.fromstring(xml_data)
        
        record = {'sample_accession': acc}
        for attr in root.iter('SAMPLE_ATTRIBUTE'):
            tag = attr.findtext('TAG', '').strip()
            val = attr.findtext('VALUE', '').strip()
            if tag in ['COUNTRY', 'Variety_Group(Tree)', 'DNA_VARNAME', 
                       'DNA_Designation', 'Source', 'Genetic_Stock_Accno']:
                record[tag] = val
        
        results.append(record)
        print(f"[{i+1}/{len(samples)}] {acc} -> {record.get('COUNTRY','?')} | {record.get('Variety_Group(Tree)','?')}")
        time.sleep(0.2)
        
    except Exception as e:
        print(f"Erreur {acc}: {e}")

# Sauvegarder
if results:
    keys = ['sample_accession','COUNTRY','Variety_Group(Tree)','DNA_VARNAME',
            'DNA_Designation','Source','Genetic_Stock_Accno']
    with open('rice_metadata.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(results)
    print(f"\nSauvegardé: rice_metadata.tsv")
