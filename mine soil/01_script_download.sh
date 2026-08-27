#!/bin/bash
# script_download.sh
set -euo pipefail

# 1. Créer la structure de répertoires
mkdir -p data/{raw,manifest} results

# 2. Générer le fichier manifest pour QIIME 2 (format Single End)
echo "sample-id,absolute-filepath" > data/manifest/single_end_manifest.csv

# Extraire les 6 premiers échantillons de "Manganese mine" ou "Abandoned site" pour un test rapide
# On filtre le CSV pour garder uniquement les runs d'intérêt (ajustez le grep selon vos besoins)
tail -n +2 SraRunTable.csv | grep -i "manganese" | head -n 6 | while IFS=, read -r run assay avglen bases proj sample model bytes center date consent filetype provider region exp country continent geo_loc instrument source latlon libname layout selection source_org platform release create ver sampname study; do
    # Nettoyer les guillemets éventuels dans les variables
    run=$(echo "$run" | tr -d '"')
    sampname=$(echo "$sampname" | tr -d '"')
    
    echo "Téléchargement de $run..."
    # Téléchargement efficace avec fasterq-dump (6 threads pour votre CPU)
    fasterq-dump --threads 6 --split-files -O data/raw/ "$run"
    
    # Ajouter au manifest (le chemin doit être absolu ou relatif au répertoire de travail)
    realpath="data/raw/${run}.fastq"
    echo "${sampname},${realpath}" >> data/manifest/single_end_manifest.csv
done

echo "Manifest généré dans data/manifest/single_end_manifest.csv"
