#!/bin/bash

ACC_FILE="SRR_Acc_List.txt"

mkdir -p data/sra_files
mkdir -p data/fastq_files

while IFS= read -r id || [ -n "$id" ]; do
    [[ -z "$id" || "$id" =~ ^# ]] && continue

    # Vérification si déjà traité
    if [ -f "data/fastq_files/${id}_1.fastq.gz" ]; then
        echo "[$id] Déjà traité — skip"
        continue
    fi

    echo "[$id] Téléchargement..."
    prefetch "$id" --output-directory data/sra_files \
        || { echo "ERREUR prefetch $id" >> errors.log; continue; }

    echo "[$id] Extraction FASTQ..."
    fasterq-dump --split-files --threads 6 \
        --outdir data/fastq_files \
        "data/sra_files/$id/$id.sra" \
        || { echo "ERREUR fasterq-dump $id" >> errors.log; continue; }

    echo "[$id] Compression..."
    pigz -p 6 data/fastq_files/${id}*.fastq
    
    
    # Ajouter dans la boucle après pigz
    rm -rf "data/sra_files/$id/"
    echo "[$id] .sra supprimé"

    echo "[$id] ✓ Terminé"

done < "$ACC_FILE"

echo "=== Pipeline terminé. Erreurs : errors.log ==="

