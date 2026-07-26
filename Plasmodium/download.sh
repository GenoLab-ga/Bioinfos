#!/bin/bash

# Le nom de ton fichier contenant les ERRxxxxxx
ACC_FILE="SRR_Acc_List.txt" 

mkdir -p data/sra_files
mkdir -p data/fastq_files

while IFS= read -r id || [ -n "$id" ]; do
    [[ -z "$id" || "$id" =~ ^# ]] && continue

    echo "-------------------------------------------------------"
    echo "Traitement de l'identifiant : $id"
    
    # Prefetch fonctionne avec ERR, SRR ou DRR sans distinction
    prefetch "$id" --output-directory data/sra_files

    # On extrait les FASTQ
    # Attention : fasterq-dump cherche le dossier créé par prefetch
    echo "Extraction des fichiers FASTQ pour $id..."
    fasterq-dump --split-files --outdir data/fastq_files "data/sra_files/$id/$id.sra"

    echo "Terminé pour $id"

done < "$ACC_FILE"

echo "_________________________________________________________________________"
echo "Tous les téléchargements et extractions sont terminés."
