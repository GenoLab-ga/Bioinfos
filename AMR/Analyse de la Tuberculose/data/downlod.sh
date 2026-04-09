#!/bin/bash
REPERTOIRE=$(pwd)
SRR_list=SRR_Acc_List.txt

# Lecture du fichier ligne par ligne
while read -r SRR; do
    echo "Téléchargement de $SRR ..."
    fasterq-dump "$SRR" --split-files -O raw_data -p
    
    # Compression
    gzip raw_data/${SRR}_1.fastq
    gzip raw_data/${SRR}_2.fastq
done < "$SRR_list"

# Le -p c'est pour suivre la progression du téléchargement.
