#!/bin/bash

RUNLIST="data/SRR_Acc_List_isolates.txt"
OUTDIR_FASTQ="data/fastq"
TMPDIR="data/sra_tmp"

mkdir -p "$OUTDIR_FASTQ" "$TMPDIR" logs

TOTAL=$(wc -l < "$RUNLIST")
COUNT=0

while read -r RUN; do
    COUNT=$((COUNT + 1))

    # Sauter si déjà traité
    if [ -f "$OUTDIR_FASTQ/${RUN}.fastq.gz" ]; then
        echo "[$COUNT/$TOTAL] $RUN déjà traité, skip."
        continue
    fi

    echo "[$COUNT/$TOTAL] Téléchargement de $RUN..."

    prefetch "$RUN" \
        --output-directory "$TMPDIR" \
        --max-size 50G \
        --progress

    fasterq-dump "$TMPDIR/$RUN/$RUN.sra" \
        --outdir "$OUTDIR_FASTQ" \
        --threads 8 \
        --progress

    pigz -p 8 "$OUTDIR_FASTQ/${RUN}.fastq"

    rm -rf "$TMPDIR/$RUN/"

    echo "[$COUNT/$TOTAL] $RUN OK → ${RUN}.fastq.gz"

done < "$RUNLIST"

echo "=== Terminé : $TOTAL samples ==="
