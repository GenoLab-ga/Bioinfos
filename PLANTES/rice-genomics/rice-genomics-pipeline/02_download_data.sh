#!/bin/bash
################################################################################
# RICE GENOMICS WGS PIPELINE — 02_download_data.sh
# Téléchargement des FASTQ depuis ENA via wget (liens FTP directs)
#
# EXÉCUTION : bash 02_download_data.sh
# DURÉE     : ~30-60 min selon la connexion (5.7 Go total)
# REPRISE   : wget -c permet de reprendre un téléchargement interrompu
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

log_message "INFO" "=== ÉTAPE 2 : TÉLÉCHARGEMENT DES DONNÉES ==="
print_config

# --- Vérifier samples.tsv ---
[ -f "$SAMPLES_FILE" ] || error_exit "samples.tsv introuvable : $SAMPLES_FILE"
SAMPLE_COUNT=$(grep -v "^#" "$SAMPLES_FILE" | grep -v "^$" | grep -v "^sample_accession" | wc -l)
log_message "INFO" "$SAMPLE_COUNT échantillons à télécharger"

# --- Récupérer le fichier ENA si absent ---
if [ ! -f "$ENA_TSV" ]; then
    log_message "INFO" "Téléchargement du fichier métadonnées ENA..."
    wget -q -O "$ENA_TSV" \
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB6180&result=read_run&fields=all&format=tsv&limit=0" \
        || error_exit "Impossible de télécharger le fichier ENA"
    log_message "INFO" "✓ ENA TSV : $ENA_TSV"
fi

mkdir -p "$RAW_FASTQ_DIR"

DOWNLOADED=0; SKIPPED=0; FAILED=0

while IFS=$'\t' read -r sample_acc run_acc country variety; do
    [[ "$sample_acc" =~ ^#.*$ ]] && continue
    [ "$sample_acc" = "sample_accession" ] && continue
    [[ -z "$sample_acc" ]] && continue

    log_message "INFO" "--- $sample_acc | $country | $variety ---"

    SAMPLE_DIR="$RAW_FASTQ_DIR/$sample_acc"
    mkdir -p "$SAMPLE_DIR"

    R1="$SAMPLE_DIR/${run_acc}_1.fastq.gz"
    R2="$SAMPLE_DIR/${run_acc}_2.fastq.gz"

    if [ -f "$R1" ] && [ -f "$R2" ] && [ -s "$R1" ] && [ -s "$R2" ]; then
        log_message "INFO" "  ✓ Déjà présent"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Extraire les URLs FTP depuis le TSV ENA (colonne fastq_ftp)
    FTP_LINKS=$(grep "^${run_acc}" "$ENA_TSV" | cut -f95)
    [ -z "$FTP_LINKS" ] && { log_message "WARN" "  ✗ Aucun lien FTP pour $run_acc"; FAILED=$((FAILED+1)); continue; }

    IFS=';' read -ra URLS <<< "$FTP_LINKS"
    for url in "${URLS[@]}"; do
        [[ ! "$url" =~ _1\.fastq\.gz$ ]] && [[ ! "$url" =~ _2\.fastq\.gz$ ]] && continue
        DEST="$SAMPLE_DIR/$(basename "$url")"
        [ -f "$DEST" ] && [ -s "$DEST" ] && continue
        log_message "INFO" "  → $(basename $url)"
        wget -c -q --show-progress "https://${url}" -O "$DEST" 2>> "$LOG_FILE" || {
            log_message "WARN" "  ✗ Échec : $(basename $url)"
            rm -f "$DEST"
            FAILED=$((FAILED+1))
            continue
        }
        log_message "INFO" "  ✓ $(basename $url) — $(du -h $DEST | cut -f1)"
    done
    DOWNLOADED=$((DOWNLOADED+1))

done < "$SAMPLES_FILE"

# --- Créer sample_info.tsv ---
echo -e "sample_accession\trun_accession\tcountry\tvariety_group\tR1\tR2" > "$SAMPLE_INFO"

while IFS=$'\t' read -r sample_acc run_acc country variety; do
    [[ "$sample_acc" =~ ^#.*$ ]] && continue
    [ "$sample_acc" = "sample_accession" ] && continue
    [[ -z "$sample_acc" ]] && continue
    R1="$RAW_FASTQ_DIR/$sample_acc/${run_acc}_1.fastq.gz"
    R2="$RAW_FASTQ_DIR/$sample_acc/${run_acc}_2.fastq.gz"
    [ -f "$R1" ] && [ -f "$R2" ] && \
        echo -e "${sample_acc}\t${run_acc}\t${country}\t${variety}\t${R1}\t${R2}" >> "$SAMPLE_INFO"
done < "$SAMPLES_FILE"

TOTAL_SIZE=$(du -sh "$RAW_FASTQ_DIR" | cut -f1)
log_message "INFO" "=== TÉLÉCHARGEMENT TERMINÉ ==="
log_message "INFO" "Téléchargés : $DOWNLOADED | Existants : $SKIPPED | Échecs : $FAILED"
log_message "INFO" "Taille totale : $TOTAL_SIZE"
log_message "INFO" "Étape suivante : bash 03_qc_preprocessing.sh"
