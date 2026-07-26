#!/bin/bash

################################################################################
# RICE GENOMICS WGS PIPELINE - 02_download_data.sh
# Téléchargement des FASTQ depuis ENA (European Nucleotide Archive)
# via wget — PAS de SRA Toolkit nécessaire
#
# POURQUOI ENA ET PAS SRA ?
# Les données du projet PRJEB6180 sont hébergées sur ENA (accessions PRJEB).
# ENA fournit des liens FTP directs récupérables avec wget, ce qui est
# plus simple et plus rapide que prefetch/fasterq-dump du SRA Toolkit.
# On a déjà identifié ces liens dans PRJEB6180_full.tsv.
#
# STRUCTURE DES FICHIERS APRÈS TÉLÉCHARGEMENT :
# raw_fastq/
# └── SAMEA2567493/          ← dossier par sample_accession
#     ├── ERR614226_1.fastq.gz   ← read 1 (forward)
#     └── ERR614226_2.fastq.gz   ← read 2 (reverse)
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

log_message "INFO" "=== Démarrage du téléchargement des données ==="
print_config

# ==============================================================================
# STEP 1 : Vérifier le fichier des échantillons
# Le fichier samples.tsv doit avoir été créé lors du setup
# Format : sample_accession <TAB> run_accession <TAB> country <TAB> variety_group
# ==============================================================================

if [ ! -f "$SAMPLES_FILE" ]; then
    error_exit "Fichier samples.tsv introuvable : $SAMPLES_FILE
    Crée-le avec les 8 échantillons sélectionnés."
fi

SAMPLE_COUNT=$(grep -v "^#" "$SAMPLES_FILE" | grep -v "^$" | wc -l)
log_message "INFO" "$SAMPLE_COUNT échantillons trouvés dans $SAMPLES_FILE"


# ==============================================================================
# STEP 2 : Récupérer les liens FTP depuis le fichier ENA téléchargé
# On utilise PRJEB6180_full.tsv (déjà présent dans ton dossier d'analyse)
# pour extraire les colonnes run_accession (col 1) et fastq_ftp (col 95)
# ==============================================================================

# Chemin vers le TSV ENA complet (adapte si nécessaire)
ENA_TSV="$HOME/PRJEB6180_full.tsv"

if [ ! -f "$ENA_TSV" ]; then
    log_message "WARN" "PRJEB6180_full.tsv introuvable dans $HOME"
    log_message "WARN" "Tentative de re-téléchargement depuis ENA..."
    wget -q -O "$ENA_TSV" \
      "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB6180&result=read_run&fields=all&format=tsv&limit=0" \
      || error_exit "Impossible de télécharger le fichier ENA"
fi

log_message "INFO" "Fichier ENA TSV : $ENA_TSV"


# ==============================================================================
# STEP 3 : Téléchargement des FASTQ
#
# Pour chaque ligne du fichier samples.tsv :
#   1. On lit le sample_accession et le run_accession
#   2. On cherche les liens FTP dans PRJEB6180_full.tsv
#   3. On filtre pour ne garder que les fichiers _1.fastq.gz et _2.fastq.gz
#   4. On télécharge avec wget -c (reprend si interrompu)
#
# NOTE : certains runs ont 3 fichiers (un single + paired). On ignore le single.
# ==============================================================================

mkdir -p "$RAW_FASTQ_DIR"

DOWNLOADED=0
SKIPPED=0
FAILED=0

while IFS=$'\t' read -r sample_acc run_acc country variety; do

    # Ignorer les commentaires et lignes vides
    [[ "$sample_acc" =~ ^#.*$ ]] && continue
    [[ -z "$sample_acc" ]] && continue

    log_message "INFO" "--- Traitement : $sample_acc ($run_acc) | $country | $variety ---"

    # Créer le dossier de destination pour ce sample
    SAMPLE_DIR="$RAW_FASTQ_DIR/$sample_acc"
    mkdir -p "$SAMPLE_DIR"

    # Vérifier si déjà téléchargé (les deux fichiers R1 et R2 existent)
    R1="$SAMPLE_DIR/${run_acc}_1.fastq.gz"
    R2="$SAMPLE_DIR/${run_acc}_2.fastq.gz"

    if [ -f "$R1" ] && [ -f "$R2" ]; then
        log_message "INFO" "  ✓ Déjà présent, on passe"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Extraire les liens FTP pour ce run depuis le TSV ENA
    # La colonne 95 contient les URLs séparées par ";"
    FTP_LINKS=$(grep "^${run_acc}" "$ENA_TSV" | cut -f95)

    if [ -z "$FTP_LINKS" ]; then
        log_message "WARN" "  ✗ Aucun lien FTP trouvé pour $run_acc"
        FAILED=$((FAILED + 1))
        continue
    fi

    # Boucler sur chaque URL (séparées par ";")
    # On ne télécharge que _1.fastq.gz et _2.fastq.gz (paired-end)
    IFS=';' read -ra URLS <<< "$FTP_LINKS"

    for url in "${URLS[@]}"; do

        # Ignorer les fichiers single-end (sans _1 ou _2)
        if [[ ! "$url" =~ _1\.fastq\.gz$ ]] && [[ ! "$url" =~ _2\.fastq\.gz$ ]]; then
            log_message "INFO" "  → Ignoré (single-end) : $(basename $url)"
            continue
        fi

        DEST_FILE="$SAMPLE_DIR/$(basename $url)"

        if [ -f "$DEST_FILE" ]; then
            log_message "INFO" "  ✓ Fichier déjà présent : $(basename $url)"
            continue
        fi

        log_message "INFO" "  → Téléchargement : $(basename $url)"

        # wget -c : reprend le téléchargement si interrompu
        # https:// devant l'URL ENA (le TSV donne des URLs sans protocole)
        wget -c -q --show-progress \
            "https://${url}" \
            -O "$DEST_FILE" \
            2>> "$LOG_FILE" || {
            log_message "WARN" "  ✗ Échec du téléchargement : $(basename $url)"
            rm -f "$DEST_FILE"   # Supprimer le fichier partiel
            FAILED=$((FAILED + 1))
            continue
        }

        log_message "INFO" "  ✓ Téléchargé : $(basename $url) ($(du -h $DEST_FILE | cut -f1))"

    done

    DOWNLOADED=$((DOWNLOADED + 1))

done < "$SAMPLES_FILE"


# ==============================================================================
# STEP 4 : Vérification de l'intégrité
# On vérifie que chaque paire de fichiers est présente et non vide
# ==============================================================================

log_message "INFO" "=== Vérification des fichiers téléchargés ==="

OK=0
MISSING=0

while IFS=$'\t' read -r sample_acc run_acc country variety; do
    [[ "$sample_acc" =~ ^#.*$ ]] && continue
    [[ -z "$sample_acc" ]] && continue

    R1="$RAW_FASTQ_DIR/$sample_acc/${run_acc}_1.fastq.gz"
    R2="$RAW_FASTQ_DIR/$sample_acc/${run_acc}_2.fastq.gz"

    if [ -f "$R1" ] && [ -f "$R2" ] && [ -s "$R1" ] && [ -s "$R2" ]; then
        SIZE_R1=$(du -h "$R1" | cut -f1)
        SIZE_R2=$(du -h "$R2" | cut -f1)
        log_message "INFO" "  ✓ $sample_acc : R1=$SIZE_R1 | R2=$SIZE_R2"
        OK=$((OK + 1))
    else
        log_message "WARN" "  ✗ Fichiers manquants ou vides pour $sample_acc ($run_acc)"
        MISSING=$((MISSING + 1))
    fi

done < "$SAMPLES_FILE"


# ==============================================================================
# STEP 5 : Créer le fichier sample_info.tsv pour les étapes suivantes
# Ce fichier sera lu par 03_qc_preprocessing.sh et les scripts suivants
# Format : sample_accession <TAB> run_accession <TAB> country <TAB> variety <TAB> R1 <TAB> R2
# ==============================================================================

SAMPLE_INFO="$DATA_DIR/sample_info.tsv"

echo -e "sample_accession\trun_accession\tcountry\tvariety_group\tR1\tR2" > "$SAMPLE_INFO"

while IFS=$'\t' read -r sample_acc run_acc country variety; do
    [[ "$sample_acc" =~ ^#.*$ ]] && continue
    [[ -z "$sample_acc" ]] && continue

    R1="$RAW_FASTQ_DIR/$sample_acc/${run_acc}_1.fastq.gz"
    R2="$RAW_FASTQ_DIR/$sample_acc/${run_acc}_2.fastq.gz"

    if [ -f "$R1" ] && [ -f "$R2" ]; then
        echo -e "${sample_acc}\t${run_acc}\t${country}\t${variety}\t${R1}\t${R2}" >> "$SAMPLE_INFO"
    fi

done < "$SAMPLES_FILE"

log_message "INFO" "Fichier d'information samples créé : $SAMPLE_INFO"


# ==============================================================================
# RÉSUMÉ
# ==============================================================================

TOTAL_SIZE=$(du -sh "$RAW_FASTQ_DIR" | cut -f1)

log_message "INFO" "=============================="
log_message "INFO" "=== TÉLÉCHARGEMENT TERMINÉ ==="
log_message "INFO" "=============================="
log_message "INFO" "Téléchargés  : $DOWNLOADED"
log_message "INFO" "Déjà présents: $SKIPPED"
log_message "INFO" "Échoués      : $FAILED"
log_message "INFO" "Vérifiés OK  : $OK / $SAMPLE_COUNT"
log_message "INFO" "Taille totale: $TOTAL_SIZE"
log_message "INFO" ""
log_message "INFO" "Étape suivante : bash 03_qc_preprocessing.sh"
