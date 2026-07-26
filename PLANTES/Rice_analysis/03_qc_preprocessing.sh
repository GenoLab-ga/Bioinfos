#!/bin/bash

################################################################################
# RICE GENOMICS WGS PIPELINE - 03_qc_preprocessing.sh
# Contrôle qualité et nettoyage des reads avec fastp
#
# POURQUOI FASTP PLUTÔT QUE FASTQC + TRIMMOMATIC ?
# - fastp fait les deux en une seule commande (QC + trimming)
# - Beaucoup plus rapide (multi-threads natif, écrit en C++)
# - Génère un rapport HTML + JSON directement (pas besoin de MultiQC séparé)
# - Détecte automatiquement les adaptateurs Illumina
# - Paramètres adaptés au WGS : on est moins agressif que pour du RNA-seq
#   car on veut conserver un maximum de couverture génomique
#
# SORTIE :
# results/01_qc/
# └── SAMEA2567493_fastp.html   ← rapport qualité visuel
# └── SAMEA2567493_fastp.json   ← rapport machine-readable
# results/02_trimmed/
# └── SAMEA2567493_R1.fastq.gz  ← reads nettoyés R1
# └── SAMEA2567493_R2.fastq.gz  ← reads nettoyés R2
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

log_message "INFO" "=== Démarrage QC et Preprocessing ==="
eval "$CONDA_ACTIVATE"


# ==============================================================================
# STEP 1 : Vérifier que fastp est installé
# Si absent, on l'installe automatiquement dans l'environnement actif
# ==============================================================================

if ! check_tool "fastp"; then
    log_message "WARN" "fastp introuvable — installation en cours..."
    mamba install -y -c bioconda fastp 2>> "$LOG_FILE" \
        || error_exit "Impossible d'installer fastp"
fi

log_message "INFO" "fastp version : $(fastp --version 2>&1 | head -1)"

# Vérifier aussi MultiQC pour agréger les rapports à la fin
if ! check_tool "multiqc"; then
    log_message "WARN" "multiqc introuvable — installation en cours..."
    mamba install -y -c bioconda multiqc 2>> "$LOG_FILE" || \
        log_message "WARN" "MultiQC non installé, rapport agrégé non disponible"
fi


# ==============================================================================
# STEP 2 : Lire le fichier sample_info.tsv
# Créé par 02_download_data.sh, il contient les chemins vers les FASTQ
# ==============================================================================

SAMPLE_INFO="$DATA_DIR/sample_info.tsv"

if [ ! -f "$SAMPLE_INFO" ]; then
    error_exit "sample_info.tsv introuvable : $SAMPLE_INFO
    Lance d'abord : bash 02_download_data.sh"
fi

# Compter les samples (on ignore la ligne d'en-tête avec tail -n +2)
SAMPLE_COUNT=$(tail -n +2 "$SAMPLE_INFO" | grep -v "^$" | wc -l)
log_message "INFO" "$SAMPLE_COUNT samples à traiter"

mkdir -p "$QC_DIR" "$TRIM_DIR"


# ==============================================================================
# STEP 3 : Lancer fastp sur chaque sample
#
# Paramètres utilisés et pourquoi :
#   --detect_adapter_for_pe   : détection auto des adaptateurs paired-end
#   --cut_front / --cut_tail  : coupe les bases de mauvaise qualité aux extrémités
#   --cut_mean_quality 20     : qualité moyenne minimale sur la fenêtre
#   --length_required 50      : longueur minimale après trimming
#   --thread $FASTP_THREADS   : parallélisation
#   --html / --json           : génère les rapports qualité
#
# NOTE WGS vs RNA-seq : pour le WGS on ne fait PAS de poly-X trimming
# et on est moins strict sur la longueur minimale pour préserver la couverture.
# ==============================================================================

PROCESSED=0

while IFS=$'\t' read -r sample_acc run_acc country variety R1 R2; do

    # Ignorer l'en-tête
    [ "$sample_acc" = "sample_accession" ] && continue
    [[ -z "$sample_acc" ]] && continue

    PROCESSED=$((PROCESSED + 1))
    log_message "INFO" "[$PROCESSED/$SAMPLE_COUNT] Traitement : $sample_acc ($country - $variety)"

    # Fichiers de sortie
    R1_TRIM="$TRIM_DIR/${sample_acc}_R1.fastq.gz"
    R2_TRIM="$TRIM_DIR/${sample_acc}_R2.fastq.gz"
    HTML_REPORT="$QC_DIR/${sample_acc}_fastp.html"
    JSON_REPORT="$QC_DIR/${sample_acc}_fastp.json"

    # Ne pas retraiter si déjà fait
    if [ -f "$R1_TRIM" ] && [ -f "$R2_TRIM" ]; then
        log_message "INFO" "  ✓ Déjà trimmé, on passe"
        continue
    fi

    # Vérifier que les fichiers d'entrée existent
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        log_message "WARN" "  ✗ Fichiers FASTQ manquants pour $sample_acc"
        log_message "WARN" "    R1 attendu : $R1"
        log_message "WARN" "    R2 attendu : $R2"
        continue
    fi

    log_message "INFO" "  → R1 : $(basename $R1)"
    log_message "INFO" "  → R2 : $(basename $R2)"

    # Lancer fastp
    fastp \
        --in1 "$R1" \
        --in2 "$R2" \
        --out1 "$R1_TRIM" \
        --out2 "$R2_TRIM" \
        --html "$HTML_REPORT" \
        --json "$JSON_REPORT" \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_tail \
        --cut_mean_quality "$FASTP_MEAN_QUALITY" \
        --length_required "$FASTP_MIN_LENGTH" \
        --thread "$FASTP_THREADS" \
        --compression 4 \
        2>> "$LOG_FILE" || {
        log_message "WARN" "  ✗ fastp a échoué pour $sample_acc"
        rm -f "$R1_TRIM" "$R2_TRIM"
        continue
    }

    # Afficher un résumé rapide depuis le JSON
    if [ -f "$JSON_REPORT" ]; then
        READS_BEFORE=$(python3 -c "import json; d=json.load(open('$JSON_REPORT')); print(d['summary']['before_filtering']['total_reads'])" 2>/dev/null || echo "?")
        READS_AFTER=$(python3 -c "import json; d=json.load(open('$JSON_REPORT')); print(d['summary']['after_filtering']['total_reads'])" 2>/dev/null || echo "?")
        Q30=$(python3 -c "import json; d=json.load(open('$JSON_REPORT')); print(round(d['summary']['after_filtering']['q30_rate']*100,1))" 2>/dev/null || echo "?")
        log_message "INFO" "  ✓ Reads avant : $READS_BEFORE | après : $READS_AFTER | Q30% : $Q30"
    fi

done < "$SAMPLE_INFO"


# ==============================================================================
# STEP 4 : Créer le tableau récapitulatif des statistiques fastp
# On parse les JSON générés pour avoir un tableau global par sample
# ==============================================================================

log_message "INFO" "Création du tableau récapitulatif QC..."

STATS_FILE="$QC_DIR/fastp_summary.tsv"
echo -e "sample\tcountry\tvariety\treads_before\treads_after\tbases_before\tbases_after\tq20_rate\tq30_rate\tgc_content\tpct_retained" > "$STATS_FILE"

while IFS=$'\t' read -r sample_acc run_acc country variety R1 R2; do
    [ "$sample_acc" = "sample_accession" ] && continue
    [[ -z "$sample_acc" ]] && continue

    JSON_REPORT="$QC_DIR/${sample_acc}_fastp.json"

    if [ ! -f "$JSON_REPORT" ]; then
        log_message "WARN" "  JSON introuvable pour $sample_acc"
        continue
    fi

    python3 << EOF >> "$STATS_FILE"
import json
d = json.load(open("$JSON_REPORT"))
b = d["summary"]["before_filtering"]
a = d["summary"]["after_filtering"]
r_before = b["total_reads"]
r_after  = a["total_reads"]
b_before = b["total_bases"]
b_after  = a["total_bases"]
q20      = round(a["q20_rate"] * 100, 2)
q30      = round(a["q30_rate"] * 100, 2)
gc       = round(a["gc_content"] * 100, 2)
pct      = round(r_after / r_before * 100, 1) if r_before > 0 else 0
print(f"$sample_acc\t$country\t$variety\t{r_before}\t{r_after}\t{b_before}\t{b_after}\t{q20}\t{q30}\t{gc}\t{pct}")
EOF

done < "$SAMPLE_INFO"

log_message "INFO" "Tableau QC : $STATS_FILE"


# ==============================================================================
# STEP 5 : MultiQC — rapport agrégé de tous les samples (optionnel mais pratique)
# MultiQC lit tous les fichiers *_fastp.json et génère un seul rapport HTML
# ==============================================================================

if check_tool "multiqc"; then
    log_message "INFO" "Génération du rapport MultiQC..."
    multiqc "$QC_DIR" \
        --outdir "$QC_DIR/multiqc" \
        --filename "multiqc_report" \
        --force \
        2>> "$LOG_FILE"
    log_message "INFO" "  ✓ Rapport MultiQC : $QC_DIR/multiqc/multiqc_report.html"
else
    log_message "WARN" "MultiQC non disponible — rapports individuels dans $QC_DIR/"
fi


# ==============================================================================
# STEP 6 : Mettre à jour sample_info.tsv avec les chemins vers les FASTQ trimmés
# Les scripts suivants (mapping) liront ces nouveaux chemins
# ==============================================================================

SAMPLE_INFO_TRIM="$DATA_DIR/sample_info_trimmed.tsv"
echo -e "sample_accession\trun_accession\tcountry\tvariety_group\tR1_trim\tR2_trim" > "$SAMPLE_INFO_TRIM"

while IFS=$'\t' read -r sample_acc run_acc country variety R1 R2; do
    [ "$sample_acc" = "sample_accession" ] && continue
    [[ -z "$sample_acc" ]] && continue

    R1_TRIM="$TRIM_DIR/${sample_acc}_R1.fastq.gz"
    R2_TRIM="$TRIM_DIR/${sample_acc}_R2.fastq.gz"

    if [ -f "$R1_TRIM" ] && [ -f "$R2_TRIM" ]; then
        echo -e "${sample_acc}\t${run_acc}\t${country}\t${variety}\t${R1_TRIM}\t${R2_TRIM}" >> "$SAMPLE_INFO_TRIM"
    fi

done < "$SAMPLE_INFO"

log_message "INFO" "Fichier mis à jour : $SAMPLE_INFO_TRIM"


# ==============================================================================
# RÉSUMÉ
# ==============================================================================

TRIM_SIZE=$(du -sh "$TRIM_DIR" 2>/dev/null | cut -f1 || echo "?")

log_message "INFO" "================================="
log_message "INFO" "=== PREPROCESSING TERMINÉ     ==="
log_message "INFO" "================================="
log_message "INFO" "Samples traités  : $PROCESSED / $SAMPLE_COUNT"
log_message "INFO" "FASTQ trimmés    : $TRIM_DIR/"
log_message "INFO" "Taille totale    : $TRIM_SIZE"
log_message "INFO" "Rapports QC      : $QC_DIR/"
log_message "INFO" "Tableau résumé   : $QC_DIR/fastp_summary.tsv"
log_message "INFO" ""
log_message "INFO" "Étape suivante : bash 04_mapping.sh"
