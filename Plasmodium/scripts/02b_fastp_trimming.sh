#!/bin/bash

# ====================================================================
# TRIMMING: Nettoyage et QC avec fastp
# Suppression adaptateurs + filtrage qualité
# ====================================================================

# ====================================================================
# 0. VARIABLES ET CHEMINS
# ====================================================================

PROJECT_ROOT=~/Documents/github_projet/Plasmodium
input_dir="$PROJECT_ROOT/data"
trimmed_dir="$PROJECT_ROOT/data_trimmed"
fastp_reports="$PROJECT_ROOT/qc_reports/fastp"
log_dir="$PROJECT_ROOT/logs"

# Échantillons
sample_list=(
    ERR11471971
    ERR11471975
    ERR11471979
    ERR11471985
)

# Paramètres fastp
THREADS=12
MIN_LENGTH=30
QUALITY_THRESHOLD=20  # Phred score minimum

mkdir -p "$trimmed_dir" "$fastp_reports" "$log_dir"

echo "======================================================"
echo "TRIMMING: Fastp - Nettoyage des reads"
echo "======================================================"
echo ""

# ====================================================================
# 1. VÉRIFICATIONS
# ====================================================================

echo "======================================================"
echo "ÉTAPE 1: Vérifications préalables"
echo "======================================================"
echo ""

# Vérifier fastp
if ! command -v fastp &> /dev/null; then
    echo "❌ Erreur: fastp non trouvé"
    echo "Installation: mamba install fastp"
    exit 1
fi

echo "Version fastp: $(fastp --version)"
echo ""

# Vérifier les fichiers
echo "Vérification fichiers FASTQ brutes:"
for sample in "${sample_list[@]}"; do
    if [ -f "$input_dir/${sample}_1.fastq.gz" ] && [ -f "$input_dir/${sample}_2.fastq.gz" ]; then
        r1_size=$(du -h "$input_dir/${sample}_1.fastq.gz" | cut -f1)
        r2_size=$(du -h "$input_dir/${sample}_2.fastq.gz" | cut -f1)
        echo "  ✓ $sample: R1=$r1_size, R2=$r2_size"
    else
        echo "  ❌ $sample: fichiers manquants"
        exit 1
    fi
done

echo ""
echo "Paramètres fastp:"
echo "  Longueur minimum: $MIN_LENGTH bp"
echo "  Qualité minimum: Q$QUALITY_THRESHOLD"
echo "  Threads: $THREADS"
echo "  Auto-détection adaptateurs: ON"
echo ""

# ====================================================================
# 2. TRIMMING PAR ÉCHANTILLON
# ====================================================================

echo "======================================================"
echo "ÉTAPE 2: Trimming et filtrage avec fastp"
echo "======================================================"
echo ""

for sample in "${sample_list[@]}"; do
    R1="$input_dir/${sample}_1.fastq.gz"
    R2="$input_dir/${sample}_2.fastq.gz"
    R1_TRIMMED="$trimmed_dir/${sample}_1.trimmed.fastq.gz"
    R2_TRIMMED="$trimmed_dir/${sample}_2.trimmed.fastq.gz"
    REPORT_JSON="$fastp_reports/${sample}.fastp.json"
    REPORT_HTML="$fastp_reports/${sample}.fastp.html"
    
    # Vérifier si déjà trimmé
    if [ -f "$R1_TRIMMED" ] && [ -f "$R2_TRIMMED" ]; then
        echo "⚠ $sample: déjà trimmé (ignoré)"
        continue
    fi
    
    echo "=========================================="
    echo "Traitement: $sample"
    echo "=========================================="
    echo "  R1: $(basename $R1)"
    echo "  R2: $(basename $R2)"
    echo ""
    
    # Lancer fastp
    # --detect_adapter_for_pe: auto-détection adaptateurs (paired-end)
    # --qualified_quality_phred: Q score minimum
    # --length_required: longueur minimum après trimming
    # --json/--html: rapports de sortie
    fastp -i "$R1" -I "$R2" \
          -o "$R1_TRIMMED" -O "$R2_TRIMMED" \
          --detect_adapter_for_pe \
          --qualified_quality_phred "$QUALITY_THRESHOLD" \
          --length_required "$MIN_LENGTH" \
          --thread "$THREADS" \
          --json "$REPORT_JSON" \
          --html "$REPORT_HTML" \
          2>&1 | tee "$log_dir/fastp_${sample}.log"
    
    if [ $? -eq 0 ]; then
        echo "✓ Trimming réussi"
        
        # Afficher stats depuis le JSON
        echo ""
        echo "Statistiques $sample:"
        if [ -f "$REPORT_JSON" ]; then
            # Extraire quelques stats du JSON
            reads_before=$(grep -o '"total_reads":[0-9]*' "$REPORT_JSON" | head -1 | cut -d: -f2)
            reads_after=$(grep -o '"passed_filter_reads":[0-9]*' "$REPORT_JSON" | head -1 | cut -d: -f2)
            
            if [ ! -z "$reads_before" ] && [ ! -z "$reads_after" ]; then
                pct=$(echo "scale=1; $reads_after * 100 / $reads_before" | bc)
                echo "  Reads avant: $reads_before"
                echo "  Reads après: $reads_after ($pct%)"
            fi
        fi
    else
        echo "❌ Erreur trimming $sample"
        exit 1
    fi
    
    echo ""
done

echo ""

# ====================================================================
# 3. RÉSUMÉ
# ====================================================================

echo "======================================================"
echo "RÉSUMÉ TRIMMING"
echo "======================================================"
echo ""

echo "Fichiers trimmés générés:"
ls -lh "$trimmed_dir"/*trimmed.fastq.gz 2>/dev/null | awk '{print "  " $9, "(" $5 ")"}'

echo ""
echo "Rapports fastp (HTML):"
ls -lh "$fastp_reports"/*.html 2>/dev/null | awk '{print "  " $9}'

echo ""
echo "Réduction de taille (R1):"
for sample in "${sample_list[@]}"; do
    if [ -f "$input_dir/${sample}_1.fastq.gz" ] && [ -f "$trimmed_dir/${sample}_1.trimmed.fastq.gz" ]; then
        r1_orig=$(du -h "$input_dir/${sample}_1.fastq.gz" | cut -f1)
        r1_trim=$(du -h "$trimmed_dir/${sample}_1.trimmed.fastq.gz" | cut -f1)
        echo "  $sample: $r1_orig → $r1_trim"
    fi
done

echo ""

# ====================================================================
# 4. LOG ET FINALISATION
# ====================================================================

echo "======================================================"
echo "RAPPORT FINAL"
echo "======================================================"
echo ""

echo "✓ Trimming fastp complété"
echo ""
echo "Fichiers de sortie:"
echo "  Localisation: $trimmed_dir/"
echo "  Format: ${sample}_1.trimmed.fastq.gz et _2.trimmed.fastq.gz"
echo ""
echo "Rapports interactifs:"
echo "  Localisation: $fastp_reports/"
echo "  Format: ${sample}.fastp.html"
echo ""
echo "Consulter les rapports:"
echo "  firefox $fastp_reports/*.html &"
echo ""

# Log
echo "[$(date)] Fastp trimming completed" >> "$log_dir/analysis.log"

echo "======================================================"
echo "✓ ÉTAPE TRIMMING FASTP COMPLÉTÉE"
echo "======================================================"
echo ""
echo "PROCHAINES ÉTAPES:"
echo "1. Consulter rapports HTML pour vérifier trimming"
echo "2. Lancer quantification: bash 03_salmon_quantification.sh"
echo ""
