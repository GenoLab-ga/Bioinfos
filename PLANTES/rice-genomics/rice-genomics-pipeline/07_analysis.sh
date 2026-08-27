#!/bin/bash
################################################################################
# RICE GENOMICS WGS PIPELINE — 07_analysis.sh
# Génétique des populations : PCA, arbre phylogénétique, Fst, π
#
# EXÉCUTION : bash 07_analysis.sh
# DURÉE     : ~30-45 min (calcul Fst sur 3.3M variants)
#
# ANALYSES :
#   1. Extraction matrice génotypique
#   2. PCA (plink2)
#   3. Arbre phylogénétique UPGMA (scipy)
#   4. Fst entre groupes variétaux (Python)
#   5. Diversité nucléotidique π par chromosome (vcftools)
#   6. Manhattan plot + densité de variants
#
# SORTIES :
#   results/06_analysis/rice_pca_plot.png
#   results/06_analysis/rice_phylo_tree.png
#   results/06_analysis/rice_fst_plot.png
#   results/06_analysis/rice_nucleotide_diversity.png
#   results/06_analysis/rice_manhattan_plot.png
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

log_message "INFO" "=== ÉTAPE 7 : ANALYSES DE GÉNÉTIQUE DES POPULATIONS ==="
eval "$CONDA_ACTIVATE"

VCF_MAF="$VCF_DIR/all_samples.maf${MAF_THRESHOLD}.vcf.gz"
[ -f "$VCF_MAF" ] && [ -s "$VCF_MAF" ] || error_exit "VCF introuvable — lancer 05_variant_calling.sh"

mkdir -p "$ANALYSIS_DIR"
cd "$ANALYSIS_DIR"

# ==============================================================================
# ÉTAPE 7.1 — Matrice génotypique
# ==============================================================================
log_message "INFO" "7.1 Extraction de la matrice génotypique..."

if [ ! -s "genotypes_raw.txt" ]; then
    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' \
        "$VCF_MAF" \
        2>/dev/null > "genotypes_raw.txt"
    log_message "INFO" "✓ Génotypes : $(wc -l < genotypes_raw.txt) variants"
else
    log_message "INFO" "✓ Matrice déjà présente"
fi

# ==============================================================================
# ÉTAPE 7.2 — PCA avec plink2
# ==============================================================================
log_message "INFO" "7.2 PCA..."

if [ ! -f "rice_pca.eigenval" ]; then
    # Conversion VCF → PLINK
    plink2 --vcf "$VCF_MAF" --make-pgen --out rice_plink \
        --threads "$MAX_THREADS" --bad-ld 2>> "$LOG_FILE"

    # Assigner des IDs aux variants
    plink2 --pfile rice_plink \
        --max-alleles 2 \
        --set-all-var-ids '@:#:$r:$a' \
        --new-id-max-allele-len 50 missing \
        --make-pgen --out rice_plink_ids \
        --threads "$MAX_THREADS" --bad-ld 2>> "$LOG_FILE"

    # Extraire les IDs valides
    grep -v "^#" rice_plink_ids.pvar | awk '$3!="." {print $3}' > valid_ids.txt

    # Fréquences alléliques
    plink2 --pfile rice_plink_ids \
        --extract valid_ids.txt \
        --freq --out rice_freq \
        --threads "$MAX_THREADS" --bad-ld 2>> "$LOG_FILE"

    # PCA
    plink2 --pfile rice_plink_ids \
        --extract valid_ids.txt \
        --read-freq rice_freq.afreq \
        --pca 7 --out rice_pca \
        --threads "$MAX_THREADS" --bad-ld 2>> "$LOG_FILE"

    log_message "INFO" "✓ PCA calculée"
else
    log_message "INFO" "✓ PCA déjà présente"
fi

# Visualisation PCA
python3 "$SCRIPT_DIR/scripts/plot_pca.py" "$ANALYSIS_DIR"
log_message "INFO" "✓ rice_pca_plot.png"

# ==============================================================================
# ÉTAPE 7.3 — Matrice de distances + arbre phylogénétique
# ==============================================================================
log_message "INFO" "7.3 Matrice de distances IBS + arbre phylogénétique..."

if [ ! -f "distance_matrix.tsv" ]; then
    python3 "$SCRIPT_DIR/scripts/calc_distances.py" "$ANALYSIS_DIR"
    log_message "INFO" "✓ Matrice de distances calculée"
else
    log_message "INFO" "✓ Matrice déjà présente"
fi

python3 "$SCRIPT_DIR/scripts/plot_tree.py" "$ANALYSIS_DIR"
log_message "INFO" "✓ rice_phylo_tree.png"

# ==============================================================================
# ÉTAPE 7.4 — Fst entre groupes
# ==============================================================================
log_message "INFO" "7.4 Calcul du Fst (~10 min)..."

if [ ! -f "fst_results.tsv" ]; then
    python3 "$SCRIPT_DIR/scripts/calc_fst.py" "$ANALYSIS_DIR"
    log_message "INFO" "✓ Fst calculé"
else
    log_message "INFO" "✓ Fst déjà calculé"
fi

python3 "$SCRIPT_DIR/scripts/plot_fst.py" "$ANALYSIS_DIR"
log_message "INFO" "✓ rice_fst_plot.png"

# ==============================================================================
# ÉTAPE 7.5 — Diversité nucléotidique π
# ==============================================================================
log_message "INFO" "7.5 Diversité nucléotidique π..."

if [ ! -f "rice_diversity.windowed.pi" ]; then
    vcftools \
        --gzvcf "$VCF_MAF" \
        --window-pi 100000 \
        --window-pi-step 50000 \
        --out rice_diversity \
        2>> "$LOG_FILE"
    log_message "INFO" "✓ π calculé"
else
    log_message "INFO" "✓ π déjà calculé"
fi

python3 "$SCRIPT_DIR/scripts/plot_diversity.py" "$ANALYSIS_DIR"
log_message "INFO" "✓ rice_nucleotide_diversity.png"

# ==============================================================================
# ÉTAPE 7.6 — Manhattan plot
# ==============================================================================
log_message "INFO" "7.6 Manhattan plot..."
python3 "$SCRIPT_DIR/scripts/plot_manhattan.py" "$ANALYSIS_DIR"
log_message "INFO" "✓ rice_manhattan_plot.png"

log_message "INFO" "=== ANALYSES TERMINÉES ==="
log_message "INFO" "Figures dans : $ANALYSIS_DIR"
ls -lh "$ANALYSIS_DIR"/*.png 2>/dev/null
