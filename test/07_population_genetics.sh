#!/bin/bash
# =============================================================
# Étape 7 : Génétique des populations
# Projet  : PRJNA1465284 – P. falciparum drug resistance
# Auteur  : GenoLabGab
# =============================================================

set -euo pipefail

# ── Paramètres ───────────────────────────────────────────────
VCF="results/vcf/pass_variants.vcf.gz"
VCF_RENAMED="results/annotation/renamed_variants.vcf.gz"
POP_DIR="results/population_genetics"
LOG="logs/07_popgen.log"

mkdir -p "${POP_DIR}/pca" \
         "${POP_DIR}/frequencies" \
         "${POP_DIR}/relatedness" \
         "${POP_DIR}/haplotypes" \
         logs

echo "=======================================" | tee "${LOG}"
echo " Pop. genetics démarrée : $(date)"       | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"

# ── Vérification outils ───────────────────────────────────────
for tool in vcftools plink Rscript; do
    if ! command -v "${tool}" &>/dev/null; then
        echo "ERREUR : ${tool} introuvable." | tee -a "${LOG}"
        exit 1
    fi
done

# =============================================================
# MODULE 1 — Fréquences alléliques par locus
# =============================================================
echo "" | tee -a "${LOG}"
echo "--- MODULE 1 : Fréquences alléliques ---" | tee -a "${LOG}"

vcftools \
    --gzvcf "${VCF}" \
    --freq \
    --out "${POP_DIR}/frequencies/all_loci" \
    2>>"${LOG}"

vcftools \
    --gzvcf "${VCF}" \
    --counts \
    --out "${POP_DIR}/frequencies/all_loci" \
    2>>"${LOG}"

echo "Fréquences calculées." | tee -a "${LOG}"

# Fréquences par locus cible
for REGION in \
    "Pf3D7_04_v3:747928-749845:dhfr" \
    "Pf3D7_05_v3:957881-962149:mdr1" \
    "Pf3D7_07_v3:403223-404903:crt" \
    "Pf3D7_08_v3:549994-555749:dhps" \
    "Pf3D7_13_v3:1725260-1726923:kelch13"; do

    CHR=$(echo "${REGION}" | cut -d: -f1)
    COORDS=$(echo "${REGION}" | cut -d: -f2)
    GENE=$(echo "${REGION}" | cut -d: -f3)
    START=$(echo "${COORDS}" | cut -d- -f1)
    END=$(echo "${COORDS}" | cut -d- -f2)

    vcftools \
        --gzvcf "${VCF}" \
        --chr "${CHR}" \
        --from-bp "${START}" \
        --to-bp "${END}" \
        --freq \
        --out "${POP_DIR}/frequencies/${GENE}" \
        2>>"${LOG}"
done

echo "Fréquences par gène calculées." | tee -a "${LOG}"

# =============================================================
# MODULE 2 — Statistiques de diversité
# =============================================================
echo "" | tee -a "${LOG}"
echo "--- MODULE 2 : Diversité génétique ---" | tee -a "${LOG}"

# Pi (diversité nucléotidique) par fenêtre
vcftools \
    --gzvcf "${VCF}" \
    --window-pi 500 \
    --out "${POP_DIR}/frequencies/nucleotide_diversity" \
    2>>"${LOG}"

# Tajima's D
vcftools \
    --gzvcf "${VCF}" \
    --TajimaD 500 \
    --out "${POP_DIR}/frequencies/tajima_d" \
    2>>"${LOG}"

# Hétérozygosité par échantillon
vcftools \
    --gzvcf "${VCF}" \
    --het \
    --out "${POP_DIR}/frequencies/heterozygosity" \
    2>>"${LOG}"

echo "Statistiques de diversité calculées." | tee -a "${LOG}"

# =============================================================
# MODULE 3 — Matrice de distance + PCA
# =============================================================
echo "" | tee -a "${LOG}"
echo "--- MODULE 3 : PCA ---" | tee -a "${LOG}"

# Filtrer SNPs bialléliques uniquement pour PCA
vcftools \
    --gzvcf "${VCF}" \
    --min-alleles 2 \
    --max-alleles 2 \
    --mac 1 \
    --recode \
    --recode-INFO-all \
    --out "${POP_DIR}/pca/biallelic_snps" \
    2>>"${LOG}"

bgzip -f "${POP_DIR}/pca/biallelic_snps.recode.vcf"
tabix -p vcf "${POP_DIR}/pca/biallelic_snps.recode.vcf.gz"

N_SNPS=$(bcftools view "${POP_DIR}/pca/biallelic_snps.recode.vcf.gz" | \
         grep -v "^#" | wc -l)
echo "SNPs bialléliques pour PCA : ${N_SNPS}" | tee -a "${LOG}"

# Conversion VCF → PLINK
plink \
    --vcf "${POP_DIR}/pca/biallelic_snps.recode.vcf.gz" \
    --vcf-half-call m \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --make-bed \
    --out "${POP_DIR}/pca/plink_data" \
    2>>"${LOG}"

# PCA avec PLINK
plink \
    --bfile "${POP_DIR}/pca/plink_data" \
    --allow-extra-chr \
    --pca 10 \
    --out "${POP_DIR}/pca/pca_results" \
    2>>"${LOG}"

echo "PCA terminée." | tee -a "${LOG}"

# =============================================================
# MODULE 4 — Relatedness (IBD)
# =============================================================
echo "" | tee -a "${LOG}"
echo "--- MODULE 4 : Relatedness ---" | tee -a "${LOG}"

plink \
    --bfile "${POP_DIR}/pca/plink_data" \
    --allow-extra-chr \
    --genome \
    --out "${POP_DIR}/relatedness/ibd_matrix" \
    2>>"${LOG}"

echo "Relatedness calculé." | tee -a "${LOG}"

# =============================================================
# MODULE 5 — Haplotypes kelch13
# =============================================================
echo "" | tee -a "${LOG}"
echo "--- MODULE 5 : Haplotypes kelch13 ---" | tee -a "${LOG}"

# Extraire les génotypes kelch13 pour tous les échantillons
bcftools view \
    -r Pf3D7_13_v3:1725260-1726923 \
    "${VCF}" | \
bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
    > "${POP_DIR}/haplotypes/kelch13_genotypes.tsv" \
    2>>"${LOG}"

# Extraire la liste des samples
bcftools query \
    -l "${VCF}" \
    > "${POP_DIR}/haplotypes/sample_list.txt" \
    2>>"${LOG}"

echo "Haplotypes kelch13 extraits." | tee -a "${LOG}"

echo "" | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"
echo " Pop. genetics terminée : $(date)"       | tee -a "${LOG}"
echo "Fichiers générés :"                       | tee -a "${LOG}"
ls -lh "${POP_DIR}/pca/"*.eigenvec \
       "${POP_DIR}/frequencies/"*.frq \
       "${POP_DIR}/relatedness/"*.genome \
       2>>"${LOG}" | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"
