#!/bin/bash
# =============================================================
# Étape 5c : GenomicsDBImport + GenotypeGVCFs + Filtrage
# Projet   : PRJNA1465284 – P. falciparum drug resistance
# =============================================================

set -euo pipefail

# ── Paramètres ───────────────────────────────────────────────
GVCF_DIR="data/gvcf_files"
DB_DIR="data/genomicsdb"
VCF_DIR="results/vcf"
REF="data/reference/Pf3D7.fasta"
INTERVALS="data/reference/resistance_loci.interval_list"
LOG="logs/05c_genotyping.log"

mkdir -p "${VCF_DIR}" logs

echo "=======================================" | tee -a "${LOG}"
echo " Génotypage démarré : $(date)"           | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"

# ── Construction arguments -V pour GenomicsDBImport ───────────
V_ARGS=""
while IFS= read -r gvcf; do
    V_ARGS="${V_ARGS} -V ${gvcf}"
done < "${GVCF_DIR}/gvcf_list.txt"

N_SAMPLES=$(wc -l < "${GVCF_DIR}/gvcf_list.txt")
echo "Échantillons à intégrer : ${N_SAMPLES}" | tee -a "${LOG}"

# ── GenomicsDBImport ──────────────────────────────────────────
echo "" | tee -a "${LOG}"
echo "GenomicsDBImport..." | tee -a "${LOG}"

rm -rf "${DB_DIR}"

gatk GenomicsDBImport \
    ${V_ARGS} \
    --genomicsdb-workspace-path "${DB_DIR}" \
    --intervals "${INTERVALS}" \
    --reader-threads 4 \
    --batch-size 50 \
    --tmp-dir /tmp \
    --verbosity ERROR \
    2>>"${LOG}" \
|| {
    echo "ERREUR : GenomicsDBImport échoué" | tee -a "${LOG}"
    exit 1
}

echo "GenomicsDBImport terminé : $(date)" | tee -a "${LOG}"

# ── GenotypeGVCFs ─────────────────────────────────────────────
echo "" | tee -a "${LOG}"
echo "GenotypeGVCFs..." | tee -a "${LOG}"

RAW_VCF="${VCF_DIR}/raw_variants.vcf.gz"

gatk GenotypeGVCFs \
    --reference "${REF}" \
    --variant "gendb://${DB_DIR}" \
    --output "${RAW_VCF}" \
    --intervals "${INTERVALS}" \
    --sample-ploidy 1 \
    --tmp-dir /tmp \
    --verbosity ERROR \
    2>>"${LOG}" \
|| {
    echo "ERREUR : GenotypeGVCFs échoué" | tee -a "${LOG}"
    exit 1
}

echo "GenotypeGVCFs terminé : $(date)" | tee -a "${LOG}"

# ── Filtrage des variants ─────────────────────────────────────
# Pas de VQSR sur amplicon → VariantFiltration par hard filters
echo "" | tee -a "${LOG}"
echo "Filtrage variants (hard filters)..." | tee -a "${LOG}"

FILTERED_VCF="${VCF_DIR}/filtered_variants.vcf.gz"

gatk VariantFiltration \
    --reference "${REF}" \
    --variant "${RAW_VCF}" \
    --output "${FILTERED_VCF}" \
    \
    --filter-expression "QD < 2.0"          --filter-name "QD2" \
    --filter-expression "FS > 60.0"         --filter-name "FS60" \
    --filter-expression "MQ < 40.0"         --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPos8" \
    --filter-expression "DP < 10"           --filter-name "LowDepth" \
    \
    --verbosity ERROR \
    2>>"${LOG}" \
|| {
    echo "ERREUR : VariantFiltration échouée" | tee -a "${LOG}"
    exit 1
}

# ── Extraction des variants PASS uniquement ───────────────────
PASS_VCF="${VCF_DIR}/pass_variants.vcf.gz"

gatk SelectVariants \
    --reference "${REF}" \
    --variant "${FILTERED_VCF}" \
    --output "${PASS_VCF}" \
    --exclude-filtered \
    2>>"${LOG}"

# ── Statistiques finales ──────────────────────────────────────
echo "" | tee -a "${LOG}"
echo "=== Statistiques VCF ===" | tee -a "${LOG}"

TOTAL_VAR=$(bcftools stats "${RAW_VCF}" | grep "^SN" | \
            grep "number of SNPs" | awk '{print $NF}')
PASS_VAR=$(bcftools stats "${PASS_VCF}"  | grep "^SN" | \
           grep "number of SNPs" | awk '{print $NF}')

echo "Variants bruts : ${TOTAL_VAR}" | tee -a "${LOG}"
echo "Variants PASS  : ${PASS_VAR}"  | tee -a "${LOG}"
echo "" | tee -a "${LOG}"

# Comptage par locus
echo "Variants PASS par locus :" | tee -a "${LOG}"
for LOCUS in kelch13 crt mdr1 dhfr dhps; do
    COUNT=$(bcftools view "${PASS_VCF}" | \
            grep -v "^#" | grep -c "${LOCUS}" 2>/dev/null || echo "0")
    echo "  ${LOCUS} : ${COUNT}" | tee -a "${LOG}"
done

echo "" | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"
echo " Génotypage terminé : $(date)"           | tee -a "${LOG}"
echo " VCF final : ${PASS_VCF}"                | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"
