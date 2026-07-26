#!/bin/bash
# =============================================================
# Étape 6 : Annotation SnpEff – version corrigée
# Projet  : PRJNA1465284 – P. falciparum drug resistance
# =============================================================

set -euo pipefail

VCF_IN="results/vcf/pass_variants.vcf.gz"
ANNOT_DIR="results/annotation"
DB="Plasmodium_falciparum"
LOG="logs/06_annotation.log"

mkdir -p "${ANNOT_DIR}" logs

echo "=======================================" | tee "${LOG}"
echo " Annotation démarrée : $(date)"          | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"

# ── Étape 1 : Renommer les chromosomes Pf3D7_XX_v3 → 1..14 ──
# SnpEff Plasmodium_falciparum utilise les noms courts 1-14
echo "Renommage chromosomes pour SnpEff..." | tee -a "${LOG}"

cat > /tmp/chr_rename.txt << 'EOF'
Pf3D7_01_v3	1
Pf3D7_02_v3	2
Pf3D7_03_v3	3
Pf3D7_04_v3	4
Pf3D7_05_v3	5
Pf3D7_06_v3	6
Pf3D7_07_v3	7
Pf3D7_08_v3	8
Pf3D7_09_v3	9
Pf3D7_10_v3	10
Pf3D7_11_v3	11
Pf3D7_12_v3	12
Pf3D7_13_v3	13
Pf3D7_14_v3	14
Pf3D7_API_v3	API
Pf3D7_MIT_v3	MIT
EOF

bcftools annotate \
    --rename-chrs /tmp/chr_rename.txt \
    "${VCF_IN}" \
    -O z -o "${ANNOT_DIR}/renamed_variants.vcf.gz" \
    2>>"${LOG}"

tabix -p vcf "${ANNOT_DIR}/renamed_variants.vcf.gz"

echo "Renommage terminé." | tee -a "${LOG}"

# ── Vérifier les noms après renommage ────────────────────────
echo "Chromosomes dans VCF renommé :" | tee -a "${LOG}"
bcftools view "${ANNOT_DIR}/renamed_variants.vcf.gz" | \
    grep -v "^#" | awk '{print $1}' | sort -u | tee -a "${LOG}"

# ── Étape 2 : Annotation SnpEff ──────────────────────────────
echo "" | tee -a "${LOG}"
echo "Annotation SnpEff..." | tee -a "${LOG}"

snpEff ann \
    -v \
    -noStats \
    "${DB}" \
    "${ANNOT_DIR}/renamed_variants.vcf.gz" \
    2>>"${LOG}" \
| bgzip > "${ANNOT_DIR}/annotated_variants.vcf.gz"

tabix -p vcf "${ANNOT_DIR}/annotated_variants.vcf.gz"

echo "Annotation terminée : $(date)" | tee -a "${LOG}"

# ── Étape 3 : Extraction TSV via INFO/ANN (pas CSQ) ──────────
echo "" | tee -a "${LOG}"
echo "Extraction TSV depuis INFO/ANN..." | tee -a "${LOG}"

# En-tête
echo -e "CHROM\tPOS\tREF\tALT\tEFFECT\tIMPACT\tGENE\tHGVS_c\tHGVS_p" \
    > "${ANNOT_DIR}/variants_annotated.tsv"

# Parser le champ ANN directement avec awk
bcftools view "${ANNOT_DIR}/annotated_variants.vcf.gz" | \
    grep -v "^#" | \
    awk 'BEGIN{OFS="\t"} {
        # Extraire le champ ANN de la colonne INFO
        match($8, /ANN=([^;]+)/, arr)
        if (arr[1] == "") next
        # Prendre la première annotation (peut être multiple)
        split(arr[1], anns, ",")
        split(anns[1], f, "|")
        # f[2]=effet f[3]=impact f[4]=gène f[10]=HGVSc f[11]=HGVSp
        print $1, $2, $4, $5, f[2], f[3], f[4], f[10], f[11]
    }' >> "${ANNOT_DIR}/variants_annotated.tsv"

echo "TSV généré." | tee -a "${LOG}"
wc -l "${ANNOT_DIR}/variants_annotated.tsv" | tee -a "${LOG}"

# ── Étape 4 : Filtrage missense ───────────────────────────────
echo "" | tee -a "${LOG}"
echo "Filtrage variants missense..." | tee -a "${LOG}"

{
    head -1 "${ANNOT_DIR}/variants_annotated.tsv"
    grep -E "missense_variant|stop_gained|frameshift_variant" \
        "${ANNOT_DIR}/variants_annotated.tsv" || true
} > "${ANNOT_DIR}/missense_variants.tsv"

N_MISSENSE=$(( $(wc -l < "${ANNOT_DIR}/missense_variants.tsv") - 1 ))
echo "Missense variants : ${N_MISSENSE}" | tee -a "${LOG}"

# ── Étape 5 : Résumé par gène ─────────────────────────────────
echo "" | tee -a "${LOG}"
echo "=== Missense par gène ===" | tee -a "${LOG}"
awk 'NR>1 && $7!="" {print $7}' \
    "${ANNOT_DIR}/missense_variants.tsv" | \
    sort | uniq -c | sort -rn | tee -a "${LOG}"

# ── Étape 6 : Mutations WHO connues ───────────────────────────
echo "" | tee -a "${LOG}"
echo "=== Mutations WHO résistance ===" | tee -a "${LOG}"

declare -A WHO_MUTS=(
    ["kelch13"]="622 441 580 539"
    ["crt"]="76"
    ["mdr1"]="86 184 1246"
    ["dhfr"]="108 51 59 164"
    ["dhps"]="437 540 581"
)

FOUND=0
while IFS= read -r line; do
    [[ "${line}" == CHROM* ]] && continue
    GENE=$(echo "${line}" | awk '{print $7}')
    HGVSP=$(echo "${line}" | awk '{print $9}')
    POS=$(echo "${line}" | awk '{print $2}')

    for GNAME in "${!WHO_MUTS[@]}"; do
        if echo "${GENE}" | grep -qi "${GNAME}"; then
            for CODON in ${WHO_MUTS[$GNAME]}; do
                if echo "${HGVSP}" | grep -q "${CODON}"; then
                    echo "  ✓ ${GENE} — pos ${POS} — ${HGVSP}" \
                        | tee -a "${LOG}"
                    FOUND=$((FOUND + 1))
                fi
            done
        fi
    done
done < "${ANNOT_DIR}/missense_variants.tsv"

[ "${FOUND}" -eq 0 ] && \
    echo "  Vérifier manuellement le TSV — noms gènes à confirmer" \
    | tee -a "${LOG}"

# ── Nettoyage ─────────────────────────────────────────────────
mv snpEff_summary.html \
    "${ANNOT_DIR}/snpEff_summary_PRJNA1465284.html" 2>/dev/null || true
mv snpEff_genes.txt \
    "${ANNOT_DIR}/snpEff_genes_PRJNA1465284.txt" 2>/dev/null || true

echo "" | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"
echo " Annotation terminée : $(date)"          | tee -a "${LOG}"
echo "=======================================" | tee -a "${LOG}"
