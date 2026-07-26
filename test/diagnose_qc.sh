#!/bin/bash
# diagnose_qc.sh — détail des modules FAIL par échantillon

FASTQC_DIR="results/qc/fastqc_reports"
OUT="results/qc/fail_summary.tsv"

echo -e "Sample\tModule\tStatus" > "${OUT}"

for zip in "${FASTQC_DIR}"/*_fastqc.zip; do
    SAMPLE=$(basename "${zip}" _fastqc.zip)
    unzip -p "${zip}" "*/summary.txt" 2>/dev/null | \
    while IFS=$'\t' read -r status module _; do
        if [[ "${status}" == "FAIL" ]]; then
            echo -e "${SAMPLE}\t${module}\t${status}" >> "${OUT}"
        fi
    done
done

echo "Résumé généré : ${OUT}"
echo ""
echo "=== Modules FAIL les plus fréquents ==="
cut -f2 "${OUT}" | sort | uniq -c | sort -rn

