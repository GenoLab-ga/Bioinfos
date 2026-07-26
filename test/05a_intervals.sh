#!/bin/bash
# =============================================================
# Étape 5a : Création des intervalles cibles – 5 loci résistance
# Coordonnées Pf3D7 – PlasmoDB release 68
# =============================================================

set -euo pipefail

REF_DIR="data/reference"
INTERVALS="${REF_DIR}/resistance_loci.interval_list"
BED="${REF_DIR}/resistance_loci.bed"
LOG="logs/05a_intervals.log"

mkdir -p logs

echo "Création des intervalles cibles..." | tee "${LOG}"

# ── Fichier BED des 5 loci ────────────────────────────────────
# Format : chr  start(0-based)  end  nom
# Coordonnées Pf3D7 avec marges ±500 bp autour du CDS
cat > "${BED}" << 'EOF'
Pf3D7_04_v3	747,927	749,845	dhfr
Pf3D7_05_v3	957,880	962,149	mdr1
Pf3D7_07_v3	403,222	404,903	crt
Pf3D7_08_v3	549,993	555,749	dhps
Pf3D7_13_v3	1,725,259	1,726,923	kelch13
EOF

# ── Conversion BED → interval_list GATK ───────────────────────
gatk BedToIntervalList \
    -I data/reference/resistance_loci.bed \
    -O data/reference/resistance_loci.interval_list \
    -SD data/reference/Pf3D7.dict \
    2>&1 | tee logs/05a_intervals.log

echo "" | tee -a "${LOG}"
echo "Intervalles générés :" | tee -a "${LOG}"
grep -v "^@" "${INTERVALS}" | tee -a "${LOG}"
echo "" | tee -a "${LOG}"
echo "Taille totale couverte :" | tee -a "${LOG}"
grep -v "^@" "${INTERVALS}" | \
    awk '{sum += $2-$1} END {print sum " bp"}' | tee -a "${LOG}"

echo "Étape 5a terminée." | tee -a "${LOG}"
