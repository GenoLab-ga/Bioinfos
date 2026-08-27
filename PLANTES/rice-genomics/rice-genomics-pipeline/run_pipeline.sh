#!/bin/bash
################################################################################
# RICE GENOMICS WGS PIPELINE — run_pipeline.sh
# Script maître — exécute le pipeline complet dans l'ordre
#
# USAGE :
#   bash run_pipeline.sh           — pipeline complet
#   bash run_pipeline.sh 3 5       — étapes 3 à 5 seulement
#   bash run_pipeline.sh 4         — étape 4 seulement
#
# ÉTAPES :
#   1 — Setup (génome de référence, indexation)
#   2 — Téléchargement des données (ENA)
#   3 — QC et preprocessing (fastp)
#   4 — Mapping (BWA-MEM)
#   5 — Variant calling (bcftools)
#   6 — Annotation (SnpEff)
#   7 — Analyses (PCA, arbre, Fst, π, Manhattan)
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

START_STEP=${1:-1}
END_STEP=${2:-7}

echo "================================================================"
echo "  RICE GENOMICS WGS PIPELINE"
echo "  Étapes : $START_STEP → $END_STEP"
echo "  Démarrage : $(date '+%Y-%m-%d %H:%M:%S')"
echo "================================================================"

for step in $(seq "$START_STEP" "$END_STEP"); do
    script="$SCRIPT_DIR/0${step}_*.sh"
    script=$(ls $script 2>/dev/null | head -1)

    if [ -z "$script" ]; then
        echo "[WARN] Étape $step : script introuvable"
        continue
    fi

    echo ""
    echo ">>> Lancement étape $step : $(basename $script)"
    echo "    Début : $(date '+%H:%M:%S')"
    SECONDS=0

    bash "$script"

    echo "    Durée : ${SECONDS}s"
    echo "    Fin   : $(date '+%H:%M:%S')"
done

echo ""
echo "================================================================"
echo "  PIPELINE TERMINÉ"
echo "  Fin : $(date '+%Y-%m-%d %H:%M:%S')"
echo "================================================================"
