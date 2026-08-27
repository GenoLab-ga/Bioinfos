#!/usr/bin/env bash
set -euo pipefail

# Création d'un environnement isolé avec les outils modernes
mamba create -n vina-pymol -c conda-forge -c bioconda \
    python=3.12 \
    vina \
    meeko \
    rdkit \
    pymol-open-source \
    pandas -y

mamba activate vina-pymol
echo "Environnement activé avec succès."
