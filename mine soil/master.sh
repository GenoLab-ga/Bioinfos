#!/bin/bash
# ==============================================================================
# MASTER PIPELINE : Analyse du Microbiome Édaphique (Proof of Concept)
# Auteur : Karl
# Description : Orchestre l'ensemble du pipeline, du téléchargement SRA à la visualisation R.
# ==============================================================================

# Mode strict : le script s'arrête immédiatement si une commande échoue
set -e 

echo "=========================================================="
echo "🚀 DÉMARRAGE DU PIPELINE MASTER D'ANALYSE MICROBIOME"
echo "=========================================================="
echo "Date de début : $(date)"
echo "Répertoire de travail : $(pwd)"

# 1. Vérification et activation de l'environnement Mamba
echo -e "\n[INFO] Vérification de l'environnement Mamba..."
if [[ "$CONDA_DEFAULT_ENV" != "qiime2-amplicon-2024.10" ]]; then
    echo "⚠️ Activation de l'environnement 'qiime2-amplicon-2024.10'..."
    # Source mamba/conda (ajustez le chemin si nécessaire)
    if [ -f ~/miniforge3/etc/profile.d/conda.sh ]; then
        source ~/miniforge3/etc/profile.d/conda.sh
    elif [ -f ~/miniconda3/etc/profile.d/conda.sh ]; then
        source ~/miniconda3/etc/profile.d/conda.sh
    fi
    mamba activate qiime2-amplicon-2024.10
fi
echo "✅ Environnement actif : $CONDA_DEFAULT_ENV"

# 2. Téléchargement des données brutes
echo -e "\n[1/6]  Étape 1 : Téléchargement des données brutes (SRA)..."
if [ -f "01_script_download.sh" ]; then
    bash 01_script_download.sh
    echo "✅ Étape 1 terminée."
else
    echo "❌ Erreur: 01_script_download.sh introuvable!"
    exit 1
fi

# 3. Génération dynamique du manifeste
echo -e "\n[2/6]  Étape 2 : Génération du manifeste QIIME 2..."
if [ -f "02_generate_manifest.py" ]; then
    python3 02_generate_manifest.py
    echo "✅ Étape 2 terminée."
else
    echo "❌ Erreur: 02_generate_manifest.py introuvable!"
    exit 1
fi

# 4. Traitement QIIME 2 (Import + DADA2)
echo -e "\n[3/6]  Étape 3 : Traitement QIIME 2 (Import + DADA2)..."
if [ -f "03_qiime_analyse.sh" ]; then
    bash 03_qiime_analyse.sh
    echo "✅ Étape 3 terminée."
else
    echo "❌ Erreur: 03_qiime_analyse.sh introuvable!"
    exit 1
fi

# 5. Classification taxonomique et Export
echo -e "\n[4/6] 🏷️ Étape 4 : Classification taxonomique et préparation pour R..."
if [ -f "03_taxonomie.sh" ]; then
    bash 03_taxonomie.sh
    echo "✅ Étape 4 terminée."
else
    echo "❌ Erreur: 03_taxonomie.sh introuvable!"
    exit 1
fi

# 6. Analyse statistique et visualisation avec R (via phyloseq.sh)
echo -e "\n[5/6] 📊 Étape 5 : Analyse phyloseq (shell wrapper)..."
if [ -f "04_phyloseq.sh" ]; then
    bash 04_phyloseq.sh
    echo "✅ Étape 5 terminée."
else
    echo "❌ Erreur: 04_phyloseq.sh introuvable!"
    exit 1
fi

# 7. Exécution directe du script R (si 04_phyloseq.sh ne le fait pas)
echo -e "\n[6/6] 📈 Étape 6 : Exécution du script R phyloseq..."
if [ -f "05_phyloseq.R" ]; then
    Rscript 05_phyloseq.R
    echo "✅ Étape 6 terminée."
else
    echo "❌ Erreur: 05_phyloseq.R introuvable!"
    exit 1
fi

# Fin du pipeline
echo -e "\n=========================================================="
echo "🎉 PIPELINE TERMINÉ AVEC SUCCÈS !"
echo "📅 Date de fin : $(date)"
echo "📁 Résultats disponibles dans le dossier 'results/'"
echo "   - Graphiques PNG/PDF"
echo "   - Tables d'abondance"
echo "   - Taxonomie"
echo "=========================================================="

# Résumé rapide
echo -e "\n📊 RÉSUMÉ DU PIPELINE:"
echo "   - Données brutes: $(ls -1 data/raw/*.fastq 2>/dev/null | wc -l) fichiers FASTQ"
echo "   - Artefacts QIIME 2: $(ls -1 results/*.qza 2>/dev/null | wc -l) fichiers"
echo "   - Visualisations: $(ls -1 results/*.qzv 2>/dev/null | wc -l) fichiers"
echo "   - Graphiques R: $(ls -1 results/*.png 2>/dev/null | wc -l) fichiers PNG"
