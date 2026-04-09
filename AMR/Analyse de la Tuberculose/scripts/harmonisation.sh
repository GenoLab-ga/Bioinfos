#!/bin/bash

# Nettoyer les noms de colonnes dans counts.txt
sed 's|mapping/||g; s|.sorted.bam||g' comptage/counts.txt > comptage/counts_clean.txt

# Vérifier
echo "=== Avant ==="
head -2 comptage/counts.txt | tail -1 | cut -f7-10

echo ""
echo "=== Après ==="
head -2 comptage/counts_clean.txt | tail -1 | cut -f7-10

# Remplacer l'original
mv comptage/counts_clean.txt comptage/counts.txt

echo ""
echo "✅ Fichier nettoyé !"
