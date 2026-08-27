echo "📦 Préparation des exports pour R..."
mkdir -p results/exports

# Export de la table d'abondance (format BIOM)
qiime tools export \
  --input-path results/table.qza \
  --output-path results/exports/table

# Export de la taxonomie
qiime tools export \
  --input-path results/taxonomy.qza \
  --output-path results/exports/taxonomy

# Conversion du fichier BIOM en TSV (format lisible par R)
biom convert \
  -i results/exports/table/feature-table.biom \
  -o results/exports/feature-table.tsv \
  --to-tsv

echo "🎉 Pipeline QIIME 2 terminé avec succès !"
echo "📁 Les fichiers pour l'analyse R sont dans : results/exports/"
