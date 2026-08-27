# Classification taxonomique

echo "🚀 Début de la classification taxonomique..."

# 1. Classification taxonomique avec le classificateur téléchargé manuellement
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads results/rep_seqs.qza \
  --o-classification results/taxonomy.qza \
  --p-n-jobs 6

echo "✅ Classification terminée."

# 2. Visualisation de la taxonomie (pour vérification rapide)
qiime metadata tabulate \
  --m-input-file results/taxonomy.qza \
  --o-visualization results/taxonomy.qzv

echo "✅ Visualisation de la taxonomie générée."
