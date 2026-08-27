# 1. Importation des données Paired-End
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path data/manifest/paired_end_manifest.csv \
  --output-path data/manifest/paired_end_demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

# 2. Résumé du contrôle qualité (génère un fichier .qzv interactif)
qiime demux summarize \
  --i-data data/manifest/paired_end_demux.qza \
  --o-visualization results/demux_summary.qzv


# Dénisage avec troncatures réduites pour améliorer la fusion
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs data/manifest/paired_end_demux.qza \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 200 \
  --p-n-threads 6 \
  --p-n-reads-learn 1000000 \
  --p-min-overlap 12 \
  --o-representative-sequences results/rep_seqs.qza \
  --o-table results/table.qza \
  --o-denoising-stats results/denoising_stats_v2.qza


