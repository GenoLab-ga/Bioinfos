#!/usr/bin/env Rscript

# ====================================================================
# SEURAT QC: Plasmodium falciparum Gamétocytes
# ====================================================================
# CORRECTIONS:
#   1. Filtrage génique: features= attend des noms, pas un vecteur logique
#   2. GetAssayData: slot= → layer= (Seurat v5)
#   3. Détection automatique version Seurat pour compatibilité v4/v5
# ====================================================================

library(Seurat)
library(ggplot2)

# Détecter la version Seurat
seurat_version <- packageVersion("Seurat")
cat(sprintf("Version Seurat détectée: %s\n", seurat_version))
is_v5 <- as.numeric(major(seurat_version)) >= 5

# ====================================================================
# 0. VARIABLES
# ====================================================================

project_root <- "~/Documents/github_projet/Plasmodium"
data_dir     <- file.path(project_root, "quantification/matrices")
output_dir   <- file.path(project_root, "analysis/01_qc")
log_dir      <- file.path(project_root, "logs")

# Seuils QC (à ajuster selon vos données)
min_features <- 200
max_features <- 5000
min_counts   <- 500      # abaissé de 1000 pour données bulk/pseudo-bulk
max_counts   <- 50000
max_mt_pct   <- 15
min_cells_per_gene <- 3

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir,    showWarnings = FALSE, recursive = TRUE)

cat("======================================================\n")
cat("SEURAT QC & PREPROCESSING\n")
cat("======================================================\n\n")

# ====================================================================
# 1. CHARGEMENT
# ====================================================================

cat("======================================================\n")
cat("ÉTAPE 1: Chargement des données\n")
cat("======================================================\n\n")

matrix_file <- file.path(data_dir, "combined_counts_matrix.csv")
if (!file.exists(matrix_file)) {
    cat("❌ Fichier matrice non trouvé:", matrix_file, "\n")
    quit(status = 1)
}

cat(sprintf("Chargement: %s\n", matrix_file))
counts_df <- read.csv(matrix_file, row.names = 1, check.names = FALSE)
counts_matrix <- as.matrix(counts_df)

cat(sprintf("  Dimensions brutes: %d gènes × %d échantillons\n",
            nrow(counts_matrix), ncol(counts_matrix)))

# Avertissement si peu d'échantillons (mode bulk vs single-cell)
if (ncol(counts_matrix) < 10) {
    cat(sprintf("\n  ⚠ ATTENTION: %d colonnes seulement.\n", ncol(counts_matrix)))
    cat("    Si ce sont des échantillons BULK (pas des cellules individuelles),\n")
    cat("    les métriques Seurat par 'cellule' reflètent les échantillons,\n")
    cat("    pas des vraies cellules. Interprétez les plots en conséquence.\n\n")
}

cat("\n")

# Créer objet Seurat
cat("Création objet Seurat...\n")
seurat_obj <- CreateSeuratObject(
    counts      = counts_matrix,
    min.cells   = min_cells_per_gene,
    min.features = min_features
)

cat(sprintf("  ✓ Après filtres initiaux: %d gènes × %d cellules/échantillons\n\n",
            nrow(seurat_obj), ncol(seurat_obj)))

# ====================================================================
# 2. MÉTRIQUES QC
# ====================================================================

cat("======================================================\n")
cat("ÉTAPE 2: Calcul des métriques QC\n")
cat("======================================================\n\n")

# Gènes mitochondriaux (Plasmodium: chromosome mitochondrial = chrM ou PF3D7_MIT)
mt_pattern <- "^PF3D7_M"
mt_genes <- grep(mt_pattern, rownames(seurat_obj), value = TRUE)
cat(sprintf("Gènes mitochondriaux trouvés (pattern '%s'): %d\n\n", mt_pattern, length(mt_genes)))

seurat_obj <- PercentageFeatureSet(seurat_obj,
                                   pattern  = mt_pattern,
                                   col.name = "percent.mt")

seurat_obj$log10_counts   <- log10(seurat_obj$nCount_RNA + 1)
seurat_obj$log10_features <- log10(seurat_obj$nFeature_RNA + 1)

# Résumé statistiques
stats_summary <- function(x) sprintf("Min=%.0f | Médian=%.0f | Max=%.0f", min(x), median(x), max(x))
cat("Statistiques QC (pré-filtrage):\n")
cat(sprintf("  nCount_RNA   : %s\n", stats_summary(seurat_obj$nCount_RNA)))
cat(sprintf("  nFeature_RNA : %s\n", stats_summary(seurat_obj$nFeature_RNA)))
cat(sprintf("  percent.mt   : médian=%.2f%%\n\n", median(seurat_obj$percent.mt)))

# ====================================================================
# 3. VISUALISATION PRÉ-FILTRAGE
# ====================================================================

cat("======================================================\n")
cat("ÉTAPE 3: Visualisation pré-filtrage\n")
cat("======================================================\n\n")

pdf(file.path(output_dir, "01_qc_pre_filtering.pdf"), width = 14, height = 10)

tryCatch({
    p1 <- VlnPlot(seurat_obj,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  ncol     = 3,
                  pt.size  = 0.5)
    print(p1)
}, error = function(e) cat("⚠ VlnPlot:", conditionMessage(e), "\n"))

tryCatch({
    p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    print(p2 + p3)
}, error = function(e) cat("⚠ FeatureScatter:", conditionMessage(e), "\n"))

dev.off()
cat(sprintf("✓ PDF sauvegardé: %s/01_qc_pre_filtering.pdf\n\n", output_dir))

# ====================================================================
# 4. FILTRAGE
# ====================================================================

cat("======================================================\n")
cat("ÉTAPE 4: Filtrage\n")
cat("======================================================\n\n")

cat(sprintf("Seuils cellulaires:\n"))
cat(sprintf("  nFeature_RNA : [%d, %d]\n", min_features, max_features))
cat(sprintf("  nCount_RNA   : [%d, %d]\n", min_counts, max_counts))
cat(sprintf("  percent.mt   : < %.1f%%\n\n", max_mt_pct))

cells_before <- ncol(seurat_obj)

seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA >= min_features &
                              nFeature_RNA <= max_features &
                              nCount_RNA   >= min_counts   &
                              nCount_RNA   <= max_counts   &
                              percent.mt   <  max_mt_pct)

cells_after <- ncol(seurat_obj)
cat(sprintf("Filtrage cellulaire: %d → %d (%.1f%% retenus)\n\n",
            cells_before, cells_after,
            ifelse(cells_before > 0, 100 * cells_after / cells_before, 0)))

if (cells_after == 0) {
    cat("❌ Aucune cellule/échantillon ne passe les filtres.\n")
    cat("   Allégez les seuils (min_counts, min_features) en haut du script.\n")
    quit(status = 1)
}

# ====================================================================
# CORRECTION BUG 4: filtrage génique
# features= attend un vecteur de noms, pas un vecteur logique
# ====================================================================

genes_before <- nrow(seurat_obj)

# Récupérer la matrice de counts (compatible v4 et v5)
counts_mat <- tryCatch(
    GetAssayData(seurat_obj, assay = "RNA", layer = "counts"),   # Seurat v5
    error = function(e)
        GetAssayData(seurat_obj, assay = "RNA", slot = "counts") # Seurat v4
)

# Calculer quels gènes sont détectés dans > min_cells_per_gene cellules
genes_to_keep <- names(which(rowSums(counts_mat > 0) > min_cells_per_gene))

# Filtrer sur les noms de gènes (CORRECTION)
seurat_obj <- subset(seurat_obj, features = genes_to_keep)

genes_after <- nrow(seurat_obj)
cat(sprintf("Filtrage génique (>%d cellules): %d → %d gènes (%.1f%% retenus)\n\n",
            min_cells_per_gene, genes_before, genes_after,
            ifelse(genes_before > 0, 100 * genes_after / genes_before, 0)))

# ====================================================================
# 5. VISUALISATION POST-FILTRAGE
# ====================================================================

cat("======================================================\n")
cat("ÉTAPE 5: Visualisation post-filtrage\n")
cat("======================================================\n\n")

pdf(file.path(output_dir, "02_qc_post_filtering.pdf"), width = 12, height = 6)
tryCatch({
    p_post <- VlnPlot(seurat_obj,
                      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                      ncol     = 3,
                      pt.size  = 0.5)
    print(p_post)
}, error = function(e) cat("⚠ VlnPlot post:", conditionMessage(e), "\n"))
dev.off()

cat(sprintf("✓ PDF sauvegardé: %s/02_qc_post_filtering.pdf\n\n", output_dir))

# ====================================================================
# 6. NORMALISATION
# ====================================================================

cat("======================================================\n")
cat("ÉTAPE 6: Normalisation (LogNormalize)\n")
cat("======================================================\n\n")

seurat_obj <- NormalizeData(seurat_obj,
                            normalization.method = "LogNormalize",
                            scale.factor         = 10000)

cat("✓ Normalisation complétée\n\n")

# ====================================================================
# 7. SAUVEGARDE
# ====================================================================

cat("======================================================\n")
cat("ÉTAPE 7: Sauvegarde\n")
cat("======================================================\n\n")

output_rds <- file.path(output_dir, "seurat_qc_filtered.rds")
saveRDS(seurat_obj, file = output_rds)
cat(sprintf("✓ Objet Seurat: %s\n", output_rds))

# Récupérer les données normalisées — compatible v4/v5
# CORRECTION BUG 5: slot= → layer= pour Seurat v5
norm_data <- tryCatch(
    GetAssayData(seurat_obj, assay = "RNA", layer = "data"),   # Seurat v5
    error = function(e)
        GetAssayData(seurat_obj, assay = "RNA", slot  = "data") # Seurat v4
)

output_csv <- file.path(output_dir, "normalized_counts.csv")
write.csv(as.matrix(norm_data), file = output_csv)
cat(sprintf("✓ Matrice normalisée: %s\n\n", output_csv))

# ====================================================================
# 8. RAPPORT FINAL
# ====================================================================

cat("======================================================\n")
cat("RAPPORT FINAL\n")
cat("======================================================\n\n")

cat(sprintf("Dimensions finales: %d gènes × %d cellules/échantillons\n\n",
            nrow(seurat_obj), ncol(seurat_obj)))

cat("Fichiers générés:\n")
cat(sprintf("  %s\n", output_rds))
cat(sprintf("  %s\n", output_csv))
cat(sprintf("  %s/01_qc_pre_filtering.pdf\n",  output_dir))
cat(sprintf("  %s/02_qc_post_filtering.pdf\n", output_dir))
cat("\n")

cat(sprintf("[%s] Seurat QC completed\n", Sys.time()),
    file   = file.path(log_dir, "analysis.log"),
    append = TRUE)

cat("======================================================\n")
cat("✓ ÉTAPE 5 COMPLÉTÉE\n")
cat("======================================================\n\n")
cat("PROCHAINES ÉTAPES:\n")
cat("  1. Vérifier les PDFs\n")
cat("  2. Ajuster les seuils si trop de cellules filtrées\n")
cat("  3. Clustering: FindVariableFeatures → ScaleData → RunPCA → FindClusters\n\n")
