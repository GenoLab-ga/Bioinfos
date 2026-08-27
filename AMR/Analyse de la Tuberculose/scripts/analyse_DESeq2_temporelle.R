# ============================================================
# Analyse d'Expression Différentielle - M. tuberculosis
# Effet de la Rifampicine (2 concentrations) à 3 temps
# ============================================================

# 1) Chargement des librairies
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)

# ============================================================
# A) Importation et Préparation des données
# ============================================================

# Lire le fichier counts.txt
counts_data <- read.delim("/home/karl/Documents/Projet_mtb/Mtb/comptage/counts.txt",
                          skip = 1,
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors = FALSE)

# Isoler la matrice de comptage
counts_matrix <- counts_data[, 7:ncol(counts_data)]
rownames(counts_matrix) <- counts_data$Geneid
counts_matrix <- as.matrix(counts_matrix)

# Importation du fichier metadata.csv
colData <- read.csv("/home/karl/Documents/Projet_mtb/Mtb/metadata.csv", stringsAsFactors = FALSE)
rownames(colData) <- colData$SRR

# Factoriser
colData$Condition <- factor(colData$Condition, levels = c("Control", "Treated"))
colData$Time <- factor(colData$Time, levels = c("4h", "24h", "72h"))
colData$Concentration <- as.numeric(colData$Concentration)

# Créer un groupe combiné pour l'analyse
colData$Group <- factor(paste(colData$Condition, colData$Concentration, colData$Time, sep = "_"))

# Vérification
print("=== Vérification ===")
print(colnames(counts_matrix)[1:5])
print(rownames(colData)[1:5])
print(paste("Correspondance:", all(rownames(colData) %in% colnames(counts_matrix))))

# Réordonnancement
counts_matrix <- counts_matrix[, rownames(colData)]

# Filtration
keep <- rowSums(counts_matrix >= 10) >= 9  # Au moins 9 échantillons (1/3 des 27)
counts_matrix <- counts_matrix[keep, ]
print(paste("Gènes après filtrage:", nrow(counts_matrix)))

# ============================================================
# B) Analyse DESeq2
# ============================================================

# Créer l'objet DESeqDataSet avec le groupe
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = colData,
  design = ~ Group
)

dds <- DESeq(dds)

# Transformation pour visualisation
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# ============================================================
# C) Comparaisons : 6 comparaisons spécifiques
# ============================================================

print("=== Comparaisons ===")

# Liste des comparaisons
comparisons <- list(
  # Rif 0.008 vs Control
  list(name = "Rif008_vs_Control_4h", 
       contrast = c("Group", "Treated_0.008_4h", "Control_0_4h")),
  list(name = "Rif008_vs_Control_24h", 
       contrast = c("Group", "Treated_0.008_24h", "Control_0_24h")),
  list(name = "Rif008_vs_Control_72h", 
       contrast = c("Group", "Treated_0.008_72h", "Control_0_72h")),
  
  # Rif 0.02 vs Control
  list(name = "Rif02_vs_Control_4h", 
       contrast = c("Group", "Treated_0.02_4h", "Control_0_4h")),
  list(name = "Rif02_vs_Control_24h", 
       contrast = c("Group", "Treated_0.02_24h", "Control_0_24h")),
  list(name = "Rif02_vs_Control_72h", 
       contrast = c("Group", "Treated_0.02_72h", "Control_0_72h"))
)

# Exécuter les comparaisons
results_list <- list()

for (comp in comparisons) {
  cat(paste("\nTraitement de", comp$name, "...\n"))
  
  res <- results(dds, contrast = comp$contrast, alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$Geneid <- rownames(res_df)
  res_df <- res_df %>% na.omit() %>% arrange(padj)
  
  # Statistiques
  n_up <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE)
  n_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE)
  cat(paste("  Up:", n_up, "| Down:", n_down, "\n"))
  
  results_list[[comp$name]] <- res_df
}

# ============================================================
# D) Visualisations
# ============================================================

# Créer le dossier pour les résultats
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)
dir.create("results/DEGs", showWarnings = FALSE)

# --- PCA ---
pca_data <- plotPCA(vsd, intgroup = c("Condition", "Time"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Ajouter la concentration pour la couleur
pca_data$Concentration <- colData[rownames(pca_data), "Concentration"]

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = factor(Concentration), shape = Time)) +
  geom_point(size = 4, alpha = 0.7) +
  labs(
    title = "PCA - M. tuberculosis (Rifampicine)",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Concentration (µg/mL)",
    shape = "Temps"
  ) +
  scale_color_manual(values = c("0" = "gray30", "0.008" = "#0072B2", "0.02" = "#D55E00")) +
  scale_shape_manual(values = c("4h" = 16, "24h" = 17, "72h" = 15)) +
  theme_bw() +
  theme(legend.position = "right")

ggsave("results/figures/PCA_Rifampicine.png", plot = pca_plot, width = 9, height = 6, dpi = 300)

# --- Volcano plots ---
create_volcano <- function(res_df, title, padj_thresh = 0.05, lfc_thresh = 1) {
  volcano_df <- res_df %>%
    mutate(
      Expression = case_when(
        padj < padj_thresh & log2FoldChange > lfc_thresh  ~ "Up",
        padj < padj_thresh & log2FoldChange < -lfc_thresh ~ "Down", 
        TRUE ~ "NS"
      ),
      negLogPval = -log10(padj)
    )
  
  n_up <- sum(volcano_df$Expression == "Up", na.rm = TRUE)
  n_down <- sum(volcano_df$Expression == "Down", na.rm = TRUE)
  
  ggplot(volcano_df, aes(x = log2FoldChange, y = negLogPval, color = Expression)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NS" = "gray60")) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "gray50") +
    labs(
      title = title,
      subtitle = paste("Up:", n_up, "| Down:", n_down, "| Seuils: padj <", padj_thresh, ", |LFC| >", lfc_thresh),
      x = expression(Log[2]~"Fold Change"),
      y = expression(-Log[10]~"P-value ajustée")
    ) +
    theme_bw() +
    theme(plot.title = element_text(size = 12))
}

for (comp_name in names(results_list)) {
  volcano <- create_volcano(results_list[[comp_name]], comp_name)
  ggsave(paste0("results/figures/Volcano_", comp_name, ".png"), 
         plot = volcano, width = 7, height = 7, dpi = 300)
}

# --- Heatmaps ---
create_heatmap <- function(res_df, vsd_obj, colData_df, group1, group2, title, n_genes = 20) {
  top_genes <- res_df %>% head(n_genes)
  top_gene_ids <- top_genes$Geneid
  
  samples <- rownames(colData_df)[colData_df$Group %in% c(group1, group2)]
  
  heatmap_matrix <- assay(vsd_obj)[top_gene_ids, samples]
  heatmap_matrix <- t(scale(t(heatmap_matrix)))
  
  annotation_col <- as.data.frame(colData_df[samples, c("Condition", "Time")])
  colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
  
  pheatmap(
    heatmap_matrix,
    color = colors,
    cluster_rows = TRUE,
    show_rownames = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_col,
    main = paste("Top", n_genes, "DEGs -", title),
    fontsize_row = 7,
    fontsize_col = 8
  )
}

for (comp in comparisons) {
  png(paste0("results/figures/Heatmap_", comp$name, ".png"), 
      width = 800, height = 1000, res = 100)
  
  group1 <- comp$contrast[2]
  group2 <- comp$contrast[3]
  
  create_heatmap(results_list[[comp$name]], vsd, colData, group1, group2, comp$name)
  dev.off()
}

# ============================================================
# E) Sauvegarde des résultats
# ============================================================

for (comp_name in names(results_list)) {
  write.csv(results_list[[comp_name]], 
            paste0("/home/karl/Documents/Projet_mtb/Mtb/results/results/DEGs/DEGs_", comp_name, ".csv"), 
            row.names = FALSE)
}

# Résumé
summary_df <- data.frame(
  Comparison = names(results_list),
  Total_DEGs_padj005 = sapply(results_list, function(x) sum(x$padj < 0.05, na.rm = TRUE)),
  Up_regulated_LFC1 = sapply(results_list, function(x) sum(x$padj < 0.05 & x$log2FoldChange > 1, na.rm = TRUE)),
  Down_regulated_LFC1 = sapply(results_list, function(x) sum(x$padj < 0.05 & x$log2FoldChange < -1, na.rm = TRUE))
)

write.csv(summary_df, "/home/karl/Documents/Projet_mtb/Mtb/results/results/Summary_DEGs.csv", row.names = FALSE)
print(summary_df)

print("=== Analyse terminée ! ===")
print("Résultats sauvegardés dans le dossier 'results/'")