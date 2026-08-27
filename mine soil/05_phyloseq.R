#!/usr/bin/env Rscript
# ==============================================================================
# Pipeline d'analyse microbiome avec phyloseq
# Projet : Biosurveillance des sols miniers (Cas d'étude COMILOG)
# Auteur : Karl (Ingénieur en Bioinformatique)
# ==============================================================================

# 1. GESTION DES DÉPENDANCES ----------------------------------------------------
cat("📦 Vérification des dépendances...\n")

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Ajout de 'stringr' ici pour la fonction str_split_fixed
cran_packages <- c("vegan", "ggplot2", "dplyr", "tidyr", "readr", "scales", "gridExtra", "stringr")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE, quietly = TRUE)
}

if (!require("phyloseq", character.only = TRUE, quietly = TRUE)) {
  BiocManager::install("phyloseq", update = FALSE, ask = FALSE)
}
library(phyloseq, quietly = TRUE)

cat("✅ Toutes les dépendances sont chargées.\n\n")


# 2. IMPORTATION DES DONNÉES ----------------------------------------------------
cat("📥 Importation des données QIIME 2...\n")

# a) Table d'abondance
otu_table_raw <- read_tsv("results/exports/feature-table.tsv", skip = 1, col_names = TRUE, show_col_types = FALSE)
colnames(otu_table_raw)[1] <- "FeatureID"
otu_matrix <- as.matrix(otu_table_raw[, -1])
rownames(otu_matrix) <- otu_table_raw$FeatureID

# b) Taxonomie (Nettoyage robuste pour éviter les erreurs tax_glom)
tax_df <- read_tsv("results/exports/taxonomy/taxonomy.tsv", col_names = TRUE, show_col_types = FALSE)
colnames(tax_df)[1] <- "FeatureID"

# 1. Uniformiser le séparateur (remplacer "; " par ";")
tax_df$Taxon <- gsub("; ", ";", tax_df$Taxon)

# 2. Séparer en 7 niveaux
tax_split <- str_split_fixed(tax_df$Taxon, ";", 7)
colnames(tax_split) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# 3. Nettoyer les préfixes (gère "p__", "D_1__", etc.) et remplacer les vides par "Unassigned"
tax_clean <- as.data.frame(apply(tax_split, 2, function(x) {
  x <- gsub("^[a-zA-Z0-9_]+__", "", x)       # Enlève k__, p__, D_0__, etc.
  x[x == "" | is.na(x)] <- "Unassigned"      # Empêche tax_glom de supprimer le taxon
  return(x)
}))
rownames(tax_clean) <- tax_df$FeatureID
tax_ps <- tax_table(as.matrix(tax_clean))

# c) Métadonnées
meta_df <- read_tsv("metadata.tsv", col_types = cols(), show_col_types = FALSE)
meta_df <- as.data.frame(meta_df)
rownames(meta_df) <- meta_df$`sample-id`

# 3. ALIGNEMENT ET CRÉATION DE L'OBJET PHYLOSEQ ---------------------------------
cat("🔍 Alignement des échantillons entre OTU et métadonnées...\n")

# Récupérer les noms d'échantillons des deux sources
samples_otu <- colnames(otu_matrix)
samples_meta <- rownames(meta_df)

# Trouver l'intersection stricte
common_samples <- intersect(samples_otu, samples_meta)

if (length(common_samples) < 2) {
  stop("❌ Erreur critique : Aucun ou trop peu d'échantillons en commun.\n",
       "Vérifiez que les noms dans 'metadata.tsv' correspondent EXACTEMENT aux colonnes du feature-table.tsv.")
}

cat("✅", length(common_samples), "échantillons communs trouvés et alignés.\n")

# Filtrer la matrice OTU et les métadonnées pour ne garder que l'intersection (dans le même ordre)
otu_matrix <- otu_matrix[, common_samples, drop = FALSE]
meta_df <- meta_df[common_samples, , drop = FALSE]

# Création des objets phyloseq
otu_ps <- otu_table(otu_matrix, taxa_are_rows = TRUE)
meta_ps <- sample_data(meta_df)

# Assemblage final
ps <- phyloseq(otu_ps, tax_ps, meta_ps)
cat("✅ Objet phyloseq créé avec succès :", nsamples(ps), "échantillons,", ntaxa(ps), "ASVs.\n")

# Filtrage des ASVs rares
ps_filt <- filter_taxa(ps, function(x) sum(x) > 10, TRUE)
cat("✅ Après filtrage des ASVs rares :", ntaxa(ps_filt), "ASVs conservés.\n\n")


# 3. CRÉATION DE L'OBJET PHYLOSEQ -----------------------------------------------
ps <- phyloseq(otu_ps, tax_ps, meta_ps)
cat("✅ Objet phyloseq créé :", nsamples(ps), "échantillons,", ntaxa(ps), "ASVs.\n")

ps_filt <- filter_taxa(ps, function(x) sum(x) > 10, TRUE)
cat("✅ Après filtrage :", ntaxa(ps_filt), "ASVs conservés.\n")

# 4. ANALYSE DE DIVERSITÉ ALPHA -------------------------------------------------
cat("📊 Calcul de la diversité Alpha...\n")

# Calcul des indices (vegan::diversity pour Shannon, vegan::specnumber pour Richesse observée)
alpha_df <- data.frame(
  sample_data(ps_filt),
  Shannon = diversity(t(otu_table(ps_filt)), index = "shannon"),
  Observed = specnumber(t(otu_table(ps_filt)))
)

# Test statistique : Kruskal-Wallis (non-paramétrique, adapté aux petits effectifs)
kw_shannon <- kruskal.test(Shannon ~ status, data = alpha_df)
kw_obs <- kruskal.test(Observed ~ status, data = alpha_df)

cat("   Test Kruskal-Wallis (Shannon) : p-value =", format.pval(kw_shannon$p.value, digits = 3), "\n")
cat("   Test Kruskal-Wallis (Observed): p-value =", format.pval(kw_obs$p.value, digits = 3), "\n")

# Graphique Diversity Alpha
p_alpha <- ggplot(alpha_df, aes(x = status, y = Shannon, fill = status)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Diversité Alpha (Indice de Shannon)",
       x = "Statut du site", y = "Indice de Shannon") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

# 5. ANALYSE DE DIVERSITÉ BÊTA --------------------------------------------------
cat("📊 Calcul de la diversité Bêta (Bray-Curtis)...\n")

ps_rel <- transform_sample_counts(ps_filt, function(x) x / sum(x))
dist_bc <- phyloseq::distance(ps_rel, method = "bray")

set.seed(42)
permanova_res <- adonis2(dist_bc ~ status, data = as(sample_data(ps_rel), "data.frame"), permutations = 999)
cat("   Test PERMANOVA (Bray-Curtis) : R2 =", round(permanova_res$R2[1], 3),
    ", p-value =", format.pval(permanova_res$`Pr(>F)`[1], digits = 3), "\n")

ord <- ordinate(ps_rel, method = "PCoA", distance = dist_bc)

p_beta <- plot_ordination(ps_rel, ord, color = "status", shape = "status") +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(type = "t", linetype = 2, level = 0.95) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Diversité Bêta (PCoA - Bray-Curtis)",
       subtitle = paste("PERMANOVA: R2 =", round(permanova_res$R2[1], 3),
                        ", p =", format.pval(permanova_res$`Pr(>F)`[1], digits = 3))) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 12), legend.title = element_text(face = "bold"))

# 6. COMPOSITION TAXONOMIQUE ----------------------------------------------------
cat("📊 Génération du diagramme en barres taxonomique...\n")

# Vérification de sécurité avant agrégation
if (ntaxa(ps_filt) == 0) {
  stop("❌ Erreur critique : Aucun taxon conservé après le filtrage.")
}

# Agrégation au niveau du Phylum
ps_phylum <- tax_glom(ps_filt, taxrank = "Phylum")

# Fallback de sécurité si l'agrégation échoue malgré tout
if (ntaxa(ps_phylum) == 0) {
  cat("⚠️ Avertissement : L'agrégation au niveau Phylum a échoué. Utilisation des données non agglomérées.\n")
  ps_phylum <- ps_filt
}

# Transformation en abondance relative et conversion en data.frame long
df_phylum <- psmelt(ps_phylum)

df_phylum <- df_phylum %>%
  group_by(status, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(status) %>%
  mutate(Percentage = Abundance / sum(Abundance) * 100) %>%
  # Regrouper les phyla rares pour la lisibilité du graphique
  mutate(Phylum = ifelse(Percentage < 1 & status == first(status), "Autres (<1%)", Phylum))

# Graphique Barplot empilé
p_taxa <- ggplot(df_phylum, aes(x = status, y = Percentage, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(df_phylum$Phylum)))) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Composition taxonomique relative (Niveau Phylum)",
       x = "Statut du site", y = "Abondance relative (%)") +
  theme_bw() +
  theme(legend.position = "right", legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

# 7. EXPORT DES RÉSULTATS -------------------------------------------------------
cat("💾 Sauvegarde des graphiques en haute résolution...\n")

# Sauvegarde individuelle (Recommandé pour les publications)
ggsave("results/alpha_diversity.png", plot = p_alpha, width = 6, height = 5, dpi = 300)
ggsave("results/beta_diversity.png", plot = p_beta, width = 7, height = 5, dpi = 300)
ggsave("results/taxonomy_barplot.png", plot = p_taxa, width = 8, height = 6, dpi = 300)

# Sauvegarde combinée en PDF
pdf("results/combined_analysis.pdf", width = 12, height = 10)
gridExtra::grid.arrange(p_alpha, p_beta, p_taxa, ncol = 2,
                        layout_matrix = rbind(c(1, 2), c(3, 3)))
dev.off()

cat("🎉 Analyse R terminée avec succès !\n")
cat("📁 Les graphiques sont sauvegardés dans le dossier 'results/'.\n")
