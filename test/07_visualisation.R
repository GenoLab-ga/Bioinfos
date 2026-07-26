# =============================================================
# Étape 7b : Visualisation population genetics
# Projet   : PRJNA1465284 – P. falciparum drug resistance
# =============================================================

library(ggplot2)
library(dplyr)
library(reshape2)

OUTDIR <- "results/population_genetics/figures"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ── Figure 1 : PCA ───────────────────────────────────────────
pca <- read.table(
  "results/population_genetics/pca/pca_results.eigenvec",
  header = FALSE
)
colnames(pca)[1:4] <- c("FID", "IID", "PC1", "PC2")

eigenval <- read.table(
  "results/population_genetics/pca/pca_results.eigenval"
)
var_pc1 <- round(eigenval[1,1] / sum(eigenval[,1]) * 100, 1)
var_pc2 <- round(eigenval[2,1] / sum(eigenval[,1]) * 100, 1)

p_pca <- ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.7, size = 2.5, color = "#1D9E75") +
  labs(
    title = "PCA – Structure génétique des isolats",
    subtitle = "PRJNA1465284 · 318 isolats · P. falciparum Éthiopie",
    x = paste0("PC1 (", var_pc1, "% variance)"),
    y = paste0("PC2 (", var_pc2, "% variance)")
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUTDIR, "01_pca.pdf"),
       p_pca, width = 8, height = 6)
ggsave(file.path(OUTDIR, "01_pca.png"),
       p_pca, width = 8, height = 6, dpi = 300)

cat("Figure 1 PCA sauvegardée.\n")

# ── Figure 2 : Fréquences mutations clés ─────────────────────
mutations <- data.frame(
  Gene     = c("dhfr", "dhfr", "dhfr", "crt",
               "mdr1", "dhps", "kelch13"),
  Mutation = c("N51I", "C59R", "S108N", "K76T",
               "Y184F", "A581G", "P441L"),
  Drug     = c("SP", "SP", "SP", "CQ",
               "LUM", "SP", "ACT"),
  Freq     = c(0.965, 0.597, 0.994, 0.810,
               0.009, 0.009, 0.042)
)

drug_colors <- c(
  "SP"  = "#E24B4A",
  "CQ"  = "#BA7517",
  "LUM" = "#7F77DD",
  "ACT" = "#1D9E75"
)

p_freq <- ggplot(
  mutations,
  aes(x = reorder(Mutation, -Freq),
      y = Freq * 100,
      fill = Drug)
) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 25, linetype = "dashed",
             color = "gray50", linewidth = 0.5) +
  scale_fill_manual(values = drug_colors) +
  scale_y_continuous(limits = c(0, 105),
                     labels = function(x) paste0(x, "%")) +
  geom_text(
    aes(label = paste0(round(Freq * 100, 1), "%")),
    vjust = -0.5, size = 3.5
  ) +
  facet_grid(. ~ Gene, scales = "free_x", space = "free") +
  labs(
    title    = "Fréquences des mutations de résistance",
    subtitle = "318 isolats cliniques · sites de résurgence éthiopiens",
    x        = "Mutation",
    y        = "Fréquence allélique (%)",
    fill     = "Antipaludéen"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title   = element_text(face = "bold"),
    strip.text   = element_text(face = "bold", size = 12),
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(OUTDIR, "02_resistance_frequencies.pdf"),
       p_freq, width = 10, height = 6)
ggsave(file.path(OUTDIR, "02_resistance_frequencies.png"),
       p_freq, width = 10, height = 6, dpi = 300)

cat("Figure 2 Fréquences sauvegardée.\n")

# ── Figure 3 : Hétérozygosité par échantillon ────────────────
het_file <- "results/population_genetics/frequencies/heterozygosity.het"
if (file.exists(het_file)) {
  het <- read.table(het_file, header = TRUE)
  het$F_stat <- het$F
  
  p_het <- ggplot(het, aes(x = F_stat)) +
    geom_histogram(bins = 40, fill = "#378ADD",
                   color = "white", alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "#A32D2D") +
    labs(
      title    = "Distribution du coefficient d'inbreeding (F)",
      subtitle = "F < 0 = excès d'hétérozygotes · F > 0 = excès d'homozygotes",
      x        = "Coefficient F",
      y        = "Nombre d'isolats"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(OUTDIR, "03_inbreeding.pdf"),
         p_het, width = 8, height = 5)
  ggsave(file.path(OUTDIR, "03_inbreeding.png"),
         p_het, width = 8, height = 5, dpi = 300)
  
  cat("Figure 3 Hétérozygosité sauvegardée.\n")
}

# ── Figure 4 : Relatedness distribution ──────────────────────
ibd_file <- "results/population_genetics/relatedness/ibd_matrix.genome"
if (file.exists(ibd_file)) {
  ibd <- read.table(ibd_file, header = TRUE)
  
  p_ibd <- ggplot(ibd, aes(x = PI_HAT)) +
    geom_histogram(bins = 50, fill = "#534AB7",
                   color = "white", alpha = 0.85) +
    geom_vline(xintercept = 0.25, linetype = "dashed",
               color = "#A32D2D", linewidth = 0.7) +
    annotate("text", x = 0.27, y = Inf,
             label = "Clones potentiels\n(PI_HAT > 0.25)",
             vjust = 1.5, hjust = 0, size = 3.5,
             color = "#A32D2D") +
    labs(
      title    = "Distribution du relatedness (IBD)",
      subtitle = "PI_HAT = proportion allèles identiques par descendance",
      x        = "PI_HAT",
      y        = "Nombre de paires d'isolats"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  
  N_CLONES <- sum(ibd$PI_HAT > 0.25)
  cat(sprintf("Paires hautement apparentées (PI_HAT > 0.25) : %d\n",
              N_CLONES))
  
  ggsave(file.path(OUTDIR, "04_relatedness.pdf"),
         p_ibd, width = 8, height = 5)
  ggsave(file.path(OUTDIR, "04_relatedness.png"),
         p_ibd, width = 8, height = 5, dpi = 300)
  
  cat("Figure 4 Relatedness sauvegardée.\n")
}

cat("\nToutes les figures générées dans :", OUTDIR, "\n")