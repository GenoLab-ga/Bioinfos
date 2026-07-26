# Guide d'Exécution du Pipeline Single-Cell RNA-seq
## Plasmodium falciparum - PRJEB55754

---

## 📋 Structure du pipeline

```
ÉTAPE 1: FastQC
   ↓ (vérifier qualité reads)
ÉTAPE 2: Indexation Salmon
   ↓ (créer index génome, fait une seule fois)
ÉTAPE 3: Quantification Salmon
   ↓ (quantifier chaque échantillon)
ÉTAPE 4: Conversion Transcripts → Gènes
   ↓ (créer matrice de comptage)
ÉTAPE 5: Seurat QC
   ↓ (filtrer et normaliser)
RÉSULTAT FINAL: Matrice prête pour clustering
```

---

## 🚀 Avant de commencer

### 1. Vérifier l'environnement

```bash
# Activation mamba
mamba activate rnaseq_env

# Vérifier les outils
fastqc --version      # Doit afficher version
salmon --version      # Doit afficher version
R --version          # Doit afficher version
```

### 2. Copier les scripts

```bash
# Les scripts sont dans /mnt/user-data/outputs/
cp /mnt/user-data/outputs/*.sh ~/

# Les rendre exécutables
chmod +x ~/*.sh

# Vérifier
ls -la ~/*fastqc*.sh ~/*salmon*.sh ~/*seurat*.sh ~/*MASTER*.sh
```

### 3. Configuration project

```bash
# Les données doivent être dans:
# ~/Documents/github_projet/Plasmodium/data/
# ├── ERR11471971_1.fastq.gz
# ├── ERR11471971_2.fastq.gz
# └── ...

# Vérifier:
ls -lh ~/Documents/github_projet/Plasmodium/data/*.fastq.gz
```

---

## 💻 Exécution

### Option 1 : Lancer tout d'un coup

```bash
cd ~/Documents/github_projet/Plasmodium

# Lancer TOUT (recommandé pour première fois)
bash ~/MASTER_PIPELINE.sh all
```

**Durée estimée :** 45-60 minutes (selon CPU/disque)

### Option 2 : Lancer étape par étape

```bash
# Étape 1 seulement
bash ~/01_fastqc_style.sh

# Étape 2
bash ~/02_salmon_indexing.sh

# Étape 3
bash ~/03_salmon_quantification.sh

# Étape 4
python ~/04_salmon_to_genecounts.sh

# Étape 5
Rscript ~/05_seurat_qc.sh
```

### Option 3 : Utiliser le script master avec numéro

```bash
# Jusqu'à l'étape 2 (FastQC + Index)
bash ~/MASTER_PIPELINE.sh 2

# Jusqu'à l'étape 3 (+ Quantification)
bash ~/MASTER_PIPELINE.sh 3

# Tout
bash ~/MASTER_PIPELINE.sh all

# Vérifier le status
bash ~/MASTER_PIPELINE.sh status
```

---

## 📊 Vérifier les résultats

Après chaque étape, vérifiez les outputs:

### Après ÉTAPE 1 (FastQC)
```bash
firefox ~/Documents/github_projet/Plasmodium/qc_reports/multiqc_report.html &

# Ou voir les stats:
tail -20 ~/Documents/github_projet/Plasmodium/logs/analysis.log
```

### Après ÉTAPE 3 (Salmon Quant)
```bash
# Voir les fichiers quantifiés
ls -lh ~/Documents/github_projet/Plasmodium/quantification/salmon_output/*/quant.sf

# Vérifier mapping rates:
tail -5 ~/Documents/github_projet/Plasmodium/logs/salmon_ERR*.log
```

### Après ÉTAPE 4 (Matrice)
```bash
# Voir la matrice
head -5 ~/Documents/github_projet/Plasmodium/quantification/matrices/combined_counts_matrix.csv

# Dimensions
wc -l ~/Documents/github_projet/Plasmodium/quantification/matrices/combined_counts_matrix.csv
```

### Après ÉTAPE 5 (Seurat)
```bash
# Voir les plots QC
ls -lh ~/Documents/github_projet/Plasmodium/analysis/01_qc/*.pdf

# Voir l'objet Seurat
ls -lh ~/Documents/github_projet/Plasmodium/analysis/01_qc/seurat_qc_filtered.rds
```

---

## ⚠️ Troubleshooting

### Erreur: "salmon command not found"
```bash
# Solution:
mamba activate rnaseq_env
which salmon  # Doit afficher chemin
```

### Erreur: "Fichiers FASTQ manquants"
```bash
# Vérifier:
ls -lh ~/Documents/github_projet/Plasmodium/data/*.fastq.gz

# Doit afficher 8 fichiers (4 paires):
# ERR11471971_1.fastq.gz
# ERR11471971_2.fastq.gz
# ERR11471975_1.fastq.gz
# ...
```

### Erreur: "Index non trouvé"
```bash
# Vérifier que l'étape 2 a complété:
ls -lh ~/Documents/github_projet/Plasmodium/reference/Pfalciparum3D7_salmon_index/

# Doit voir: sa/, txpInfo.json, etc.
```

### Erreur: "R package not found"
```R
# Depuis terminal R:
Rscript -
> install.packages("Seurat")
> install.packages("ggplot2")
# ou via conda:
mamba install r-seurat r-ggplot2
```

---

## 📝 Structure de sortie

Après completion complète, vous aurez:

```
~/Documents/github_projet/Plasmodium/
├── data/
│   ├── ERR11471971_1.fastq.gz  (données brutes)
│   └── ...
├── reference/
│   ├── Pfalciparum3D7_transcripts.fasta
│   └── Pfalciparum3D7_salmon_index/
├── qc_reports/
│   ├── fastqc/
│   │   └── *_fastqc.html
│   └── multiqc_report.html
├── quantification/
│   ├── salmon_output/
│   │   ├── ERR11471971/quant.sf
│   │   └── ...
│   └── matrices/
│       └── combined_counts_matrix.csv  ← MAIN RESULT
├── analysis/
│   └── 01_qc/
│       ├── seurat_qc_filtered.rds  ← R OBJECT
│       ├── normalized_counts.csv
│       ├── 01_qc_pre_filtering.pdf
│       └── 02_qc_post_filtering.pdf
└── logs/
    ├── analysis.log
    └── salmon_*.log
```

---

## 🎓 Pédagogie : Comprendre chaque étape

| Étape | Ce qu'elle fait | Sortie |
|-------|-----------------|--------|
| **1** | Vérifier qualité reads brutes | Rapports HTML |
| **2** | Créer index pour aligner | Dossier index |
| **3** | Quantifier expression par échantillon | quant.sf files |
| **4** | Agréger transcripts en gènes | Matrice CSV |
| **5** | Nettoyer et normaliser | Objet Seurat R |

---

## ✅ Checklist Finale

- [ ] Environnement mamba activé
- [ ] Scripts copiés et rendus exécutables
- [ ] Données FASTQ présentes (8 fichiers)
- [ ] Exécution sans erreurs
- [ ] Rapport MultiQC consulté
- [ ] Matrice CSV créée
- [ ] Objet Seurat créé

---

## 🎯 Prochaines étapes après le pipeline

Une fois le pipeline complété:

1. **Clustering**: PCA → UMAP → Louvain
2. **Identification cellulaires**: Markers genes
3. **Analyse biologique**: Annotation types cellulaires
4. **Comparaison**: Contrôle vs. unknown samples
5. **Publication**: Figures + résultats

---

**Prêt ? Exécutez :**

```bash
bash ~/MASTER_PIPELINE.sh all
```

**Et revenez me dire les résultats !** 🚀

