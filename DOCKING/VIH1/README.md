# Criblage Virtuel de Produits Naturels contre la Capside p24 du VIH-1

**Auteur : Keny Karl MOUNGUELE**  
**Année : 2026**

---

## Résumé

La protéine de capside p24 du VIH-1 constitue une cible thérapeutique prometteuse en raison de son rôle essentiel dans l'assemblage et la maturation virale. Dans cette étude, un criblage virtuel par docking moléculaire a été réalisé sur la structure cristallographique de la p24 (PDB : 4XFX) contre une bibliothèque de 100 produits naturels issus de la base de données ZINC Natural Products.

Le docking a été effectué à l'aide d'AutoDock Vina v1.2.5 (blind docking, exhaustiveness=8) en utilisant une grille de 80×80×80 Å centrée sur les coordonnées (9,2 ; -23,2 ; 0,0). Les résultats ont permis d'identifier 21 molécules présentant un score ≤ -7,0 kcal/mol, parmi lesquelles le meilleur hit est le composé CNP0013909.3 (ZINC000150341587), un triterpénoïde de type oléanane (C₄₄H₅₇N₃O₄, MW = 691,94 g/mol), avec une énergie libre de liaison de -8,99 kcal/mol.

---

## Outils utilisés

| Outil | Version | Usage |
|---|---|---|
| AutoDock Vina | v1.2.5 | Docking moléculaire |
| Open Babel | v3.x | Conversion SDF → PDBQT |
| ADMETlab 2.0 | - | Propriétés physicochimiques et toxicité |
| pkCSM | - | Prédiction absorption/distribution |
| Python | 3.x | Automatisation pipeline batch |
| ZINC Natural Products | - | Bibliothèque de ligands |
| Protein Data Bank (PDB) | - | Structure 3D cible (4XFX) |
| COCONUT | - | Identification des produits naturels |

---

## Résultats clés

- **21 hits** identifiés sur 100 molécules criblées (score ≤ -7,0 kcal/mol)
- **Meilleur hit** : CNP0013909.3 (ZINC000150341587)
- **Classe chimique** : Triterpénoïde oléanane (C₄₄H₅₇N₃O₄)
- **Score de docking** : -8,99 kcal/mol
- **Absorption intestinale (HIA)** : 100%

---

## Structure du repository

---

## Auteur & Contact

**Keny Karl MOUNGUELE**  
Ingénieur Biotechnologiste & Bioinformaticien  
Fondateur de GenoLabGab

🌐 Site web : [https://genolabgab.vercel.app](https://genolabgab.vercel.app)  
💼 LinkedIn : [linkedin.com/in/karl-mounguele](https://www.linkedin.com/in/karl-mounguele)  
📧 Email : genolabgab@proton.me
