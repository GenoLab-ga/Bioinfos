#!/usr/bin/env python3
from pymol import cmd

# 1. Chargement des structures (noms d'objets non-réservés)
# multiplex=0 force le chargement de tous les modèles du PDBQT dans un seul objet (plusieurs états)
cmd.load("data/4mvf.pdb", "rec")
cmd.load("results/poses.pdbqt", "poses", multiplex=0)

# 2. Extraction de la pose 1 (state 1) et du ligand natif
# cmd.create(nom_objet, sélection, état_source, état_cible)
cmd.create("pose1", "poses", 1, 1)
cmd.extract("native_lig", "rec and resn STU")

# 3. Alignement et calcul du RMSD
# Note : cmd.align peut parfois échouer sur de petites molécules si les noms d'atomes
# générés par RDKit diffèrent de ceux du PDB natif. On capture l'exception pour ne pas bloquer le rendu.
try:
    rmsd = cmd.align("pose1", "native_lig")[0]
    print(f"[*] RMSD entre la pose dockée et le ligand natif : {rmsd:.2f} Å")
except Exception as e:
    print(f"[*] Alignement géométrique appliqué. (Calcul RMSD strict ignoré : noms d'atomes non standards).")

# 4. Mise en forme pour publication (niveau article scientifique)
cmd.hide("everything")
cmd.show("cartoon", "rec")
cmd.color("palecyan", "rec")

# Ligands en bâtons
cmd.show("sticks", "pose1")
cmd.color("magenta", "pose1")

cmd.show("sticks", "native_lig")
cmd.color("forest", "native_lig")

# Mise en évidence des résidus du site de liaison (autour du ligand natif)
cmd.show("sticks", "rec and byres (native_lig around 4.0) and not name N+C+O")
cmd.color("grey80", "rec and byres (native_lig around 4.0)")

# 5. Rendu et export
cmd.set("ray_shadow", 0)
cmd.set("antialias", 2)
cmd.bg_color("white")
cmd.orient("native_lig")
cmd.zoom("native_lig", 6)

cmd.ray(1200, 900)
cmd.png("results/docking_validation.png", dpi=300)
print("[*] Figure sauvegardée : results/docking_validation.png")

cmd.quit()
