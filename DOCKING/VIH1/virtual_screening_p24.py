#!/usr/bin/env python3
"""
Pipeline de criblage virtuel - Protéine p24 (4XFX)
Auteur : généré pour votre projet VIH
Usage  : python3 virtual_screening_p24.py
"""

import os
import subprocess
import glob
import csv
import re
from pathlib import Path

# ============================================================
# CONFIGURATION — Modifiez ces valeurs selon votre setup
# ============================================================

VINA_PATH      = "/usr/bin/vina"               # ou chemin absolu : "/usr/bin/vina"
OBABEL_PATH    = "/usr/bin/obabel"             # Open Babel

# Dossiers
BASE_DIR       = Path("docking_p24")
PROTEIN_DIR    = BASE_DIR / "proteine"
LIGANDS_SDF    = BASE_DIR / "ligands" / "natural_products.sdf"  # votre fichier SDF
LIGANDS_DIR    = BASE_DIR / "ligands" / "pdbqt"
RESULTS_DIR    = BASE_DIR / "resultats"
LOGS_DIR       = BASE_DIR / "logs"

# Fichier récepteur
RECEPTOR       = PROTEIN_DIR / "p24_clean.pdbqt"

# Coordonnées de la grille (blind docking p24)
CENTER_X       =   9.2
CENTER_Y       = -23.2
CENTER_Z       =   0.0
SIZE_X         =  80
SIZE_Y         =  80
SIZE_Z         =  80

# Paramètres Vina
EXHAUSTIVENESS =  8
NUM_MODES      =  9
ENERGY_RANGE   =  3
CPU            =  8    # nombre de cœurs à utiliser

# ============================================================
# ÉTAPE 0 — Créer l'arborescence des dossiers
# ============================================================

def create_directories():
    print("\n[1/4] Création des dossiers...")
    for d in [PROTEIN_DIR, LIGANDS_DIR, RESULTS_DIR, LOGS_DIR]:
        d.mkdir(parents=True, exist_ok=True)
    print("    ✓ Arborescence créée")
    print(f"""
    docking_p24/
    ├── proteine/P24.pdbqt        ← placez 4XFX_receptor.pdbqt ici
    ├── ligands/
    │   ├── natural_products.sdf  ← placez votre SDF ici
    │   └── pdbqt/       ← fichiers convertis (auto)
    ├── resultats/       ← poses de docking
    └── logs/            ← scores et logs
    """)

# ============================================================
# ÉTAPE 1 — Conversion SDF → PDBQT (Open Babel)
# ============================================================

def convert_sdf_to_pdbqt():
    print("\n[2/4] Conversion SDF → PDBQT avec Open Babel...")

    if not LIGANDS_SDF.exists():
        print(f"    ✗ Fichier SDF introuvable : {LIGANDS_SDF}")
        print("      Placez votre fichier SDF dans docking_p24/ligands/")
        return False

    cmd = [
        OBABEL_PATH,
        str(LIGANDS_SDF),
        "-O", str(LIGANDS_DIR / "ligand_.pdbqt"),
        "-m",                        # un fichier par molécule
        "--partialcharge", "gasteiger",
        "-h",                        # ajouter hydrogènes
        #"--gen3d",                   # générer coordonnées 3D si absentes
        "2>/dev/null"
    ]

    print(f"    Commande : {' '.join(cmd[:-1])}")

    try:
        result = subprocess.run(
            cmd[:-1],
            capture_output=True,
            text=True
        )
        ligands = list(LIGANDS_DIR.glob("ligand_*.pdbqt"))
        print(f"    ✓ {len(ligands)} ligands convertis en PDBQT")

        if len(ligands) == 0:
            print("    ✗ Aucun ligand converti — vérifiez votre fichier SDF")
            return False
        return True

    except FileNotFoundError:
        print("    ✗ Open Babel introuvable — vérifiez l'installation")
        print("      sudo apt install openbabel")
        return False

# ============================================================
# ÉTAPE 2 — Docking AutoDock Vina (boucle sur tous les ligands)
# ============================================================

def run_docking():
    print("\n[3/4] Lancement du docking AutoDock Vina...")
    if not RECEPTOR.exists():
        print(f"    ✗ Récepteur introuvable : {RECEPTOR}")
        return []

    ligands = sorted(LIGANDS_DIR.glob("ligand_*.pdbqt"))
    if not ligands:
        print("    ✗ Aucun ligand PDBQT trouvé")
        return []

    print(f"    Récepteur   : {RECEPTOR}")
    print(f"    Ligands     : {len(ligands)} molécules")
    print(f"    Grille      : center=({CENTER_X}, {CENTER_Y}, {CENTER_Z})")
    print(f"                  size=({SIZE_X}, {SIZE_Y}, {SIZE_Z}) Å")
    print(f"    Exhaustiveness : {EXHAUSTIVENESS}")
    print(f"    CPUs        : {CPU}")
    print()

    scores = []
    total = len(ligands)

    for i, ligand_path in enumerate(ligands, 1):
        ligand_name = ligand_path.stem
        output_path = RESULTS_DIR / f"{ligand_name}_out.pdbqt"
        log_path    = LOGS_DIR    / f"{ligand_name}.log"

        print(f"    [{i:3d}/{total}] Docking {ligand_name}...", end=" ", flush=True)

        cmd = [
            VINA_PATH,
            "--receptor",        str(RECEPTOR),
            "--ligand",          str(ligand_path),
            "--center_x",        str(CENTER_X),
            "--center_y",        str(CENTER_Y),
            "--center_z",        str(CENTER_Z),
            "--size_x",          str(SIZE_X),
            "--size_y",          str(SIZE_Y),
            "--size_z",          str(SIZE_Z),
            "--exhaustiveness",  str(EXHAUSTIVENESS),
            "--num_modes",       str(NUM_MODES),
            "--energy_range",    str(ENERGY_RANGE),
            "--cpu",             str(CPU),
            "--out",             str(output_path),
            # ✅ --log supprimé
        ]

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300
            )

            # Sauvegarder stdout comme log
            with open(log_path, "w") as f:
                f.write(result.stdout)

            # Extraire le score depuis stdout
            best_score = extract_best_score_from_stdout(result.stdout)

            if best_score is not None:
                scores.append({
                    "ligand":  ligand_name,
                    "score":   best_score,
                    "log":     str(log_path),
                    "output":  str(output_path)
                })
                print(f"score = {best_score:.2f} kcal/mol")
            else:
                print("score non trouvé")

        except subprocess.TimeoutExpired:
            print("timeout !")
        except FileNotFoundError:
            print("\n    ✗ Vina introuvable")
            break

    return scores
# ============================================================
# ÉTAPE 3 — Extraction du meilleur score depuis le log
# ============================================================

def extract_best_score_from_stdout(stdout_text):
    """
    Extrait le meilleur score depuis stdout Vina 1.2.5
    Format : '   1       -8.002          0          0'
    """
    for line in stdout_text.splitlines():
        # Format avec pipes : "   1 |    -8.002 | ..."
        match = re.match(r"^\s+1\s*\|\s*(-?\d+\.?\d*)\s*\|", line)
        if match:
            return float(match.group(1))
        # Format sans pipes : "   1       -8.002      0.000      0.000"
        match2 = re.match(r"^\s+1\s+(-?\d+\.?\d*)\s+", line)
        if match2:
            return float(match2.group(1))
    return None
# ============================================================
# ÉTAPE 4 — Analyse et tri des résultats
# ============================================================

def analyze_results(scores):
    print("\n[4/4] Analyse des résultats...")

    if not scores:
        print("    ✗ Aucun résultat à analyser")
        return

    # Trier par score croissant (plus négatif = meilleure affinité)
    scores_sorted = sorted(scores, key=lambda x: x["score"])

    # Sauvegarder en CSV
    csv_path = BASE_DIR / "resultats_docking.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["ligand", "score", "output", "log"])
        writer.writeheader()
        writer.writerows(scores_sorted)

    print(f"    ✓ Résultats sauvegardés : {csv_path}")
    print()
    print("=" * 55)
    print("    TOP 10 — Meilleures affinités de liaison")
    print("=" * 55)
    print(f"    {'#':<4} {'Ligand':<25} {'Score (kcal/mol)'}")
    print("-" * 55)

    for rank, entry in enumerate(scores_sorted[:10], 1):
        print(f"    {rank:<4} {entry['ligand']:<25} {entry['score']:.2f}")

    print("=" * 55)
    print()
    print(f"    Total dockés     : {len(scores)}")
    print(f"    Score moyen      : {sum(s['score'] for s in scores)/len(scores):.2f} kcal/mol")
    print(f"    Meilleur score   : {scores_sorted[0]['score']:.2f} kcal/mol ({scores_sorted[0]['ligand']})")
    print(f"    Score seuil hit  : ≤ -7.0 kcal/mol (recommandé)")
    print()

    # Identifier les hits
    hits = [s for s in scores_sorted if s["score"] <= -7.0]
    print(f"    Hits candidats (≤ -7.0 kcal/mol) : {len(hits)} molécules")
    print()
    print("    Prochaines étapes :")
    print("    1. Visualiser les hits dans PyMol ou ChimeraX")
    print("    2. Analyser les interactions dans Discovery Studio")
    print("    3. Filtrer par ADMET avec SwissADME")

# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("=" * 55)
    print("   CRIBLAGE VIRTUEL — Protéine p24 VIH (4XFX)")
    print("   AutoDock Vina — Produits naturels")
    print("=" * 55)

    create_directories()

    # Vérifier que le récepteur existe
    if not RECEPTOR.exists():
        print(f"\n⚠️  Placez votre fichier 4XFX_receptor.pdbqt dans :")
        print(f"    {PROTEIN_DIR}/")
        print("\nEnsuite relancez le script.")
        exit(1)

    # Lancer le pipeline
    ok = convert_sdf_to_pdbqt()
    if ok:
        scores = run_docking()
        analyze_results(scores)
        print("\n✓ Pipeline terminé !")
    else:
        print("\n✗ Correction nécessaire avant de continuer.")
