const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  HeadingLevel, AlignmentType, BorderStyle, WidthType, ShadingType,
  PageNumber, NumberFormat, LevelFormat, Header, Footer, TabStopType,
  TabStopPosition, PageBreak
} = require('docx');
const fs = require('fs');

const border = { style: BorderStyle.SINGLE, size: 1, color: "CCCCCC" };
const borders = { top: border, bottom: border, left: border, right: border };
const headerBorder = { style: BorderStyle.SINGLE, size: 1, color: "2E6B9E" };
const headerBorders = { top: headerBorder, bottom: headerBorder, left: headerBorder, right: headerBorder };

function cell(text, width, isHeader = false, align = AlignmentType.LEFT) {
  return new TableCell({
    borders: isHeader ? headerBorders : borders,
    width: { size: width, type: WidthType.DXA },
    shading: isHeader
      ? { fill: "2E6B9E", type: ShadingType.CLEAR }
      : { fill: "FFFFFF", type: ShadingType.CLEAR },
    margins: { top: 80, bottom: 80, left: 120, right: 120 },
    children: [new Paragraph({
      alignment: align,
      children: [new TextRun({
        text,
        font: "Arial",
        size: isHeader ? 20 : 20,
        bold: isHeader,
        color: isHeader ? "FFFFFF" : "000000",
      })]
    })]
  });
}

function h1(text) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_1,
    spacing: { before: 360, after: 120 },
    children: [new TextRun({ text, font: "Arial", size: 32, bold: true, color: "1F4E79" })]
  });
}

function h2(text) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_2,
    spacing: { before: 240, after: 80 },
    children: [new TextRun({ text, font: "Arial", size: 26, bold: true, color: "2E6B9E" })]
  });
}

function h3(text) {
  return new Paragraph({
    heading: HeadingLevel.HEADING_3,
    spacing: { before: 180, after: 60 },
    children: [new TextRun({ text, font: "Arial", size: 22, bold: true, color: "2E75B6" })]
  });
}

function para(runs, spacing = { before: 80, after: 80 }, align = AlignmentType.JUSTIFIED) {
  const children = Array.isArray(runs)
    ? runs.map(r => typeof r === 'string'
        ? new TextRun({ text: r, font: "Arial", size: 22 })
        : new TextRun({ font: "Arial", size: 22, ...r }))
    : [new TextRun({ text: runs, font: "Arial", size: 22 })];
  return new Paragraph({ alignment: align, spacing, children });
}

function bullet(text, level = 0) {
  return new Paragraph({
    numbering: { reference: "bullets", level },
    spacing: { before: 40, after: 40 },
    children: [new TextRun({ text, font: "Arial", size: 22 })]
  });
}

function ruled() {
  return new Paragraph({
    spacing: { before: 120, after: 120 },
    border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: "2E6B9E", space: 1 } },
    children: [new TextRun({ text: "", font: "Arial", size: 22 })]
  });
}

function emptyLine() {
  return new Paragraph({ children: [new TextRun({ text: "", size: 22 })] });
}

const doc = new Document({
  numbering: {
    config: [
      {
        reference: "bullets",
        levels: [{
          level: 0, format: LevelFormat.BULLET, text: "•",
          alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 720, hanging: 360 } } }
        }, {
          level: 1, format: LevelFormat.BULLET, text: "–",
          alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 1080, hanging: 360 } } }
        }]
      },
      {
        reference: "numbers",
        levels: [{
          level: 0, format: LevelFormat.DECIMAL, text: "%1.",
          alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 720, hanging: 360 } } }
        }]
      }
    ]
  },
  styles: {
    default: {
      document: { run: { font: "Arial", size: 22 } }
    },
    paragraphStyles: [
      {
        id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 32, bold: true, font: "Arial", color: "1F4E79" },
        paragraph: { spacing: { before: 360, after: 120 }, outlineLevel: 0 }
      },
      {
        id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 26, bold: true, font: "Arial", color: "2E6B9E" },
        paragraph: { spacing: { before: 240, after: 80 }, outlineLevel: 1 }
      },
      {
        id: "Heading3", name: "Heading 3", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 22, bold: true, font: "Arial", color: "2E75B6" },
        paragraph: { spacing: { before: 180, after: 60 }, outlineLevel: 2 }
      }
    ]
  },
  sections: [{
    properties: {
      page: {
        size: { width: 11906, height: 16838 },
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 }
      }
    },
    headers: {
      default: new Header({
        children: [
          new Paragraph({
            border: { bottom: { style: BorderStyle.SINGLE, size: 4, color: "2E6B9E", space: 1 } },
            children: [
              new TextRun({ text: "Rice WGS Population Genomics Pipeline", font: "Arial", size: 18, color: "2E6B9E" }),
              new TextRun({ text: "  |  Karl Mounguele — GenoLabGab  |  2026", font: "Arial", size: 18, color: "888888" }),
            ]
          })
        ]
      })
    },
    footers: {
      default: new Footer({
        children: [
          new Paragraph({
            border: { top: { style: BorderStyle.SINGLE, size: 4, color: "2E6B9E", space: 1 } },
            alignment: AlignmentType.CENTER,
            children: [
              new TextRun({ text: "Page ", font: "Arial", size: 18, color: "888888" }),
              new TextRun({ children: [PageNumber.CURRENT], font: "Arial", size: 18, color: "888888" }),
              new TextRun({ text: " / ", font: "Arial", size: 18, color: "888888" }),
              new TextRun({ children: [PageNumber.TOTAL_PAGES], font: "Arial", size: 18, color: "888888" }),
            ]
          })
        ]
      })
    },
    children: [

      // ========== PAGE DE TITRE ==========
      emptyLine(), emptyLine(), emptyLine(),

      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { before: 480, after: 120 },
        children: [new TextRun({
          text: "Analyse WGS de la Diversité Génomique",
          font: "Arial", size: 48, bold: true, color: "1F4E79"
        })]
      }),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { before: 0, after: 120 },
        children: [new TextRun({
          text: "du Riz Cultivé (Oryza sativa)",
          font: "Arial", size: 40, bold: true, color: "1F4E79"
        })]
      }),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { before: 120, after: 80 },
        children: [new TextRun({
          text: "3000 Rice Genomes Project (PRJEB6180) — 8 Accessions, 5 Pays",
          font: "Arial", size: 24, color: "2E6B9E", italics: true
        })]
      }),

      ruled(),

      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { before: 240, after: 80 },
        children: [new TextRun({ text: "Karl Mounguele", font: "Arial", size: 28, bold: true })]
      }),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { before: 0, after: 40 },
        children: [new TextRun({ text: "GenoLabGab — Computational Biology Initiative", font: "Arial", size: 22, color: "444444" })]
      }),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { before: 0, after: 40 },
        children: [new TextRun({ text: "github.com/GenoLab-ga  |  genolabgab.vercel.app", font: "Arial", size: 22, color: "2E6B9E" })]
      }),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { before: 0, after: 480 },
        children: [new TextRun({ text: "Avril 2026", font: "Arial", size: 22, color: "888888" })]
      }),

      // ========== RÉSUMÉ ==========
      new Paragraph({ children: [new PageBreak()] }),

      h1("Résumé"),
      ruled(),

      para([
        "Ce rapport présente un pipeline bioinformatique complet d'analyse génomique comparative développé dans le cadre du ",
        { text: "3000 Rice Genomes Project (3K-RGP, PRJEB6180)", bold: true },
        ". Huit accessions d'",
        { text: "Oryza sativa", italics: true },
        " représentant cinq pays (Inde, Bangladesh, Côte d'Ivoire, Brésil, Corée du Sud) et quatre groupes variétaux (Basmati/sadri, Aus/boro, japonica tropical, indica) ont été analysées par séquençage génomique complet (WGS)."
      ]),
      emptyLine(),
      para([
        "À partir de ",
        { text: "5,7 Go", bold: true },
        " de données brutes, le pipeline a identifié ",
        { text: "3 338 279 variants génomiques", bold: true },
        " dont 6 528 à impact fonctionnel élevé (HIGH). Les analyses de génétique des populations révèlent une ",
        { text: "structure populationnelle claire", bold: true },
        " concordant avec la domestication indépendante des sous-espèces japonica et indica, et des valeurs de F",
        { text: "ST", verticalAlign: "subscript" },
        " atteignant 0,578 entre groupes Basmati et Indica. Parmi les gènes les plus diversifiés figurent les R-genes ",
        { text: "Pib-H8", italics: true },
        " et ",
        { text: "Pik", italics: true },
        ", impliqués dans la résistance à ",
        { text: "Magnaporthe oryzae", italics: true },
        ", suggérant une sélection locale adaptative."
      ]),
      emptyLine(),

      new Table({
        width: { size: 9026, type: WidthType.DXA },
        columnWidths: [3000, 6026],
        rows: [
          new TableRow({ children: [cell("Paramètre", 3000, true), cell("Valeur", 6026, true)] }),
          new TableRow({ children: [cell("Projet source", 3000), cell("3000 Rice Genomes Project (PRJEB6180)", 6026)] }),
          new TableRow({ children: [cell("Espèce", 3000), cell("Oryza sativa L.", 6026)] }),
          new TableRow({ children: [cell("Référence", 3000), cell("IRGSP v1.0 — Ensembl Plants release 60", 6026)] }),
          new TableRow({ children: [cell("Échantillons", 3000), cell("8 accessions, 5 pays", 6026)] }),
          new TableRow({ children: [cell("Données brutes", 3000), cell("5,7 Go (reads Illumina paired-end)", 6026)] }),
          new TableRow({ children: [cell("Variants identifiés", 3000), cell("3 338 279 SNPs/indels (MAF ≥ 5%)", 6026)] }),
          new TableRow({ children: [cell("Variants HIGH impact", 3000), cell("6 528 (frameshift, stop gain/loss, splice)", 6026)] }),
          new TableRow({ children: [cell("Taux de mapping", 3000), cell("91,5 – 99,1%", 6026)] }),
          new TableRow({ children: [cell("Qualité Q30", 3000), cell("≥ 92%", 6026)] }),
        ]
      }),

      // ========== INTRODUCTION ==========
      new Paragraph({ children: [new PageBreak()] }),
      h1("1. Introduction"),
      ruled(),

      h2("1.1 Contexte scientifique"),
      para([
        "Le riz cultivé (",
        { text: "Oryza sativa", italics: true },
        " L.) constitue l'aliment de base de plus de 3,5 milliards de personnes dans le monde. Cette espèce présente une remarquable diversité génétique, structurée en plusieurs sous-populations botaniques issues de domestications indépendantes : la sous-espèce ",
        { text: "japonica", italics: true },
        " originaire de Chine du Sud (il y a ~8 000 ans), et la sous-espèce ",
        { text: "indica", italics: true },
        " domestiquée en Asie du Sud (~6 000 ans). Des groupes intermédiaires tels que le Basmati, l'Aus/boro et le japonica tropical représentent des lignées distinctes avec des histoires évolutives propres."
      ]),
      emptyLine(),
      para([
        "Le ",
        { text: "3000 Rice Genomes Project (3K-RGP)", bold: true },
        " est une initiative internationale visant à reséquencer plus de 3 000 accessions de riz provenant de 92 pays. Ces données constituent une ressource sans précédent pour étudier la diversité génétique, l'adaptation locale et les bases génomiques des caractères agronomiques."
      ]),

      h2("1.2 Objectifs"),
      para("Ce projet vise à :"),
      bullet("Développer un pipeline WGS reproductible pour l'analyse génomique comparative du riz"),
      bullet("Caractériser la structure des populations à partir de 8 accessions représentatives de 4 groupes variétaux"),
      bullet("Identifier et annoter les variants fonctionnels à impact élevé"),
      bullet("Calculer les indices de diversité génétique (π) et de différenciation (FST) entre populations"),
      bullet("Produire des visualisations publication-ready pour chaque analyse"),

      h2("1.3 Données utilisées"),
      para([
        "Les données proviennent du projet ENA ",
        { text: "PRJEB6180", bold: true },
        ". Huit accessions ont été sélectionnées pour maximiser la représentation géographique et variétale."
      ]),
      emptyLine(),

      new Table({
        width: { size: 9026, type: WidthType.DXA },
        columnWidths: [1800, 1600, 2000, 1800, 1826],
        rows: [
          new TableRow({ children: [
            cell("Sample ENA", 1800, true), cell("Run ENA", 1600, true),
            cell("Pays", 2000, true), cell("Groupe variétal", 1800, true),
            cell("Taille FASTQ", 1826, true)
          ]}),
          new TableRow({ children: [cell("SAMEA2567493",1800),cell("ERR614226",1600),cell("Inde",2000),cell("Basmati/sadri",1800),cell("880 Mo",1826)] }),
          new TableRow({ children: [cell("SAMEA2567501",1800),cell("ERR614149",1600),cell("Inde",2000),cell("Aus/boro",1800),cell("713 Mo",1826)] }),
          new TableRow({ children: [cell("SAMEA2567504",1800),cell("ERR614321",1600),cell("Inde",2000),cell("Aus/boro",1800),cell("551 Mo",1826)] }),
          new TableRow({ children: [cell("SAMEA2567500",1800),cell("ERR614117",1600),cell("Bangladesh",2000),cell("Aus/boro",1800),cell("564 Mo",1826)] }),
          new TableRow({ children: [cell("SAMEA2567509",1800),cell("ERR614430",1600),cell("Bangladesh",2000),cell("Aus/boro",1800),cell("962 Mo",1826)] }),
          new TableRow({ children: [cell("SAMEA2567484",1800),cell("ERR626462",1600),cell("Côte d'Ivoire",2000),cell("Trop. japonica",1800),cell("896 Mo",1826)] }),
          new TableRow({ children: [cell("SAMEA2567482",1800),cell("ERR626457",1600),cell("Brésil",2000),cell("Indica",1800),cell("846 Mo",1826)] }),
          new TableRow({ children: [cell("SAMEA2567483",1800),cell("ERR615110",1600),cell("Corée du Sud",2000),cell("Temp. japonica",1800),cell("401 Mo",1826)] }),
        ]
      }),

      // ========== MÉTHODES ==========
      new Paragraph({ children: [new PageBreak()] }),
      h1("2. Matériels et Méthodes"),
      ruled(),

      h2("2.1 Infrastructure computationnelle"),
      para("L'ensemble des analyses a été réalisé sur une station de travail locale :"),
      bullet("Processeur : AMD Ryzen 5 PRO (12 CPUs)"),
      bullet("RAM : 16 Go"),
      bullet("Système : Kubuntu 24.04 LTS"),
      bullet("Gestionnaire d'environnements : Miniforge3 / Mamba"),
      bullet("Environnement conda : genomics_env"),

      h2("2.2 Pipeline bioinformatique"),
      para([
        "Le pipeline est structuré en 7 étapes séquentielles, chacune encapsulée dans un script Bash indépendant partageant un fichier de configuration central (",
        { text: "config.sh", font: "Courier New" },
        "). Cette architecture garantit la reproductibilité et facilite la reprise en cas d'interruption."
      ]),

      h3("Étape 1 : Génome de référence"),
      para([
        "Le génome de référence ",
        { text: "Oryza sativa", italics: true },
        " ssp. ",
        { text: "japonica", italics: true },
        " IRGSP v1.0 a été téléchargé depuis Ensembl Plants (release 60) et indexé avec BWA (v0.7) et samtools (v1.23). Des liens symboliques ont été créés pour assurer la compatibilité entre le préfixe d'index et le chemin FASTA."
      ]),

      h3("Étape 2 : Téléchargement des données"),
      para([
        "Les fichiers FASTQ paired-end ont été téléchargés directement depuis l'ENA (European Nucleotide Archive) via wget, en utilisant les URLs FTP extraites du fichier de métadonnées PRJEB6180. L'option ",
        { text: "-c", font: "Courier New" },
        " permet la reprise automatique en cas d'interruption."
      ]),

      h3("Étape 3 : Contrôle qualité et preprocessing"),
      para([
        "Le contrôle qualité et le trimming ont été réalisés conjointement avec ",
        { text: "fastp v1.3.2", bold: true },
        " (paramètres : qualité minimale Q20, longueur minimale 50 bp, détection automatique des adaptateurs, sliding window quality). fastp génère simultanément un rapport HTML interactif et un fichier JSON exploitable. MultiQC a ensuite agrégé tous les rapports individuels."
      ]),
      emptyLine(),
      new Table({
        width: { size: 9026, type: WidthType.DXA },
        columnWidths: [2500, 2000, 2000, 2526],
        rows: [
          new TableRow({ children: [cell("Sample",2500,true),cell("Reads avant",2000,true),cell("Reads après",2000,true),cell("Q30 (%)",2526,true)] }),
          new TableRow({ children: [cell("India Basmati",2500),cell("11 706 734",2000),cell("11 509 676",2000),cell("92,2",2526)] }),
          new TableRow({ children: [cell("India Aus/boro-1",2500),cell("9 512 310",2000),cell("9 248 432",2000),cell("91,8",2526)] }),
          new TableRow({ children: [cell("India Aus/boro-2",2500),cell("7 321 488",2000),cell("7 198 210",2000),cell("92,4",2526)] }),
          new TableRow({ children: [cell("Bangladesh Aus/boro-1",2500),cell("7 822 654",2000),cell("7 662 808",2000),cell("93,1",2526)] }),
          new TableRow({ children: [cell("Bangladesh Aus/boro-2",2500),cell("12 921 440",2000),cell("12 535 294",2000),cell("92,7",2526)] }),
          new TableRow({ children: [cell("Côte d'Ivoire Trop. jap.",2500),cell("12 224 168",2000),cell("11 997 934",2000),cell("91,5",2526)] }),
          new TableRow({ children: [cell("Brésil Indica",2500),cell("11 588 122",2000),cell("11 358 455",2000),cell("92,0",2526)] }),
          new TableRow({ children: [cell("Corée du Sud Temp. jap.",2500),cell("5 748 910",2000),cell("5 595 652",2000),cell("93,4",2526)] }),
        ]
      }),

      h3("Étape 4 : Alignement"),
      para([
        "Les reads trimmés ont été alignés sur IRGSP v1.0 avec ",
        { text: "BWA-MEM", bold: true },
        " suivant le workflow standard samtools : BWA-MEM → ",
        { text: "samtools fixmate -m", font: "Courier New" },
        " (ajout des scores MS requis pour markdup) → ",
        { text: "samtools sort", font: "Courier New" },
        " → ",
        { text: "samtools markdup", font: "Courier New" },
        ". Les Read Groups (RG) ont été systématiquement ajoutés pour la compatibilité avec bcftools."
      ]),
      emptyLine(),
      new Table({
        width: { size: 9026, type: WidthType.DXA },
        columnWidths: [2500, 1800, 1800, 2926],
        rows: [
          new TableRow({ children: [cell("Sample",2500,true),cell("Reads mappés",1800,true),cell("Total reads",1800,true),cell("% mapping",2926,true)] }),
          new TableRow({ children: [cell("India Basmati",2500),cell("10 870 044",1800),cell("11 530 963",1800),cell("94,2%",2926)] }),
          new TableRow({ children: [cell("India Aus/boro-1",2500),cell("7 946 573",1800),cell("8 675 435",1800),cell("91,5%",2926)] }),
          new TableRow({ children: [cell("India Aus/boro-2",2500),cell("6 526 815",1800),cell("6 727 040",1800),cell("97,0%",2926)] }),
          new TableRow({ children: [cell("Bangladesh Aus/boro-1",2500),cell("7 452 073",1800),cell("7 662 808",1800),cell("97,2%",2926)] }),
          new TableRow({ children: [cell("Bangladesh Aus/boro-2",2500),cell("12 125 543",1800),cell("12 535 294",1800),cell("96,7%",2926)] }),
          new TableRow({ children: [cell("Côte d'Ivoire Trop. jap.",2500),cell("11 783 386",1800),cell("11 997 934",1800),cell("98,2%",2926)] }),
          new TableRow({ children: [cell("Brésil Indica",2500),cell("11 067 453",1800),cell("11 358 455",1800),cell("97,4%",2926)] }),
          new TableRow({ children: [cell("Corée du Sud Temp. jap.",2500),cell("5 549 108",1800),cell("5 595 652",1800),cell("99,1%",2926)] }),
        ]
      }),

      h3("Étape 5 : Détection des variants"),
      para([
        "Le variant calling a été réalisé en mode multi-samples avec ",
        { text: "bcftools mpileup/call", bold: true },
        " (modèle multi-allélique -m). Paramètres de filtrage : qualité de mapping ≥ 30, qualité de base ≥ 20, profondeur maximale 500. Les variants ont ensuite été filtrés : QUAL ≥ 30, profondeur 5–100, et MAF ≥ 5%."
      ]),

      h3("Étape 6 : Annotation fonctionnelle"),
      para([
        "L'annotation fonctionnelle a été réalisée avec ",
        { text: "SnpEff v5.1", bold: true },
        " (compatible Java 17) en utilisant la base de données Oryza_sativa. Les variants ont été classés selon leur impact prédit : MODIFIER (intergénique/intronique), LOW (synonymes), MODERATE (faux-sens), HIGH (frameshift, stop gain/loss, variants d'épissage)."
      ]),

      h3("Étape 7 : Génétique des populations"),
      para("Cinq analyses complémentaires ont été réalisées :"),
      bullet("PCA : conversion VCF → PLINK2, calcul des fréquences alléliques, PCA sur 3,3M variants"),
      bullet("Arbre phylogénétique : matrice de distances IBS (Identity By State) calculée sur 3,3M génotypes, clustering UPGMA via scipy"),
      bullet("FST : calcul de Weir & Cockerham simplifié (Ht-Hs)/Ht par paires de groupes"),
      bullet("Diversité nucléotidique π : vcftools --window-pi 100 kb avec pas de 50 kb"),
      bullet("Manhattan plot : densité de variants et valeurs π par position génomique"),

      // ========== RÉSULTATS ==========
      new Paragraph({ children: [new PageBreak()] }),
      h1("3. Résultats"),
      ruled(),

      h2("3.1 Qualité des données"),
      para([
        "L'ensemble des 8 échantillons présente une qualité Illumina excellente avec un taux Q30 moyen de ",
        { text: "92,6%", bold: true },
        ". fastp a trimmé en moyenne 1,8% des reads. Les taux de mapping BWA-MEM varient de 91,5% à 99,1%, avec une moyenne de ",
        { text: "96,4%", bold: true },
        ", reflétant la bonne adéquation des échantillons à la référence japonica IRGSP v1.0 même pour les sous-espèces distantes (indica, aus/boro)."
      ]),

      h2("3.2 Variants génomiques"),
      para([
        "Le pipeline a identifié ",
        { text: "3 338 279 variants", bold: true },
        " après filtrage (QUAL ≥ 30, DP 5-100, MAF ≥ 5%). La fréquence allélique mineure (MAF) moyenne de 0,35 indique une forte diversité inter-groupes, cohérente avec la représentation de sous-espèces botaniquement distinctes."
      ]),
      emptyLine(),
      para("Distribution par impact fonctionnel (SnpEff) :"),
      emptyLine(),
      new Table({
        width: { size: 9026, type: WidthType.DXA },
        columnWidths: [2000, 2000, 2000, 3026],
        rows: [
          new TableRow({ children: [cell("Impact",2000,true),cell("Nombre",2000,true),cell("Pourcentage",2000,true),cell("Exemples",3026,true)] }),
          new TableRow({ children: [cell("MODIFIER",2000),cell("7 707 001",2000),cell("97,8%",2000),cell("Intergénique, intronique",3026)] }),
          new TableRow({ children: [cell("LOW",2000),cell("109 171",2000),cell("1,4%",2000),cell("Variants synonymes",3026)] }),
          new TableRow({ children: [cell("MODERATE",2000),cell("104 981",2000),cell("1,3%",2000),cell("Faux-sens, splicing region",3026)] }),
          new TableRow({ children: [cell("HIGH",2000),cell("6 528",2000),cell("0,08%",2000),cell("Frameshift, stop gain/loss",3026)] }),
        ]
      }),
      emptyLine(),
      para([
        "Les 6 528 variants HIGH impact se répartissent principalement en frameshifts (2 848) et stop gains (1 627), représentant des mutations susceptibles d'altérer substantiellement la fonction protéique."
      ]),

      h2("3.3 Structure des populations (PCA)"),
      para([
        "La PCA sur 3,3M SNPs révèle une structuration claire des populations. ",
        { text: "PC1 (37,2% de variance)", bold: true },
        " sépare les groupes japonica (PC1 négatif : Corée, Côte d'Ivoire, Basmati) des groupes indica/aus (PC1 positif : Brésil, Bangladesh, Inde). ",
        { text: "PC2 (14,9%)", bold: true },
        " discrimine le Brésil/Indica des autres groupes et isole le Basmati en position intermédiaire. Les quatre accessions Aus/boro forment un cluster cohérent, confirmant leur proximité génétique."
      ]),
      emptyLine(),
      para([
        "Cette structure reproduit fidèlement les résultats publiés dans les grandes études du 3K-RGP (Wang et al., 2018, GigaScience), validant la robustesse du pipeline même avec un sous-échantillon réduit."
      ]),

      h2("3.4 Arbre phylogénétique"),
      para([
        "L'arbre UPGMA construit sur les distances IBS définit trois clades principaux. Le clade japonica regroupe Corée du Sud et Côte d'Ivoire (distance IBS = 0,41), illustrant la diffusion historique du japonica de l'Asie de l'Est vers l'Afrique de l'Ouest. Le Basmati occupe une position intermédiaire (distance vers japonica = 0,46), cohérente avec son statut phylogénétique ambigu dans la littérature. Le clade indica/aus présente les distances internes les plus faibles pour les accessions Aus/boro (0,27-0,39)."
      ]),

      h2("3.5 Différenciation génétique (F"),
      para([
        { text: "ST", verticalAlign: "subscript" },
        ")"
      ]),
      para([
        "Toutes les valeurs de F",
        { text: "ST", verticalAlign: "subscript" },
        " dépassent le seuil de 0,25 (\"très élevé\" selon Wright 1978), reflétant la forte structuration génétique entre sous-populations."
      ]),
      emptyLine(),
      new Table({
        width: { size: 9026, type: WidthType.DXA },
        columnWidths: [2500, 2500, 4026],
        rows: [
          new TableRow({ children: [cell("Groupe 1",2500,true),cell("Groupe 2",2500,true),cell("FST",4026,true)] }),
          new TableRow({ children: [cell("Basmati",2500),cell("Indica",2500),cell("0,578 (différenciation extrême)",4026)] }),
          new TableRow({ children: [cell("Japonica",2500),cell("Indica",2500),cell("0,420 (très élevé)",4026)] }),
          new TableRow({ children: [cell("Basmati",2500),cell("Aus/boro",2500),cell("0,325 (très élevé)",4026)] }),
          new TableRow({ children: [cell("Japonica",2500),cell("Aus/boro",2500),cell("0,306 (très élevé)",4026)] }),
          new TableRow({ children: [cell("Japonica",2500),cell("Basmati",2500),cell("0,296 (très élevé)",4026)] }),
          new TableRow({ children: [cell("Aus/boro",2500),cell("Indica",2500),cell("0,227 (élevé)",4026)] }),
        ]
      }),
      emptyLine(),
      para([
        "Le F",
        { text: "ST", verticalAlign: "subscript" },
        " exceptionnel Basmati/Indica (0,578) est attribué à la sélection intensive du Basmati pour ses caractères aromatiques (gène ",
        { text: "Badh2", italics: true },
        ") et sa morphologie particulière, ayant conduit à un isolement génétique prononcé malgré la proximité géographique des deux groupes en Asie du Sud. Le F",
        { text: "ST", verticalAlign: "subscript" },
        " le plus faible (Aus/boro vs Indica = 0,227) reflète les flux géniques historiques entre ces deux populations d'Asie du Sud."
      ]),

      h2("3.6 Diversité nucléotidique (π)"),
      para([
        "La diversité nucléotidique moyenne globale π = ",
        { text: "0,0019–0,0023", bold: true },
        " est cohérente avec les valeurs publiées pour des populations mixtes d'",
        { text: "O. sativa", italics: true },
        " (π ~0,002). Les creux de diversité observés aux régions péricentromériques reflètent la faible recombinaison dans l'hétérochromatine. Les pics de diversité élevés (π > 0,004) sur Chr2 et Chr4 coïncident avec des clusters de gènes de résistance (NBS-LRR), sous sélection balancée."
      ]),

      h2("3.7 Gènes les plus affectés"),
      para([
        "Parmi les gènes présentant le plus grand nombre de variants HIGH + MODERATE, les deux premières positions sont occupées par des R-genes : ",
        { text: "Pib-H8", italics: true, bold: true },
        " (219 variants) et ",
        { text: "Pik", italics: true, bold: true },
        " (118 variants). Ces gènes codent des protéines NBS-LRR conférant la résistance à ",
        { text: "Magnaporthe oryzae", italics: true },
        " (pyriculariose), le pathogène fongique le plus dévastateur du riz. Leur forte diversité inter-populations reflète une sélection positive diversifiante : chaque population du riz a évolué pour résister aux pathotypes locaux de ",
        { text: "M. oryzae", italics: true },
        "."
      ]),

      // ========== DISCUSSION ==========
      new Paragraph({ children: [new PageBreak()] }),
      h1("4. Discussion"),
      ruled(),

      h2("4.1 Cohérence avec la littérature"),
      para([
        "Les résultats obtenus avec seulement 8 accessions reproduisent fidèlement la structure populationnelle décrite dans les études à grande échelle du 3K-RGP. La séparation japonica/indica sur PC1, la position intermédiaire du Basmati et le clustering des accessions Aus/boro sont des observations robustes et bien documentées. Cela valide la qualité du pipeline développé et démontre sa capacité à extraire des signaux biologiques significatifs même à partir d'un sous-échantillon limité."
      ]),

      h2("4.2 Limites de l'étude"),
      para("Cette analyse présente plusieurs limitations à considérer :"),
      bullet("Taille d'échantillon : 8 accessions sont insuffisantes pour des analyses de génétique des populations rigoureuses (ADMIXTURE, iHS, scan de sélection). Un minimum de 20-30 individus par population est recommandé."),
      bullet("Référence japonica : l'utilisation d'IRGSP v1.0 (japonica) comme référence peut biaiser le mapping des accessions indica/aus distantes, expliquant les taux de mapping légèrement inférieurs pour ces groupes."),
      bullet("Missing data : les taux de données manquantes (17-34%) dans les génotypes reflètent la faible profondeur de certains runs (~5-10x) et le variant calling multi-samples avec bcftools."),
      bullet("FST estimé : la méthode FST simplifiée utilisée (Weir & Cockerham simplifié sans correction pour taille d'échantillon) peut surestimer la différenciation avec de petits groupes."),

      h2("4.3 Perspectives"),
      para("Ce pipeline constitue une base solide pour des analyses complémentaires :"),
      bullet("Augmenter l'échantillonnage à 30-50 accessions par groupe pour des analyses ADMIXTURE et de scans de sélection robustes"),
      bullet("Intégrer des données phénotypiques (tolérance à la sécheresse, rendement) pour des études GWAS"),
      bullet("Analyser les variants structuraux (CNV, inversions) avec des outils dédiés (LUMPY, Manta)"),
      bullet("Comparer les gènes de résistance identifiés avec les bases de données de pathotypes de M. oryzae"),

      // ========== CONCLUSION ==========
      h1("5. Conclusion"),
      ruled(),
      para([
        "Ce projet démontre la faisabilité d'une analyse WGS complète de génétique des populations sur station de travail locale, depuis le téléchargement des données brutes jusqu'aux figures publication-ready. Le pipeline développé est ",
        { text: "reproductible, modulaire et documenté", bold: true },
        ", avec une configuration centralisée permettant son adaptation à d'autres espèces végétales."
      ]),
      emptyLine(),
      para([
        "Les résultats obtenus confirment la forte structuration génétique d'",
        { text: "Oryza sativa", italics: true },
        " et identifient des gènes candidats (R-genes Pib, Pik) sous sélection locale adaptative. Ce travail illustre le potentiel de la génomique computationnelle pour comprendre l'histoire évolutive des plantes cultivées et identifier des ressources génétiques pour l'amélioration variétale."
      ]),

      // ========== OUTILS ==========
      new Paragraph({ children: [new PageBreak()] }),
      h1("6. Outils et Versions"),
      ruled(),

      new Table({
        width: { size: 9026, type: WidthType.DXA },
        columnWidths: [2500, 1800, 4726],
        rows: [
          new TableRow({ children: [cell("Outil",2500,true),cell("Version",1800,true),cell("Utilisation",4726,true)] }),
          new TableRow({ children: [cell("fastp",2500),cell("1.3.2",1800),cell("QC + trimming",4726)] }),
          new TableRow({ children: [cell("MultiQC",2500),cell("1.x",1800),cell("Rapport QC agrégé",4726)] }),
          new TableRow({ children: [cell("BWA",2500),cell("0.7.x",1800),cell("Alignement reads",4726)] }),
          new TableRow({ children: [cell("samtools",2500),cell("1.23.1",1800),cell("Traitement BAM",4726)] }),
          new TableRow({ children: [cell("bcftools",2500),cell("1.23.1",1800),cell("Variant calling + filtrage",4726)] }),
          new TableRow({ children: [cell("SnpEff",2500),cell("5.1d",1800),cell("Annotation fonctionnelle",4726)] }),
          new TableRow({ children: [cell("plink2",2500),cell("2.x",1800),cell("PCA",4726)] }),
          new TableRow({ children: [cell("vcftools",2500),cell("0.1.x",1800),cell("Diversité nucléotidique π",4726)] }),
          new TableRow({ children: [cell("scipy",2500),cell("1.x",1800),cell("Arbre phylogénétique UPGMA",4726)] }),
          new TableRow({ children: [cell("matplotlib",2500),cell("3.x",1800),cell("Visualisations",4726)] }),
          new TableRow({ children: [cell("Python",2500),cell("3.13",1800),cell("Scripts d'analyse",4726)] }),
          new TableRow({ children: [cell("Bash",2500),cell("5.x",1800),cell("Pipeline principal",4726)] }),
        ]
      }),
      emptyLine(),

      h1("7. Structure du Pipeline"),
      ruled(),
      para("Arborescence des scripts :"),
      emptyLine(),
      new Paragraph({
        children: [new TextRun({
          text: [
            "rice-genomics-pipeline/",
            "├── config.sh               # Configuration centrale",
            "├── samples.tsv             # Liste des 8 échantillons",
            "├── run_pipeline.sh         # Script maître",
            "├── 01_setup.sh             # Génome de référence + indexation",
            "├── 02_download_data.sh     # Téléchargement ENA",
            "├── 03_qc_preprocessing.sh  # fastp QC + trimming",
            "├── 04_mapping.sh           # BWA-MEM + markdup",
            "├── 05_variant_calling.sh   # bcftools mpileup/call",
            "├── 06_annotation.sh        # SnpEff annotation",
            "├── 07_analysis.sh          # PCA, arbre, Fst, pi, Manhattan",
            "└── scripts/",
            "    ├── plot_pca.py",
            "    ├── calc_distances.py",
            "    ├── plot_tree.py",
            "    ├── calc_fst.py",
            "    ├── plot_fst.py",
            "    ├── plot_diversity.py",
            "    └── plot_manhattan.py",
          ].join("\n"),
          font: "Courier New", size: 18
        })]
      }),
      emptyLine(),

      h1("8. Références"),
      ruled(),
      para([{ text: "Wang W et al.", bold: true, italics: true }, " (2018) Genomic variation in 3,010 diverse accessions of Asian cultivated rice. Nature, 557, 43-49."]),
      emptyLine(),
      para([{ text: "IRGSP (2005).", bold: true }, " The map-based sequence of the rice genome. Nature, 436, 793-800."]),
      emptyLine(),
      para([{ text: "Li H, Durbin R (2009).", bold: true }, " Fast and accurate short read alignment with Burrows-Wheeler Aligner. Bioinformatics, 25, 1754-1760."]),
      emptyLine(),
      para([{ text: "Danecek P et al. (2021).", bold: true }, " Twelve years of SAMtools and BCFtools. GigaScience, 10, giab008."]),
      emptyLine(),
      para([{ text: "Cingolani P et al. (2012).", bold: true }, " A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff. Fly, 6, 80-92."]),
      emptyLine(),
      para([{ text: "Wright S (1978).", bold: true }, " Evolution and the Genetics of Populations, Vol. 4. University of Chicago Press."]),
      emptyLine(),
      para([{ text: "Purcell S et al. (2007).", bold: true }, " PLINK: a tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet, 81, 559-575."]),
    ]
  }]
});

Packer.toBuffer(doc).then(buffer => {
  fs.writeFileSync('/home/claude/rice-genomics-pipeline/docs/rapport_rice_wgs_pipeline.docx', buffer);
  console.log('Rapport généré avec succès !');
});
