# Identification of Novel PPARγ Modulators from Marine-Derived Chemical Space

## Project Overview

This repository contains the computational workflow developed for my MSc research project in **Technologies and Analytics in Precision Medicine at RCSI**.

The project investigated whether metabolites derived from brown seaweed could represent potential novel modulators of **peroxisome proliferator-activated receptor gamma (PPARγ)**, a nuclear receptor involved in lipid metabolism, glucose homeostasis and inflammatory signalling.

PPARγ is highly expressed in colonic epithelial cells and represents a potential therapeutic target in inflammatory bowel disease (IBD) and colorectal cancer (CRC). However, long-term use of synthetic PPARγ agonists can be limited by adverse effects, motivating interest in alternative chemical scaffolds.

This study therefore explored metabolites from:

- *Saccharina latissima*
- *Laminaria digitata*
- *Undaria pinnatifida*

using an integrated computational drug-discovery workflow combining cheminformatics, structural analysis, ligand-based screening and molecular docking.

The complete MSc thesis is included in this repository for full scientific background, methodology, results and discussion.

---

## Workflow Summary

```text
PPARγ reference ligand curation + seaweed metabolite curation
        ↓
RDKit molecular descriptors and Morgan fingerprints
        ↓
Butina clustering + chemotype/scaffold classification
        ↓
PPARγ binding-pocket and ligand interaction analysis
        ↓
Structure–activity relationship analysis
        ↓
Consensus pharmacophore development
        ↓
Tanimoto similarity + pharmacophore + AI-assisted screening
        ↓
Molecular docking and integrated candidate prioritisation
```
The workflow combined both ligand-based and structure-based approaches so that candidate compounds could be assessed from multiple independent perspectives rather than relying on a single screening metric.

---

## Datasets

### PPARγ Agonist Reference Database

Human PPARγ ligand-bound structures were curated from the **RCSB Protein Data Bank** using UniProt accession **P37231**.

Structures were filtered to remove crystallographic artefacts, small molecules below 150 Da and irrelevant entries. Biological activity information was incorporated from sources including ChEMBL, BindingDB and primary literature.

The final reference dataset contained:

- **226 co-crystallised PPARγ agonist complexes**
- **218 wild-type complexes**
- **8 mutant complexes**
- **66 complexes with experimentally reported EC50 values**

### Seaweed Metabolite Database

Natural products associated with the three brown seaweed species were obtained from the **LOTUS Natural Products database**.

Compounds lacking usable SMILES structures or below 150 Da were excluded.

The final database contained **51 metabolites**:

- 31 *Undaria pinnatifida*
- 16 *Saccharina latissima*
- 4 *Laminaria digitata*

---

## Key Analyses

The project included:

- Molecular descriptor calculation using RDKit
- Morgan fingerprint generation and Tanimoto similarity
- Butina clustering and chemical-space analysis
- UMAP visualisation
- SMARTS-based chemotype and scaffold classification
- PPARγ binding-site interaction analysis
- Structure–activity relationship analysis using EC50/pEC50
- Pharmacophore development and virtual screening
- AI-assisted descriptor-based candidate prioritisation
- Molecular docking in MOE
- Integration of multiple screening outputs for final candidate ranking

---

## Key Findings

- The PPARγ reference database contained diverse chemotypes including non-TZD synthetic compounds, TZDs, natural-product-derived ligands and fatty-acid/lipid-like agonists.
- *Undaria pinnatifida* showed the greatest scaffold diversity among the seaweed species.
- Conserved ligand-recognition residues included **CYS285, ARG288, SER289, HIS323, MET364 and HIS449**.
- Fatty-acid metabolites were consistently prioritised across structural similarity, pharmacophore screening, AI-assisted assessment and molecular docking.
- One seaweed metabolite showed an exact fingerprint match to the known endogenous PPARγ ligand arachidonic acid.
- Different screening approaches showed only partial agreement, supporting the value of combining complementary methods.
- All three compounds prioritised consistently across Tanimoto similarity, pharmacophore screening and AI-assisted screening also ranked among the ten highest-scoring docked compounds.

These findings support further investigation of selected seaweed-derived metabolites, while not constituting evidence of biological agonist activity without experimental validation.

---

## Key Skills Demonstrated

### Cheminformatics

- RDKit molecular descriptor calculation
- SMILES handling and molecular structure processing
- Morgan fingerprints
- Tanimoto similarity
- Butina clustering
- Chemical-space analysis
- SMARTS substructure matching
- Chemotype and scaffold classification

### Structural Bioinformatics

- Protein–ligand interaction analysis
- Binding-pocket characterisation
- Residue-level interaction mapping
- Structural interpretation of ligand recognition
- Pharmacophore modelling
- Molecular docking
- Docking validation and pose comparison

### Data Analysis and Programming

- Python
- pandas
- NumPy
- SciPy
- Matplotlib
- Seaborn
- Data cleaning and integration
- Statistical analysis
- Data visualisation
- Multi-source dataset curation

### Drug Discovery and Scientific Analysis

- Structure–activity relationship analysis
- Ligand-based virtual screening
- Structure-based virtual screening
- Candidate prioritisation
- Integration of complementary screening methods
- Critical interpretation of computational predictions
- Scientific literature integration
- Research workflow design

---

## Technologies and Resources

### Programming

- Python
- pandas
- NumPy
- SciPy
- Matplotlib
- Seaborn

### Cheminformatics

- RDKit
- Morgan fingerprints
- Tanimoto similarity
- Butina clustering
- SMARTS

### Structural and Virtual Screening Tools

- Molecular Operating Environment (MOE)
- Ligand interaction analysis
- Molecular surface mapping
- Pharmacophore modelling
- Molecular docking

### Databases

- RCSB Protein Data Bank
- UniProt
- ChEMBL
- BindingDB
- LOTUS Natural Products database

---

## Repository Structure

```text
.
├── README.md
│
├── thesis/
│   └── Aminat_Coster_MSc_Thesis.pdf
│
├── scripts/
│   ├── 01_curate_pparg_database.py
│   ├── 02_curate_seaweed_database.py
│   ├── 03_calculate_molecular_descriptors.py
│   ├── 04_chemotype_clustering.py
│   ├── 05_binding_interaction_analysis.py
│   ├── 06_structure_activity_analysis.py
│   ├── 07_tanimoto_screening.py
│   └── 08_screening_integration.py
│
├── data/
│   ├── raw/
│   │   ├── PPARγ ligand and affinity exports
│   │   ├── LOTUS seaweed metabolite data
│   │   └── MOE-exported screening and docking files
│   │
│   └── processed/
│       ├── curated PPARγ agonist database
│       ├── curated seaweed metabolite database
│       ├── molecular descriptor datasets
│       ├── clustered and classified datasets
│       ├── ligand interaction data
│       └── screening result files
│
└── results/
    ├── figures/
    │   ├── PPARγ agonist chemotype distribution
    │   ├── seaweed metabolite scaffold distribution
    │   ├── ligand interaction frequency analysis
    │   ├── structure–activity relationship figures
    │   ├── Tanimoto similarity screening figures
    │   ├── AI-assisted prioritisation results
    │   └── integrated virtual-screening analysis
    │
    └── tables/
        ├── residue interaction frequencies
        ├── PPARγ activity classifications
        ├── SAR functional-group analysis
        ├── Tanimoto screening results
        ├── pharmacophore screening hits
        └── integrated screening and docking results
```

The `scripts/` directory follows the order of the computational workflow presented in the thesis, progressing from database curation and molecular characterisation through structural analysis, virtual screening and final candidate prioritisation.

The `results/` directory contains the final figures and tables used in the submitted MSc thesis rather than every exploratory output generated during development.

Some structural analyses, including binding-pocket surface mapping, pharmacophore generation, pharmacophore overlays, docking pose visualisation and structural superposition, were performed directly in **Molecular Operating Environment (MOE)**. These final outputs are included in the repository where relevant but are not reproduced programmatically in Python.

The scripts have been reorganised from the original research notebooks into cleaner, concise versions while preserving the analytical workflow and methods used in the project.

---

## Results Included in the Repository

The repository includes the final outputs used in the thesis, including:

- PPARγ agonist chemotype distribution
- Seaweed metabolite scaffold distribution
- Conserved ligand-binding residue analysis
- Protein–ligand interaction summaries
- Structure–activity relationship figures
- Tanimoto similarity screening results
- Pharmacophore screening outputs
- AI-assisted prioritisation results
- Integrated virtual-screening figures
- Molecular docking rankings

Figures generated directly in MOE, such as molecular surface representations, pharmacophore visualisations and structural overlays, are retained as final outputs but are not reproduced programmatically in Python.

---

## Reproducibility Notes

- Python scripts used for data curation, cheminformatics analysis and visualisation are provided in the `scripts/` directory.
- Analyses performed directly in MOE are documented in the thesis and corresponding outputs are included in the repository.
- Some stages of database curation involved manual literature review and annotation and are therefore not fully automated.
- The project was designed for computational candidate prioritisation rather than prediction of confirmed biological activity.
- Experimental validation would be required to confirm PPARγ binding and agonist activity.

---

## MSc Thesis

**The Identification of Novel PPARγ Modulators from a Marine-Derived Chemical Space**  
MSc Technologies and Analytics in Precision Medicine  
Royal College of Surgeons in Ireland (RCSI), 2026

The full thesis is included in this repository and provides detailed background, methodology, results, discussion and references.

---

## Author

**Aminat Coster**

MSc Technologies and Analytics in Precision Medicine  
BSc (Hons) Biomolecular Sciences – Biotechnology and Drug Development
