# Identification-of-Novel-PPAR-Modulators-from-Marine-Derived-Chemical-Space

## Overview

This repository contains the computational workflow developed for my MSc research project in **Technologies and Analytics in Precision Medicine at RCSI**.

The project investigated whether metabolites derived from brown seaweed could represent potential novel modulators of **peroxisome proliferator-activated receptor gamma (PPARγ)**, a nuclear receptor involved in lipid metabolism, glucose homeostasis and inflammatory signalling.

PPARγ is highly expressed in colonic epithelial cells and represents a potential therapeutic target in inflammatory bowel disease (IBD) and colorectal cancer (CRC). Although synthetic PPARγ agonists have demonstrated anti-inflammatory and therapeutic effects, their long-term use can be limited by adverse effects. This project therefore explored marine-derived natural products as an alternative chemical space for identifying structurally diverse candidate PPARγ modulators.

The study combined **cheminformatics, structural bioinformatics and virtual screening approaches** to:

- curate a reference database of experimentally characterised PPARγ agonists,
- characterise the chemical diversity of known agonists and brown seaweed metabolites,
- investigate the structural features associated with PPARγ ligand recognition,
- develop a consensus PPARγ pharmacophore,
- screen seaweed metabolites using complementary ligand- and structure-based approaches,
- and prioritise candidate compounds for further experimental investigation.

The full MSc thesis is included in this repository for detailed scientific context, methodology, interpretation and discussion.

---

## Research Question

**Can structurally diverse metabolites from brown seaweed reproduce the molecular characteristics associated with known PPARγ agonists and therefore represent plausible candidates for further investigation as PPARγ modulators?**

The project focused on metabolites reported in three brown seaweed species:

- *Saccharina latissima*
- *Laminaria digitata*
- *Undaria pinnatifida*

Rather than relying on a single screening method, candidate compounds were evaluated using several complementary approaches to capture different aspects of ligand recognition.

---

## Project Workflow

The computational workflow consisted of seven main stages:

### 1. Database Curation

Two datasets were constructed for the project.

#### PPARγ reference agonist database

Experimentally characterised PPARγ ligand-bound structures were collected from the **RCSB Protein Data Bank (PDB)** using the human PPARγ UniProt accession **P37231**.

Structures were filtered to retain high-quality human X-ray crystal structures and remove:

- crystallographic artefacts,
- very small ligands below 150 Da,
- clearly irrelevant non-agonist structures,
- and duplicate records.

Available biological activity information, particularly **EC50 values**, was integrated from sources including:

- RCSB PDB,
- ChEMBL,
- BindingDB,
- and primary literature.

Ligands were additionally classified according to wild-type or mutant receptor status.

The final reference database contained **226 co-crystallised PPARγ agonist complexes**, including **218 wild-type and 8 mutant complexes**, with experimentally reported EC50 values available for 66 records.

#### Brown seaweed metabolite database

Natural products associated with *S. latissima*, *L. digitata* and *U. pinnatifida* were retrieved from the **LOTUS Natural Products database**.

Compounds were filtered to remove:

- records without a valid SMILES representation,
- compounds below 150 Da,
- and duplicate structures.

The final curated seaweed database contained **51 metabolites**:

- 31 from *Undaria pinnatifida*
- 16 from *Saccharina latissima*
- 4 from *Laminaria digitata*

---

### 2. Molecular Descriptor Calculation and Chemical-Space Characterisation

RDKit was used to calculate molecular descriptors for both the PPARγ agonist and seaweed metabolite databases.

Descriptors included:

- Molecular weight
- LogP
- Topological polar surface area (TPSA)
- Hydrogen-bond donors
- Hydrogen-bond acceptors
- Rotatable bonds
- Aromatic rings
- Heavy atoms

Morgan fingerprints were then generated from ligand SMILES representations using:

- radius = 2
- 2048 bits
- chirality enabled

Chemical similarity was calculated using the **Tanimoto coefficient**.

Ligands were grouped using **Butina clustering**, with dataset-specific distance thresholds selected to reflect differences in structural diversity between the reference agonist and seaweed datasets.

Cluster structure was evaluated using:

- two-dimensional UMAP visualisation,
- silhouette score,
- intra-cluster similarity,
- and manual inspection of representative compounds.

Substructure matching using SMARTS patterns was then used to classify compounds into broader chemotype or scaffold groups.

---

### 3. PPARγ Binding-Pocket Characterisation

Representative PPARγ–ligand complexes spanning the major agonist chemotypes were examined using the **Molecular Operating Environment (MOE)**.

Protein–ligand interactions were characterised using the MOE ligand interaction tools, including:

- hydrogen-bond donor interactions,
- hydrogen-bond acceptor interactions,
- ionic interactions,
- π-H interactions,
- and ligand proximity to binding-site residues.

Interaction data were consolidated and analysed in Python to identify residues repeatedly involved in ligand recognition across structurally different agonist classes.

Frequently observed residues included:

- CYS285
- ARG288
- SER289
- HIS323
- MET364
- HIS449

Surface mapping was also used to examine the electrostatic and lipophilic characteristics of the PPARγ ligand-binding pocket.

This analysis helped establish the structural environment that candidate ligands would need to complement.

---

### 4. Structure–Activity Relationship Analysis

Experimentally reported EC50 values were converted to **pEC50** values and used to divide reference agonists into higher- and lower-activity groups.

Structural features were then compared between these groups.

Functional groups assessed included:

- carboxylic acids,
- hydroxyl groups,
- amides,
- sulfonamides,
- ethers,
- ketones,
- aromatic and heteroaromatic rings,
- halogens,
- trifluoromethyl groups,
- flexible linkers,
- and double bonds.

This analysis was used alongside the binding-pocket interaction data to identify recurring molecular features associated with PPARγ activity.

---

### 5. Consensus Pharmacophore Development

Information from three complementary analyses was integrated:

1. binding-pocket surface characterisation,
2. conserved protein–ligand interactions,
3. structure–activity relationship analysis.

A consensus pharmacophore was then generated in MOE using a representative PPARγ crystal structure.

The final model represented the recurring structural requirements observed across experimentally characterised agonists and contained aromatic, hydrogen-bond accepting, hydrophobic and negative ionisable features together with volume exclusion regions representing sterically inaccessible parts of the binding pocket.

This pharmacophore was subsequently used to screen the brown seaweed metabolite database.

---

### 6. Virtual Screening and Candidate Prioritisation

Seaweed metabolites were evaluated using three independent screening strategies.

#### Tanimoto similarity screening

Each seaweed metabolite was compared against representative reference agonists using Morgan fingerprints and Tanimoto similarity.

Similarity results were evaluated across four broad PPARγ agonist classes:

- Non-TZD synthetic
- TZD synthetic
- Natural product-derived
- Fatty acid / lipid-like

#### Pharmacophore screening

Seaweed metabolites were conformationally prepared in MOE and screened against the consensus PPARγ pharmacophore.

Pharmacophore RMSD was used to assess geometric agreement between each metabolite and the consensus feature arrangement.

#### AI-assisted prioritisation

Candidate structures were also assessed using a descriptor-based AI-assisted prioritisation approach.

For each compound, information including:

- SMILES structure,
- molecular weight,
- LogP,
- TPSA,
- hydrogen-bond donor and acceptor counts,
- and rotatable bond count

was used to assign a **High, Medium or Low priority** together with a structural rationale.

These approaches were treated as complementary rather than interchangeable, allowing candidates supported by multiple independent methods to be identified with greater confidence.

---

### 7. Molecular Docking and Integrated Candidate Evaluation

Molecular docking was performed in MOE using the PPARγ crystal structure **5Z6S**.

The docking site was defined using residues within 4.5 Å of the co-crystallised ligand.

The docking protocol was first validated by re-docking the crystallographic reference ligand and assessing:

- docking score,
- recovery of the experimental binding orientation,
- pose RMSD,
- and conservation of known protein–ligand interactions.

Following validation, all 51 seaweed metabolites were docked using the same receptor preparation and docking parameters.

Docking results were then compared with:

- Tanimoto similarity,
- pharmacophore RMSD,
- and AI-assisted prioritisation.

This final integration step was used to identify compounds consistently supported across multiple computational approaches.

---

## Key Findings

Several broad findings emerged from the integrated workflow.

### Chemical diversity

The reference PPARγ agonist database contained a broad range of chemotypes, with non-TZD synthetic ligands forming the largest group alongside TZDs, natural product-derived agonists and fatty-acid/lipid-like ligands.

Among the seaweed species, ***Undaria pinnatifida*** displayed the greatest scaffold diversity.

### Conserved ligand recognition

Despite substantial structural variation among reference agonists, several PPARγ binding-site residues repeatedly contributed to ligand recognition.

This suggested that conserved molecular interactions may be more important than overall scaffold similarity alone.

### Seaweed-derived fatty acids emerged consistently

Fatty-acid metabolites from *Undaria pinnatifida* were repeatedly prioritised across:

- structural similarity screening,
- pharmacophore screening,
- AI-assisted prioritisation,
- and molecular docking.

One metabolite showed an exact structural fingerprint match to the known endogenous PPARγ ligand arachidonic acid.

### Different screening methods captured different information

Pharmacophore RMSD and docking score showed only a weak relationship, while some compounds performed strongly in docking despite relatively low structural similarity.

This supported the use of multiple complementary approaches rather than relying on a single screening metric.

### Consensus candidate prioritisation

All three compounds prioritised consistently across Tanimoto similarity, pharmacophore screening and AI-assisted assessment also ranked among the ten highest-scoring compounds during molecular docking.

The results therefore provided computational support for further investigation of several seaweed-derived metabolites, particularly lipid-like candidates, while **not constituting evidence of biological agonist activity without experimental validation**.

---

## Repository Structure

```text
.
├── README.md
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
│   ├── 08_screening_integration.py
│   └── ...
│
├── results/
│   ├── figures/
│   └── tables/
│
└── data/
    ├── raw/
    └── processed/
