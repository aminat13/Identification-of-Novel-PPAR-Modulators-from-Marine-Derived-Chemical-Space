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
