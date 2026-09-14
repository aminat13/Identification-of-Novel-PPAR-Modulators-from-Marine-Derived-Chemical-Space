# 03_calculate_molecular_descriptors.py

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


# ---------------------------------------------------------
# Calculate RDKit descriptors
# ---------------------------------------------------------

def calculate_descriptors(df, smiles_column):

    molecules = df[smiles_column].apply(
        Chem.MolFromSmiles
    )

    df["Molecular_Weight"] = molecules.apply(
        Descriptors.MolWt
    )

    df["LogP"] = molecules.apply(
        Descriptors.MolLogP
    )

    df["TPSA"] = molecules.apply(
        Descriptors.TPSA
    )

    df["HBD"] = molecules.apply(
        Descriptors.NumHDonors
    )

    df["HBA"] = molecules.apply(
        Descriptors.NumHAcceptors
    )

    df["Rotatable_Bonds"] = molecules.apply(
        Descriptors.NumRotatableBonds
    )

    df["Aromatic_Rings"] = molecules.apply(
        rdMolDescriptors.CalcNumAromaticRings
    )

    df["Heavy_Atoms"] = molecules.apply(
        Descriptors.HeavyAtomCount
    )

    return df


# ---------------------------------------------------------
# PPARγ agonist database
# ---------------------------------------------------------

pparg_df = pd.read_excel(
    "data/processed/PPARG_Agonists_Curated.xlsx"
)

pparg_df = calculate_descriptors(
    pparg_df,
    "Ligand_SMILES"
)

pparg_df.to_excel(
    "data/processed/PPARG_Agonists_Descriptors.xlsx",
    index=False
)

print(
    "PPARγ descriptors calculated:",
    len(pparg_df),
    "records"
)


# ---------------------------------------------------------
# Seaweed metabolite database
# ---------------------------------------------------------

seaweed_df = pd.read_excel(
    "data/processed/Seaweed_Metabolites_Curated.xlsx"
)

seaweed_df = calculate_descriptors(
    seaweed_df,
    "SMILES"
)

seaweed_df.to_excel(
    "data/processed/Seaweed_Metabolites_Descriptors.xlsx",
    index=False
)

print(
    "Seaweed descriptors calculated:",
    len(seaweed_df),
    "records"
)
