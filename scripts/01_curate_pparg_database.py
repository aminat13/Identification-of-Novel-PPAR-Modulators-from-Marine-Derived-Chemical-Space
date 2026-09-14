# 01_curate_pparg_database.py

import re
import pandas as pd
from rdkit import Chem
from chembl_webresource_client.new_client import new_client


# ---------------------------------------------------------
# File paths
# ---------------------------------------------------------

ligand_file = "data/raw/rcsb_pdb_ligands.xlsx"
citation_file = "data/raw/rcsb_pdb_primary_citations.xlsx"
affinity_file = "data/raw/rcsb_pdb_binding_affinity.xlsx"

output_file = "data/processed/PPARG_Agonists_Initial_Curation.xlsx"


# ---------------------------------------------------------
# Load PDB ligand data
# ---------------------------------------------------------

pparg_df = pd.read_excel(ligand_file)

print("Original rows:", len(pparg_df))
print("Original unique ligands:", pparg_df["Ligand ID"].nunique())


# ---------------------------------------------------------
# Remove crystallisation artefacts and small molecules
# ---------------------------------------------------------

artefacts = [
    "HOH", "DMS", "EOH", "EDO", "GOL", "MPD",
    "PEG", "PG4", "PGE", "1PE", "2PE",
    "ACT", "ACN", "IPA", "MOH",
    "MES", "HEP", "TRS", "CIT",
    "SO4", "PO4", "NO3", "SCN",
    "NA", "K", "CA", "MG", "MN", "ZN", "CL",
    "ACE", "FMT", "ACY",
    "NAG", "BMA", "MAN", "FUC",
    "DTT", "BME", "IMD", "TLA",
    "UNX"
]

pparg_df = pparg_df[
    (pparg_df["Ligand MW"] >= 150) &
    (~pparg_df["Ligand ID"].isin(artefacts))
].copy()

print("Rows after MW/artefact filtering:", len(pparg_df))
print("Unique ligands remaining:", pparg_df["Ligand ID"].nunique())


# ---------------------------------------------------------
# Add primary citation information
# ---------------------------------------------------------

citation_df = pd.read_excel(citation_file)

pparg_df = pparg_df.merge(
    citation_df,
    on="Entry ID",
    how="left"
)


# ---------------------------------------------------------
# Remove entries from clearly non-agonist studies
# ---------------------------------------------------------

non_agonist_terms = [
    "antagonist",
    "inverse agonist",
    "inverse-agonist",
    "inhibitor",
    "inhibition",
    "repression",
    "repressive",
    "corepressor",
    "corepression",
    "phosphorylation inhibition",
    "covalent inhibitor",
    "ppi inhibitors",
    "competitive ppar gamma ppi inhibitors",
    "ligand-resistance",
    "resistance syndrome",
    "ammonia"
]

pattern = "|".join(re.escape(term) for term in non_agonist_terms)

non_agonist_mask = pparg_df["Title"].astype(str).str.contains(
    pattern,
    case=False,
    na=False
)

removed_non_agonists = pparg_df[non_agonist_mask].copy()
pparg_df = pparg_df[~non_agonist_mask].copy()

print("Rows removed by title screening:", len(removed_non_agonists))
print("Rows remaining:", len(pparg_df))


# ---------------------------------------------------------
# Add EC50 values deposited in the PDB
# ---------------------------------------------------------

affinity_df = pd.read_excel(affinity_file)

pparg_df.columns = pparg_df.columns.str.strip()
affinity_df.columns = affinity_df.columns.str.strip()

affinity_df = affinity_df.rename(columns={"Ligand": "Ligand ID"})

ec50_df = affinity_df[
    affinity_df["Type"].astype(str).str.upper().str.strip() == "EC50"
].copy()

ec50_df["Value"] = pd.to_numeric(
    ec50_df["Value"],
    errors="coerce"
)

# Remove repeated copies of the same measurement
ec50_df = ec50_df.drop_duplicates(
    subset=["Entry ID", "Ligand ID", "Value"]
)

# Average multiple EC50 values reported for the same PDB ligand
ec50_df = (
    ec50_df
    .groupby(["Entry ID", "Ligand ID"], as_index=False)
    .agg(EC50=("Value", "mean"))
)

pparg_df = pparg_df.merge(
    ec50_df,
    on=["Entry ID", "Ligand ID"],
    how="left"
)

print(
    "Rows with PDB EC50:",
    pparg_df["EC50"].notna().sum()
)


# ---------------------------------------------------------
# Search ChEMBL for ligands still missing EC50 data
# ---------------------------------------------------------

ligands_with_ec50 = set(
    pparg_df.loc[
        pparg_df["EC50"].notna(),
        "Ligand ID"
    ]
)

missing_ec50 = (
    pparg_df[
        ~pparg_df["Ligand ID"].isin(ligands_with_ec50)
    ][
        ["Entry ID", "Ligand ID", "Ligand Name", "Ligand SMILES"]
    ]
    .drop_duplicates()
    .copy()
)


def get_inchikey(smiles):
    mol = Chem.MolFromSmiles(str(smiles))

    if mol is None:
        return None

    return Chem.MolToInchiKey(mol)


missing_ec50["InChIKey"] = missing_ec50["Ligand SMILES"].apply(
    get_inchikey
)

missing_ec50 = missing_ec50[
    missing_ec50["InChIKey"].notna()
].copy()


# Match PDB ligands to ChEMBL compounds
molecule_client = new_client.molecule

inchikeys = missing_ec50["InChIKey"].unique().tolist()

chembl_matches = list(
    molecule_client.filter(
        molecule_structures__standard_inchi_key__in=inchikeys
    ).only(
        "molecule_chembl_id",
        "molecule_structures"
    )
)

chembl_ids = []

for compound in chembl_matches:

    structure = compound.get("molecule_structures") or {}

    chembl_ids.append({
        "InChIKey": structure.get("standard_inchi_key"),
        "ChEMBL_ID": compound.get("molecule_chembl_id")
    })

chembl_ids = pd.DataFrame(chembl_ids).dropna()

missing_ec50 = missing_ec50.merge(
    chembl_ids,
    on="InChIKey",
    how="left"
)


# ---------------------------------------------------------
# Retrieve PPARγ EC50 activity from ChEMBL
# ---------------------------------------------------------

activity_client = new_client.activity

compound_ids = (
    missing_ec50["ChEMBL_ID"]
    .dropna()
    .unique()
    .tolist()
)

activities = list(
    activity_client.filter(
        molecule_chembl_id__in=compound_ids,
        target_chembl_id="CHEMBL235",
        standard_type="EC50"
    ).only(
        "molecule_chembl_id",
        "standard_value",
        "standard_units",
        "standard_relation"
    )
)

chembl_ec50 = pd.DataFrame(activities)

if not chembl_ec50.empty:

    chembl_ec50 = chembl_ec50.rename(columns={
        "molecule_chembl_id": "ChEMBL_ID",
        "standard_value": "ChEMBL_EC50",
        "standard_units": "ChEMBL_Units",
        "standard_relation": "ChEMBL_Relation"
    })

    chembl_results = missing_ec50.merge(
        chembl_ec50,
        on="ChEMBL_ID",
        how="inner"
    )

    chembl_results["ChEMBL_EC50"] = pd.to_numeric(
        chembl_results["ChEMBL_EC50"],
        errors="coerce"
    )

    # Use EC50 values reported in nM
    chembl_results = chembl_results[
        chembl_results["ChEMBL_Units"] == "nM"
    ].copy()

    chembl_results = chembl_results.drop_duplicates(
        subset=[
            "Entry ID",
            "Ligand ID",
            "ChEMBL_EC50"
        ]
    )

    chembl_average = (
        chembl_results
        .groupby(
            ["Entry ID", "Ligand ID"],
            as_index=False
        )
        .agg(
            ChEMBL_EC50=("ChEMBL_EC50", "mean")
        )
    )

    pparg_df = pparg_df.merge(
        chembl_average,
        on=["Entry ID", "Ligand ID"],
        how="left"
    )

    # Keep PDB EC50 where available and use ChEMBL to fill gaps
    pparg_df["EC50"] = pparg_df["EC50"].fillna(
        pparg_df["ChEMBL_EC50"]
    )


# ---------------------------------------------------------
# Remove redundant rows without EC50
# ---------------------------------------------------------

ligands_with_ec50 = set(
    pparg_df.loc[
        pparg_df["EC50"].notna(),
        "Ligand ID"
    ]
)

pparg_df = pparg_df[
    ~(
        pparg_df["Ligand ID"].isin(ligands_with_ec50) &
        pparg_df["EC50"].isna()
    )
].copy()


# ---------------------------------------------------------
# Save for final manual curation
# ---------------------------------------------------------

# Remaining EC50 values from BindingDB / primary literature and
# WT/mutant classifications were added during manual database review.

with pd.ExcelWriter(output_file) as writer:

    pparg_df.to_excel(
        writer,
        sheet_name="PPARG_candidates",
        index=False
    )

    removed_non_agonists.to_excel(
        writer,
        sheet_name="removed_non_agonists",
        index=False
    )


print("\nInitial PPARγ curation complete")
print("Rows retained:", len(pparg_df))
print("Unique ligands:", pparg_df["Ligand ID"].nunique())
print("Rows with EC50:", pparg_df["EC50"].notna().sum())
