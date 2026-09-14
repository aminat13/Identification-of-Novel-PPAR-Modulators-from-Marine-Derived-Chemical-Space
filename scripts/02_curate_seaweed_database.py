# 02_curate_seaweed_database.py

import requests
import pandas as pd
from urllib.parse import quote


# ---------------------------------------------------------
# Seaweed species included in the study
# ---------------------------------------------------------

species = [
    "Saccharina latissima",
    "Laminaria digitata",
    "Undaria pinnatifida"
]

all_compounds = []


# ---------------------------------------------------------
# Search LOTUS for each species
# ---------------------------------------------------------

for seaweed in species:

    print(f"Searching LOTUS for {seaweed}...")

    query = quote(seaweed)

    url = (
        "https://lotus.naturalproducts.net/api/search/simple"
        f"?query={query}"
    )

    response = requests.get(url)
    data = response.json()

    compounds = pd.DataFrame(
        data["naturalProducts"]
    )

    compounds = compounds[
        [
            "lotus_id",
            "traditional_name",
            "smiles",
            "molecular_formula",
            "molecular_weight"
        ]
    ].copy()

    compounds["Species"] = seaweed

    all_compounds.append(compounds)


# ---------------------------------------------------------
# Combine species into one database
# ---------------------------------------------------------

seaweed_df = pd.concat(
    all_compounds,
    ignore_index=True
)

seaweed_df = seaweed_df.rename(columns={
    "lotus_id": "LOTUS_ID",
    "traditional_name": "Compound_Name",
    "smiles": "SMILES",
    "molecular_formula": "Molecular_Formula",
    "molecular_weight": "Molecular_Weight"
})


print("\nCompounds retrieved:", len(seaweed_df))


# ---------------------------------------------------------
# Remove compounds without usable structures
# ---------------------------------------------------------

seaweed_df["Molecular_Weight"] = pd.to_numeric(
    seaweed_df["Molecular_Weight"],
    errors="coerce"
)

seaweed_df = seaweed_df[
    seaweed_df["SMILES"].notna()
].copy()

seaweed_df = seaweed_df[
    seaweed_df["SMILES"].str.strip() != ""
].copy()


# ---------------------------------------------------------
# Remove compounds below 150 Da
# ---------------------------------------------------------

seaweed_df = seaweed_df[
    seaweed_df["Molecular_Weight"] >= 150
].copy()


# ---------------------------------------------------------
# Remove duplicate compounds
# ---------------------------------------------------------

seaweed_df = seaweed_df.drop_duplicates(
    subset=["LOTUS_ID", "SMILES"],
    keep="first"
)


# ---------------------------------------------------------
# Save curated database
# ---------------------------------------------------------

seaweed_df.to_excel(
    "data/processed/Seaweed_Metabolites_Curated.xlsx",
    index=False
)


print("\nSeaweed database complete")
print("Total metabolites:", len(seaweed_df))

print("\nMetabolites by species:")
print(seaweed_df["Species"].value_counts())
