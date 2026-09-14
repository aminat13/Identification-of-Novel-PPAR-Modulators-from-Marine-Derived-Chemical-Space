# 04_chemotype_clustering.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import umap

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
from matplotlib.patches import Patch


# ---------------------------------------------------------
# File paths
# ---------------------------------------------------------

pparg_file = "data/processed/PPARG_Agonists_Descriptors.xlsx"
seaweed_file = "data/processed/Seaweed_Metabolites_Descriptors.xlsx"

pparg_output = "data/processed/PPARG_Agonists_Clustered.xlsx"
seaweed_output = "data/processed/Seaweed_Metabolites_Clustered.xlsx"


# ---------------------------------------------------------
# Morgan fingerprints
# ---------------------------------------------------------

morgan_generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=2048,
    includeChirality=True
)


def cluster_database(df, smiles_col, id_col, cutoff):
    """Generate Morgan fingerprints and perform Butina clustering."""

    compounds = []
    fingerprints = []

    for _, row in df.iterrows():

        mol = Chem.MolFromSmiles(str(row[smiles_col]))

        if mol is None:
            continue

        compounds.append((mol, row[id_col]))
        fingerprints.append(
            morgan_generator.GetFingerprint(mol)
        )

    # Butina uses a triangular distance matrix
    distance_matrix = []

    for i in range(1, len(fingerprints)):

        similarities = DataStructs.BulkTanimotoSimilarity(
            fingerprints[i],
            fingerprints[:i]
        )

        distance_matrix.extend(
            [1 - similarity for similarity in similarities]
        )

    clusters = Butina.ClusterData(
        distance_matrix,
        len(fingerprints),
        cutoff,
        isDistData=True
    )

    clusters = sorted(
        clusters,
        key=len,
        reverse=True
    )

    # Assign cluster number to each molecule
    cluster_labels = np.full(
        len(compounds),
        -1,
        dtype=int
    )

    for cluster_id, cluster in enumerate(clusters):

        for molecule_index in cluster:
            cluster_labels[molecule_index] = cluster_id

    cluster_lookup = pd.DataFrame({
        id_col: [compound[1] for compound in compounds],
        "Cluster_ID": cluster_labels
    })

    df = df.merge(
        cluster_lookup,
        on=id_col,
        how="left"
    )

    print(
        f"{id_col}: {len(clusters)} clusters from "
        f"{len(compounds)} molecules"
    )

    print(
        "Singleton clusters:",
        sum(1 for cluster in clusters if len(cluster) == 1)
    )

    return df, fingerprints, clusters


# =========================================================
# PPARγ AGONIST DATABASE
# =========================================================

pparg_df = pd.read_excel(pparg_file)
pparg_df.columns = pparg_df.columns.str.strip()

# Keep wild-type structures for clustering
if "Mutation_Status" in pparg_df.columns:
    pparg_df = pparg_df[
        pparg_df["Mutation_Status"] == "WT"
    ].copy()


pparg_df, pparg_fingerprints, pparg_clusters = cluster_database(
    pparg_df,
    smiles_col="Ligand_SMILES",
    id_col="Ligand_ID",
    cutoff=0.60
)


# ---------------------------------------------------------
# PPARγ chemotype classification
# ---------------------------------------------------------

chemotype_smarts = {
    "TZD": "O=C1NC(=O)CS1",
    "Indazole": "c1ccc2c(c1)cnn2",
    "Benzoxazole": "c1ccc2ocnc2c1",
    "Benzimidazole": "c1nc2ccccc2[nH0,nH]1",
    "Indole": "c1ccc2[nH,nH0]ccc2c1",
    "Flavonoid": "O=c1cc(-c2ccccc2)oc2ccccc12",
    "Sulfonamide": "S(=O)(=O)N",
    "Fatty Acid": "CCCC(=O)O",
    "Phenylpropanoic Acid": "c1ccccc1CC(C(=O)O)",
    "Cercosporamide": "O=C1OC2=CC=CC(=O)C2C1=O",
    "Terpenoid": "CC(C)=CCC",
    "Steroid": "C1CC2CCCC3CCCC(C1)C23"
}

priority_order = [
    "TZD",
    "Flavonoid",
    "Cercosporamide",
    "Steroid",
    "Indazole",
    "Benzoxazole",
    "Benzimidazole",
    "Indole",
    "Sulfonamide",
    "Terpenoid",
    "Phenylpropanoic Acid",
    "Fatty Acid"
]

compiled_patterns = {
    name: Chem.MolFromSmarts(smarts)
    for name, smarts in chemotype_smarts.items()
}


def classify_compound(smiles):

    mol = Chem.MolFromSmiles(str(smiles))

    if mol is None:
        return "Invalid SMILES"

    for chemotype in priority_order:

        pattern = compiled_patterns[chemotype]

        if mol.HasSubstructMatch(pattern):
            return chemotype

    return "Other/Unclassified"


pparg_df["Chemotype"] = (
    pparg_df["Ligand_SMILES"]
    .apply(classify_compound)
)


# ---------------------------------------------------------
# Collapse chemotypes into broad classes
# ---------------------------------------------------------

broad_class_map = {
    "TZD": "TZD",

    "Benzimidazole": "Non-TZD synthetic",
    "Indole": "Non-TZD synthetic",
    "Indazole": "Non-TZD synthetic",
    "Benzoxazole": "Non-TZD synthetic",
    "Phenylpropanoic Acid": "Non-TZD synthetic",
    "Sulfonamide": "Non-TZD synthetic",

    "Fatty Acid": "Fatty acids",

    "Terpenoid": "Natural product-derived",
    "Cercosporamide": "Natural product-derived",
    "Flavonoid": "Natural product-derived",
    "Steroid": "Natural product-derived",

    "Other/Unclassified": "Unclassified"
}

pparg_df["Broad_Class"] = (
    pparg_df["Chemotype"]
    .map(broad_class_map)
    .fillna("Unclassified")
)


# ---------------------------------------------------------
# PPARγ chemotype distribution - thesis Figure 1
# ---------------------------------------------------------

chemotype_counts = (
    pparg_df
    .groupby(["Chemotype", "Broad_Class"])
    .size()
    .reset_index(name="Ligand_Count")
    .sort_values("Ligand_Count")
)

broad_class_colours = {
    "TZD": "#A06AB4",
    "Non-TZD synthetic": "#4C78A8",
    "Fatty acids": "#4E9A7D",
    "Natural product-derived": "#D17A4B",
    "Unclassified": "#9A8775"
}

bar_colours = chemotype_counts["Broad_Class"].map(
    broad_class_colours
)

fig, ax = plt.subplots(figsize=(9, 6))

ax.barh(
    chemotype_counts["Chemotype"],
    chemotype_counts["Ligand_Count"],
    color=bar_colours,
    edgecolor="#4A4A4A",
    linewidth=0.8
)

ax.set_xlabel("Number of ligands")
ax.set_ylabel("Chemotype")

legend_handles = [
    Patch(
        facecolor=colour,
        edgecolor="#4A4A4A",
        label=group
    )
    for group, colour in broad_class_colours.items()
]

ax.legend(
    handles=legend_handles,
    title="Broad class",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

plt.savefig(
    "results/figures/01_pparg_chemotype_distribution.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# =========================================================
# SEAWEED DATABASE
# =========================================================

seaweed_df = pd.read_excel(seaweed_file)
seaweed_df.columns = seaweed_df.columns.str.strip()

seaweed_df, seaweed_fingerprints, seaweed_clusters = cluster_database(
    seaweed_df,
    smiles_col="SMILES",
    id_col="LOTUS_ID",
    cutoff=0.40
)


# ---------------------------------------------------------
# Scaffold classes identified from representative clusters
# ---------------------------------------------------------

# Classes assigned after inspection of representative
# compounds from each major Butina cluster.

cluster_class_map = {
    0: "Carotenoids",
    1: "Nitrogen-containing compounds",
    2: "Fatty acids",
    3: "Terpenoids",
    4: "Steroids",
    5: "Carbohydrates",
    6: "Phenols"
}

seaweed_df["Scaffold_Class"] = (
    seaweed_df["Cluster_ID"]
    .map(cluster_class_map)
    .fillna("Other")
)


# ---------------------------------------------------------
# Seaweed scaffold distribution - thesis Figure 2
# ---------------------------------------------------------

strain_pivot = (
    seaweed_df
    .groupby(["Species", "Scaffold_Class"])
    .size()
    .unstack(fill_value=0)
)

cluster_palette = [
    "#A06AB4",
    "#4C78A8",
    "#4E9A7D",
    "#D17A4B",
    "#D27C9B",
    "#B85C6A",
    "#5F9EA0",
    "#C7A64A"
]

ax = strain_pivot.plot(
    kind="barh",
    stacked=True,
    figsize=(9, 5),
    color=cluster_palette[:len(strain_pivot.columns)],
    edgecolor="#4A4A4A",
    linewidth=0.8
)

ax.set_xlabel("Number of compounds")
ax.set_ylabel("")

plt.legend(
    title="Scaffold class",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    frameon=False
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

plt.savefig(
    "results/figures/02_seaweed_scaffold_distribution.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ---------------------------------------------------------
# UMAP coordinates for cluster assessment
# ---------------------------------------------------------

for dataset, fingerprints, name in [
    (pparg_df, pparg_fingerprints, "PPARγ"),
    (seaweed_df, seaweed_fingerprints, "Seaweed")
]:

    fp_array = []

    for fp in fingerprints:

        array = np.zeros(
            fp.GetNumBits(),
            dtype=int
        )

        DataStructs.ConvertToNumpyArray(
            fp,
            array
        )

        fp_array.append(array)

    fp_array = np.array(fp_array)

    reducer = umap.UMAP(
        n_neighbors=15,
        min_dist=0.1,
        metric="jaccard",
        random_state=42
    )

    embedding = reducer.fit_transform(fp_array)

    dataset["UMAP1"] = embedding[:, 0]
    dataset["UMAP2"] = embedding[:, 1]


# ---------------------------------------------------------
# Save classified databases
# ---------------------------------------------------------

pparg_df.to_excel(
    pparg_output,
    index=False
)

seaweed_df.to_excel(
    seaweed_output,
    index=False
)

print("\nClustering and structural classification complete.")
