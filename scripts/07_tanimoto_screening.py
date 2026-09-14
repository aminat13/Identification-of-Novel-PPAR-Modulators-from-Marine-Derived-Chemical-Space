# 07_tanimoto_screening.py

import pandas as pd
import matplotlib.pyplot as plt

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from matplotlib.colors import LinearSegmentedColormap


# ---------------------------------------------------------
# Load reference and seaweed databases
# ---------------------------------------------------------

reference_df = pd.read_excel(
    "data/processed/PPARG_Agonists_Clustered.xlsx"
)

seaweed_df = pd.read_excel(
    "data/processed/Seaweed_Metabolites_Clustered.xlsx"
)

reference_df.columns = (
    reference_df.columns.str.strip()
)

seaweed_df.columns = (
    seaweed_df.columns.str.strip()
)


# ---------------------------------------------------------
# Create RDKit molecules
# ---------------------------------------------------------

reference_df["ROMol"] = (
    reference_df["Ligand_SMILES"]
    .apply(Chem.MolFromSmiles)
)

seaweed_df["ROMol"] = (
    seaweed_df["SMILES"]
    .apply(Chem.MolFromSmiles)
)

reference_df = reference_df[
    reference_df["ROMol"].notna()
].copy()

seaweed_df = seaweed_df[
    seaweed_df["ROMol"].notna()
].copy()


# ---------------------------------------------------------
# Generate Morgan fingerprints
# ---------------------------------------------------------

morgan_generator = (
    rdFingerprintGenerator
    .GetMorganGenerator(
        radius=2,
        fpSize=2048,
        includeChirality=True
    )
)

reference_df["Morgan_FP"] = (
    reference_df["ROMol"]
    .apply(
        morgan_generator.GetFingerprint
    )
)

seaweed_df["Morgan_FP"] = (
    seaweed_df["ROMol"]
    .apply(
        morgan_generator.GetFingerprint
    )
)


# ---------------------------------------------------------
# Compare every seaweed compound against reference agonists
# ---------------------------------------------------------

similarity_records = []

for _, seaweed_row in seaweed_df.iterrows():

    similarities = (
        DataStructs
        .BulkTanimotoSimilarity(
            seaweed_row["Morgan_FP"],
            reference_df["Morgan_FP"].tolist()
        )
    )

    for reference_index, score in enumerate(
        similarities
    ):

        reference_row = (
            reference_df.iloc[
                reference_index
            ]
        )

        similarity_records.append({
            "LOTUS_ID":
                seaweed_row["LOTUS_ID"],

            "Strain":
                seaweed_row["Species"],

            "Seaweed_Scaffold_Class":
                seaweed_row["Scaffold_Class"],

            "Reference_Ligand_ID":
                reference_row["Ligand_ID"],

            "Reference_Chemotype":
                reference_row["Chemotype"],

            "Reference_Broad_Class":
                reference_row["Broad_Class"],

            "Tanimoto_Similarity":
                float(score),

            "Seaweed_SMILES":
                seaweed_row["SMILES"],

            "Reference_SMILES":
                reference_row["Ligand_SMILES"]
        })


similarity_df = pd.DataFrame(
    similarity_records
)


# ---------------------------------------------------------
# Standardise broad class names used in the thesis
# ---------------------------------------------------------

class_name_map = {
    "Non-TZD synthetic":
        "Non-TZD Synthetic",

    "TZD":
        "TZD Synthetic",

    "Natural product-derived":
        "Natural Product Derived",

    "Fatty acids":
        "Fatty Acid/Lipid-Like"
}

similarity_df["Reference_Broad_Class"] = (
    similarity_df[
        "Reference_Broad_Class"
    ]
    .replace(class_name_map)
)


# ---------------------------------------------------------
# Best reference match for each seaweed metabolite
# ---------------------------------------------------------

best_matches_df = (
    similarity_df
    .sort_values(
        "Tanimoto_Similarity",
        ascending=False
    )
    .drop_duplicates(
        subset="LOTUS_ID",
        keep="first"
    )
    .reset_index(drop=True)
)

best_matches_df["Similarity_Rank"] = (
    best_matches_df[
        "Tanimoto_Similarity"
    ]
    .rank(
        method="dense",
        ascending=False
    )
    .astype(int)
)


# ---------------------------------------------------------
# Maximum similarity by seaweed species and agonist class
# ---------------------------------------------------------

heatmap_df = (
    similarity_df
    .groupby(
        [
            "Strain",
            "Reference_Broad_Class"
        ],
        as_index=False
    )[
        "Tanimoto_Similarity"
    ]
    .max()
    .pivot(
        index="Strain",
        columns="Reference_Broad_Class",
        values="Tanimoto_Similarity"
    )
)

strain_order = [
    "Laminaria digitata",
    "Saccharina latissima",
    "Undaria pinnatifida"
]

class_order = [
    "Non-TZD Synthetic",
    "TZD Synthetic",
    "Natural Product Derived",
    "Fatty Acid/Lipid-Like"
]

heatmap_df = heatmap_df.reindex(
    index=[
        strain
        for strain in strain_order
        if strain in heatmap_df.index
    ],
    columns=[
        agonist_class
        for agonist_class in class_order
        if agonist_class in heatmap_df.columns
    ]
)


# ---------------------------------------------------------
# Thesis Figure 9
# ---------------------------------------------------------

tanimoto_cmap = (
    LinearSegmentedColormap
    .from_list(
        "tanimoto_contrast",
        [
            "#F4F4F4",
            "#5F9EA0",
            "#C7A64A",
            "#B85C6A"
        ]
    )
)

fig, ax = plt.subplots(
    figsize=(9, 4.8)
)

image = ax.imshow(
    heatmap_df.values,
    cmap=tanimoto_cmap,
    vmin=0,
    vmax=1,
    aspect="auto"
)

ax.set_xticks(
    range(len(heatmap_df.columns))
)

ax.set_xticklabels(
    heatmap_df.columns,
    rotation=30,
    ha="right"
)

ax.set_yticks(
    range(len(heatmap_df.index))
)

ax.set_yticklabels(
    heatmap_df.index
)

for row_index in range(
    heatmap_df.shape[0]
):

    for column_index in range(
        heatmap_df.shape[1]
    ):

        value = heatmap_df.iloc[
            row_index,
            column_index
        ]

        if pd.isna(value):
            continue

        text_colour = (
            "white"
            if value >= 0.70
            else "black"
        )

        ax.text(
            column_index,
            row_index,
            f"{value:.2f}",
            ha="center",
            va="center",
            fontweight="bold",
            color=text_colour
        )

colour_bar = fig.colorbar(
    image,
    ax=ax,
    fraction=0.045,
    pad=0.04
)

colour_bar.set_label(
    "Maximum Tanimoto similarity"
)

for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout()

plt.savefig(
    "results/figures/09_tanimoto_similarity_by_species.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ---------------------------------------------------------
# Scaffold contribution to highest similarities
# ---------------------------------------------------------

comparisons = [
    (
        "Undaria pinnatifida",
        "Fatty Acid/Lipid-Like"
    ),
    (
        "Undaria pinnatifida",
        "Natural Product Derived"
    ),
    (
        "Saccharina latissima",
        "Natural Product Derived"
    )
]

scaffold_palette = {
    "Fatty acids": "#4E9A7D",
    "Steroids": "#B85C6A",
    "Terpenoids": "#5F9EA0",
    "Carotenoids": "#4C78A8",
    "Nitrogen-containing compounds": "#D17A4B",
    "Phenols": "#D27C9B",
    "Carbohydrates": "#C7A64A",
    "Other": "#A06AB4"
}


# ---------------------------------------------------------
# Thesis Figure 10
# ---------------------------------------------------------

fig, axes = plt.subplots(
    1,
    3,
    figsize=(18, 6),
    sharex=True
)

for ax, (
    strain,
    agonist_class
) in zip(
    axes,
    comparisons
):

    subset = similarity_df[
        (
            similarity_df["Strain"]
            == strain
        )
        &
        (
            similarity_df[
                "Reference_Broad_Class"
            ]
            == agonist_class
        )
    ]

    summary = (
        subset
        .groupby(
            "Seaweed_Scaffold_Class"
        )[
            "Tanimoto_Similarity"
        ]
        .max()
        .sort_values()
    )

    colours = [
        scaffold_palette.get(
            scaffold,
            "#A06AB4"
        )
        for scaffold in summary.index
    ]

    ax.barh(
        summary.index,
        summary.values,
        color=colours,
        edgecolor="black"
    )

    ax.set_xlim(0, 1)

    ax.set_xlabel(
        "Maximum similarity"
    )

    ax.set_title(
        f"{strain}\n{agonist_class}",
        fontsize=13
    )

    ax.grid(
        axis="x",
        alpha=0.25
    )

plt.tight_layout()

plt.savefig(
    "results/figures/10_scaffold_similarity_contributions.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ---------------------------------------------------------
# Save screening results
# ---------------------------------------------------------

with pd.ExcelWriter(
    "results/tables/tanimoto_screening_results.xlsx"
) as writer:

    similarity_df.to_excel(
        writer,
        sheet_name="All_Comparisons",
        index=False
    )

    best_matches_df.to_excel(
        writer,
        sheet_name="Best_Match_Per_Seaweed",
        index=False
    )

    heatmap_df.to_excel(
        writer,
        sheet_name="Species_Class_Maximum"
    )
