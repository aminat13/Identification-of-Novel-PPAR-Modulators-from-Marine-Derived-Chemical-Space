# 05_binding_interaction_analysis.py

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


# ---------------------------------------------------------
# Load ligand-residue interactions exported from MOE
# ---------------------------------------------------------

interaction_file = "data/processed/PPARG_Ligand_Interactions.xlsx"

interaction_df = pd.read_excel(interaction_file)
interaction_df.columns = interaction_df.columns.str.strip()


# ---------------------------------------------------------
# Create residue labels
# ---------------------------------------------------------

interaction_df = interaction_df[
    [
        "Residue_Name",
        "Residue_Number",
        "Interaction_Type"
    ]
].dropna().copy()

interaction_df["Residue"] = (
    interaction_df["Residue_Name"]
    .astype(str)
    .str.strip()
    .str.upper()
    +
    pd.to_numeric(
        interaction_df["Residue_Number"],
        errors="coerce"
    )
    .astype("Int64")
    .astype(str)
)


# ---------------------------------------------------------
# Standardise interaction names
# ---------------------------------------------------------

interaction_name_map = {
    "H-bond acceptor": "H-acceptor",
    "Hydrogen bond acceptor": "H-acceptor",

    "H-bond donor": "H-donor",
    "Hydrogen bond donor": "H-donor",

    "Pi-H": "pi-H",
    "PI-H": "pi-H",
    "π-H": "pi-H",

    "H-pi": "H-pi",
    "H-π": "H-pi",

    "ionic": "Ionic"
}

interaction_df["Interaction_Type"] = (
    interaction_df["Interaction_Type"]
    .astype(str)
    .str.strip()
    .replace(interaction_name_map)
)


# ---------------------------------------------------------
# Count interaction types for conserved residues
# ---------------------------------------------------------

residue_order = [
    "ARG288",
    "CYS285",
    "HIS323",
    "HIS449",
    "ILE341",
    "LEU340",
    "MET364",
    "SER289",
    "TYR473"
]

interaction_order = [
    "H-acceptor",
    "H-donor",
    "H-pi",
    "Ionic",
    "pi-H"
]

heatmap_data = pd.crosstab(
    interaction_df["Residue"],
    interaction_df["Interaction_Type"]
)

heatmap_data = heatmap_data.reindex(
    index=residue_order,
    columns=interaction_order,
    fill_value=0
)

print(heatmap_data)


# ---------------------------------------------------------
# Interaction-frequency heatmap - thesis Figure 5
# ---------------------------------------------------------

interaction_cmap = LinearSegmentedColormap.from_list(
    "interaction_frequency",
    [
        "#F4F4F4",
        "#5F9EA0",
        "#C7A64A",
        "#B85C6A"
    ]
)

fig, ax = plt.subplots(figsize=(8.5, 6))

image = ax.imshow(
    heatmap_data.values,
    cmap=interaction_cmap,
    aspect="auto",
    interpolation="nearest"
)

ax.set_xticks(
    range(len(heatmap_data.columns))
)

ax.set_xticklabels(
    heatmap_data.columns,
    fontsize=11
)

ax.set_yticks(
    range(len(heatmap_data.index))
)

ax.set_yticklabels(
    heatmap_data.index,
    fontsize=11
)

ax.set_xlabel(
    "Interaction type",
    fontsize=13
)

ax.set_ylabel(
    "Residue",
    fontsize=13
)

colour_bar = fig.colorbar(
    image,
    ax=ax,
    fraction=0.045,
    pad=0.04
)

colour_bar.set_label(
    "Interaction frequency",
    fontsize=11
)

for spine in ax.spines.values():
    spine.set_visible(False)

ax.tick_params(
    axis="both",
    length=0
)

plt.tight_layout()

plt.savefig(
    "results/figures/05_residue_interaction_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ---------------------------------------------------------
# Save interaction counts
# ---------------------------------------------------------

heatmap_data.to_excel(
    "results/tables/residue_interaction_frequencies.xlsx"
)
