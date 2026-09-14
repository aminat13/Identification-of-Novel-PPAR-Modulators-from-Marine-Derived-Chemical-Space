# 08_screening_integration.py

import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from rdkit import Chem
from scipy.stats import pearsonr
from matplotlib.lines import Line2D


# ---------------------------------------------------------
# File paths
# ---------------------------------------------------------

screening_file = Path(
    "data/processed/Screening_Results.xlsx"
)

docking_file = Path(
    "data/processed/full_seaweed_docking.sdf"
)

seaweed_file = Path(
    "data/processed/Seaweed_Metabolites_Clustered.xlsx"
)

tanimoto_file = Path(
    "results/tables/tanimoto_screening_results.xlsx"
)


# =========================================================
# PHARMACOPHORE SCREENING
# =========================================================

pharmacophore_df = pd.read_excel(
    screening_file,
    sheet_name="pharmacophore_results"
)

pharmacophore_df.columns = (
    pharmacophore_df.columns.str.strip()
)

pharmacophore_df["LOTUS_ID"] = (
    pharmacophore_df["LOTUS_ID"]
    .astype(str)
    .str.strip()
    .str.upper()
)

pharmacophore_df["rmsd"] = pd.to_numeric(
    pharmacophore_df["rmsd"],
    errors="coerce"
)

pharmacophore_df = (
    pharmacophore_df
    .dropna(
        subset=[
            "LOTUS_ID",
            "rmsd"
        ]
    )
    .sort_values("rmsd")
    .drop_duplicates(
        "LOTUS_ID",
        keep="first"
    )
)

pharmacophore_df = (
    pharmacophore_df.rename(
        columns={
            "rmsd":
                "Pharmacophore_RMSD"
        }
    )
)


# ---------------------------------------------------------
# Pharmacophore hits used in thesis Table 2
# ---------------------------------------------------------

pharmacophore_hits = (
    pharmacophore_df[
        pharmacophore_df[
            "Pharmacophore_RMSD"
        ] < 2.0
    ]
    .sort_values(
        "Pharmacophore_RMSD"
    )
)

pharmacophore_hits.to_excel(
    "results/tables/02_pharmacophore_hits.xlsx",
    index=False
)

print(
    "Pharmacophore hits:",
    len(pharmacophore_hits)
)


# =========================================================
# AI-ASSISTED PRIORITISATION
# =========================================================

ai_df = pd.read_excel(
    screening_file,
    sheet_name="AI_results"
)

ai_df.columns = (
    ai_df.columns.str.strip()
)

ai_df["LOTUS_ID"] = (
    ai_df["LOTUS_ID"]
    .astype(str)
    .str.strip()
    .str.upper()
)

ai_df["AI_Priority"] = (
    ai_df["PPARg_Priority"]
    .astype(str)
    .str.strip()
    .str.title()
)

ai_df = ai_df[
    ai_df["AI_Priority"].isin(
        [
            "High",
            "Medium",
            "Low"
        ]
    )
].copy()

ai_df = ai_df.drop_duplicates(
    "LOTUS_ID"
)


# ---------------------------------------------------------
# Add species information
# ---------------------------------------------------------

seaweed_df = pd.read_excel(
    seaweed_file
)

seaweed_df.columns = (
    seaweed_df.columns.str.strip()
)

seaweed_lookup = (
    seaweed_df[
        [
            "LOTUS_ID",
            "Species",
            "Scaffold_Class"
        ]
    ]
    .drop_duplicates(
        "LOTUS_ID"
    )
)

ai_df = ai_df.merge(
    seaweed_lookup,
    on="LOTUS_ID",
    how="left"
)


# ---------------------------------------------------------
# AI priority by species - thesis Figure 12
# ---------------------------------------------------------

ai_counts = (
    ai_df
    .groupby(
        [
            "Species",
            "AI_Priority"
        ]
    )
    .size()
    .unstack(fill_value=0)
)

priority_order = [
    "High",
    "Medium",
    "Low"
]

ai_counts = ai_counts.reindex(
    columns=[
        priority
        for priority in priority_order
        if priority in ai_counts.columns
    ]
)

priority_colours = {
    "High": "#B85C6A",
    "Medium": "#C7A64A",
    "Low": "#5F9EA0"
}

ax = ai_counts.plot(
    kind="bar",
    figsize=(8, 5),
    color=[
        priority_colours[column]
        for column in ai_counts.columns
    ],
    edgecolor="#2F2F2F",
    linewidth=0.8
)

ax.set_xlabel("")
ax.set_ylabel(
    "Number of metabolites"
)

ax.legend(
    title="AI priority",
    frameon=False
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.xticks(
    rotation=20,
    ha="right"
)

plt.tight_layout()

plt.savefig(
    "results/figures/12_ai_priority_by_species.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# =========================================================
# MOLECULAR DOCKING
# =========================================================

supplier = Chem.SDMolSupplier(
    str(docking_file),
    removeHs=False
)

docking_records = []

for mol in supplier:

    if mol is None:
        continue

    props = mol.GetPropsAsDict()

    lotus_id = None

    for key in [
        "LOTUS_ID",
        "LOTUSID",
        "ID",
        "Name"
    ]:

        if key in props:
            lotus_id = props[key]
            break

    if lotus_id is None:
        lotus_id = mol.GetProp(
            "_Name"
        )

    docking_score = None

    for key in [
        "S",
        "Score",
        "Docking_score",
        "Docking_Score",
        "dG",
        "Affinity"
    ]:

        if key in props:
            docking_score = props[key]
            break

    docking_records.append({
        "LOTUS_ID": lotus_id,
        "Docking_Score": docking_score
    })


docking_df = pd.DataFrame(
    docking_records
)

docking_df["LOTUS_ID"] = (
    docking_df["LOTUS_ID"]
    .astype(str)
    .str.strip()
    .str.upper()
)

docking_df["Docking_Score"] = (
    pd.to_numeric(
        docking_df["Docking_Score"],
        errors="coerce"
    )
)

# Retain the best-scoring pose for each compound
docking_df = (
    docking_df
    .dropna(
        subset=[
            "LOTUS_ID",
            "Docking_Score"
        ]
    )
    .sort_values(
        "Docking_Score"
    )
    .drop_duplicates(
        "LOTUS_ID",
        keep="first"
    )
)


# =========================================================
# TANIMOTO RESULTS
# =========================================================

tanimoto_df = pd.read_excel(
    tanimoto_file,
    sheet_name="Best_Match_Per_Seaweed"
)

tanimoto_df.columns = (
    tanimoto_df.columns.str.strip()
)

tanimoto_df["LOTUS_ID"] = (
    tanimoto_df["LOTUS_ID"]
    .astype(str)
    .str.strip()
    .str.upper()
)

tanimoto_df = (
    tanimoto_df[
        [
            "LOTUS_ID",
            "Tanimoto_Similarity"
        ]
    ]
    .drop_duplicates(
        "LOTUS_ID"
    )
)


# =========================================================
# COMBINE ALL SCREENING RESULTS
# =========================================================

combined_df = (
    seaweed_lookup
    .merge(
        tanimoto_df,
        on="LOTUS_ID",
        how="left"
    )
    .merge(
        pharmacophore_df[
            [
                "LOTUS_ID",
                "Pharmacophore_RMSD"
            ]
        ],
        on="LOTUS_ID",
        how="left"
    )
    .merge(
        ai_df[
            [
                "LOTUS_ID",
                "AI_Priority"
            ]
        ],
        on="LOTUS_ID",
        how="left"
    )
    .merge(
        docking_df,
        on="LOTUS_ID",
        how="left"
    )
)


# ---------------------------------------------------------
# Pharmacophore RMSD vs docking score
# ---------------------------------------------------------

correlation_df = (
    combined_df
    .dropna(
        subset=[
            "Pharmacophore_RMSD",
            "Docking_Score"
        ]
    )
)

pearson_r, pearson_p = pearsonr(
    correlation_df[
        "Pharmacophore_RMSD"
    ],
    correlation_df[
        "Docking_Score"
    ]
)

print(
    f"Pearson r = {pearson_r:.3f}"
)

print(
    f"p = {pearson_p:.4f}"
)


# =========================================================
# INTEGRATED SCREENING - THESIS FIGURE 14
# =========================================================

plot_df = combined_df.dropna(
    subset=[
        "Pharmacophore_RMSD",
        "Docking_Score",
        "Tanimoto_Similarity",
        "AI_Priority"
    ]
).copy()

markers = {
    "High": "o",
    "Medium": "^",
    "Low": "s"
}

fig, ax = plt.subplots(
    figsize=(8, 6)
)

scatter = None

for priority in [
    "High",
    "Medium",
    "Low"
]:

    subset = plot_df[
        plot_df["AI_Priority"]
        == priority
    ]

    if subset.empty:
        continue

    scatter = ax.scatter(
        subset[
            "Pharmacophore_RMSD"
        ],
        subset[
            "Docking_Score"
        ],
        c=subset[
            "Tanimoto_Similarity"
        ],
        cmap="plasma",
        marker=markers[priority],
        s=90,
        edgecolor="black",
        linewidth=0.6,
        alpha=0.9,
        vmin=plot_df[
            "Tanimoto_Similarity"
        ].min(),
        vmax=plot_df[
            "Tanimoto_Similarity"
        ].max()
    )


# ---------------------------------------------------------
# Tanimoto colour scale
# ---------------------------------------------------------

colour_bar = fig.colorbar(
    scatter,
    ax=ax
)

colour_bar.set_label(
    "Tanimoto similarity"
)


# ---------------------------------------------------------
# AI priority shape legend
# ---------------------------------------------------------

priority_legend = [
    Line2D(
        [0],
        [0],
        marker=markers[priority],
        linestyle="None",
        markerfacecolor="lightgrey",
        markeredgecolor="black",
        markersize=9,
        label=priority
    )
    for priority in [
        "High",
        "Medium",
        "Low"
    ]
]

ax.legend(
    handles=priority_legend,
    title="AI priority",
    frameon=False,
    loc="upper right"
)

ax.set_xlabel(
    "Pharmacophore RMSD (Å)"
)

ax.set_ylabel(
    "Docking score"
)

ax.grid(
    alpha=0.25
)


# ---------------------------------------------------------
# Add Pearson correlation
# ---------------------------------------------------------

ax.text(
    0.03,
    0.04,
    f"Pearson r = {pearson_r:.3f}\n"
    f"p = {pearson_p:.3f}",
    transform=ax.transAxes,
    fontsize=9,
    bbox={
        "facecolor": "white",
        "alpha": 0.8,
        "edgecolor": "lightgrey"
    }
)

plt.tight_layout()

plt.savefig(
    "results/figures/14_integrated_screening.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# =========================================================
# TOP 10 DOCKING CANDIDATES - THESIS TABLE 3
# =========================================================

top10_docking = (
    combined_df
    .dropna(
        subset=[
            "Docking_Score"
        ]
    )
    .sort_values(
        "Docking_Score"
    )
    .head(10)
)

top10_docking = top10_docking[
    [
        "LOTUS_ID",
        "Species",
        "Scaffold_Class",
        "Tanimoto_Similarity",
        "Pharmacophore_RMSD",
        "Docking_Score",
        "AI_Priority"
    ]
]

print(
    "\nTop 10 docking candidates:"
)

print(top10_docking)


# ---------------------------------------------------------
# Save final screening tables
# ---------------------------------------------------------

with pd.ExcelWriter(
    "results/tables/integrated_screening_results.xlsx"
) as writer:

    combined_df.to_excel(
        writer,
        sheet_name="All_Compounds",
        index=False
    )

    pharmacophore_hits.to_excel(
        writer,
        sheet_name="Pharmacophore_Hits",
        index=False
    )

    top10_docking.to_excel(
        writer,
        sheet_name="Top_10_Docking",
        index=False
    )
