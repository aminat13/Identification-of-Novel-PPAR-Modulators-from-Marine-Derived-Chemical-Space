# 06_structure_activity_analysis.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ---------------------------------------------------------
# Load EC50-annotated PPARγ agonists
# ---------------------------------------------------------

pparg_file = "data/processed/PPARG_Agonists_Clustered.xlsx"

pparg_df = pd.read_excel(pparg_file)
pparg_df.columns = pparg_df.columns.str.strip()


# ---------------------------------------------------------
# Convert EC50 to pEC50
# ---------------------------------------------------------

pparg_df["EC50_(nM)"] = pd.to_numeric(
    pparg_df["EC50_(nM)"],
    errors="coerce"
)

sar_df = pparg_df.dropna(
    subset=["EC50_(nM)"]
).copy()

sar_df["pEC50"] = -np.log10(
    sar_df["EC50_(nM)"] * 1e-9
)


# ---------------------------------------------------------
# Define activity groups
# ---------------------------------------------------------

activity_threshold = 6.5

sar_df["Activity_Status"] = np.where(
    sar_df["pEC50"] >= activity_threshold,
    "High",
    "Low"
)

print(
    sar_df["Activity_Status"]
    .value_counts()
)


# ---------------------------------------------------------
# Compare activity across major agonist classes
# ---------------------------------------------------------

class_column = (
    "Class"
    if "Class" in sar_df.columns
    else "Broad_Class"
)

high_counts = (
    sar_df[
        sar_df["Activity_Status"] == "High"
    ][class_column]
    .value_counts()
)

low_counts = (
    sar_df[
        sar_df["Activity_Status"] == "Low"
    ][class_column]
    .value_counts()
)

comparison_df = pd.concat(
    [high_counts, low_counts],
    axis=1
)

comparison_df.columns = [
    "High activity",
    "Low activity"
]

class_order = [
    "Non-TZD Synthetic",
    "TZD Synthetic",
    "Natural Product Derived",
    "Fatty Acid/Lipid-Like"
]

# Allow for slightly different labels in the processed database
class_name_map = {
    "Non-TZD synthetic": "Non-TZD Synthetic",
    "TZD": "TZD Synthetic",
    "Natural product-derived": "Natural Product Derived",
    "Fatty acids": "Fatty Acid/Lipid-Like"
}

comparison_df = comparison_df.rename(
    index=class_name_map
)

comparison_df = (
    comparison_df
    .reindex(class_order)
    .fillna(0)
    .astype(int)
)


# ---------------------------------------------------------
# Thesis Figure 6
# ---------------------------------------------------------

x = np.arange(
    len(comparison_df)
)

width = 0.35

fig, ax = plt.subplots(
    figsize=(8, 5)
)

bars_high = ax.bar(
    x - width / 2,
    comparison_df["High activity"],
    width,
    color="#A06AB4",
    edgecolor="#2F2F2F",
    linewidth=1.2,
    label="High activity"
)

bars_low = ax.bar(
    x + width / 2,
    comparison_df["Low activity"],
    width,
    color="#4E9A7D",
    edgecolor="#2F2F2F",
    linewidth=1.2,
    label="Low activity"
)

for bars in [bars_high, bars_low]:

    for bar in bars:

        height = bar.get_height()

        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            int(height),
            ha="center",
            va="bottom",
            fontsize=10
        )

ax.set_xticks(x)

ax.set_xticklabels(
    comparison_df.index,
    rotation=20,
    ha="right"
)

ax.set_ylabel(
    "Number of ligands"
)

ax.legend(
    frameon=False
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

plt.savefig(
    "results/figures/06_activity_class_comparison.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# =========================================================
# MANUALLY ANNOTATED FUNCTIONAL GROUPS
# =========================================================

# Functional-group annotations were produced during manual
# SAR review and are stored separately.

functional_group_file = (
    "data/processed/PPARG_SAR_Functional_Groups.xlsx"
)

functional_df = pd.read_excel(
    functional_group_file
)

functional_df.columns = (
    functional_df.columns.str.strip()
)


# ---------------------------------------------------------
# Add pEC50/activity classification
# ---------------------------------------------------------

activity_lookup = (
    sar_df[
        [
            "Ligand_ID",
            "pEC50",
            "Activity_Status"
        ]
    ]
    .drop_duplicates("Ligand_ID")
)

functional_df = functional_df.drop(
    columns=[
        "pEC50",
        "Activity_Status"
    ],
    errors="ignore"
)

functional_df = functional_df.merge(
    activity_lookup,
    on="Ligand_ID",
    how="left"
)


# ---------------------------------------------------------
# Functional groups assessed in the thesis
# ---------------------------------------------------------

functional_groups = [
    "Carboxylic acid",
    "Hydroxyl",
    "Amide",
    "Sulfonamide",
    "Ether",
    "Ketone",
    "Additional aromatic ring",
    "Heteroaromatic ring",
    "Halogen",
    "CF₃",
    "Flexible linker",
    "Double bond"
]

functional_df[functional_groups] = (
    functional_df[functional_groups]
    .fillna(0)
    .replace({
        "✓": 1,
        "✔": 1,
        "Yes": 1,
        "yes": 1,
        "Present": 1
    })
    .astype(int)
)


# ---------------------------------------------------------
# Functional groups by activity - thesis Figure 7
# ---------------------------------------------------------

functional_by_activity = (
    functional_df
    .groupby("Activity_Status")[
        functional_groups
    ]
    .sum()
)

print(functional_by_activity)

functional_by_activity.T.plot(
    kind="barh",
    figsize=(12, 7),
    color=[
        "#B85C6A",
        "#6C8EAD"
    ],
    edgecolor="#2F2F2F",
    linewidth=1.0
)

ax = plt.gca()

ax.set_xlabel(
    "Number of ligands"
)

ax.set_ylabel("")

ax.legend(
    title="Activity status",
    frameon=False
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

plt.savefig(
    "results/figures/07_functional_groups_by_activity.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ---------------------------------------------------------
# Save SAR results
# ---------------------------------------------------------

sar_df.to_excel(
    "results/tables/pparg_activity_classification.xlsx",
    index=False
)

functional_df.to_excel(
    "results/tables/pparg_sar_functional_groups.xlsx",
    index=False
)
