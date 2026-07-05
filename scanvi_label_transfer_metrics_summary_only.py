# ------------------------------------------------------------
# scANVI mouse-reference label-transfer metrics
# Summary tables only; no bar plots or heatmaps
#
# This block can be run after the main scANVI label-transfer script,
# after the following fields have been created:
#   adata.obsm["X_scANVI"]
#   adata.obs["species"]
#   adata.obs["scanvi_final_label"]
#   adata.obs["scanvi_confidence"]
# ------------------------------------------------------------

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

# -----------------------------
# Settings
# -----------------------------

results_path = Path("~/snRNA_seq_integ/cross_species_tanycytes_full").expanduser()
outdir = results_path / "mouse_reference_label_transfer_metrics"
outdir.mkdir(parents=True, exist_ok=True)

latent_key = "X_scANVI"
species_key = "species"
label_key = "scanvi_final_label"
confidence_key = "scanvi_confidence"

mouse_species = "mouse"
human_species = "human"
k = 30

# -----------------------------
# Basic checks
# -----------------------------

required_obs = [species_key, label_key]
missing_obs = [key for key in required_obs if key not in adata.obs.columns]

if missing_obs:
    raise ValueError(f"Missing required obs columns: {missing_obs}")

if latent_key not in adata.obsm:
    raise ValueError(f"Missing latent representation: adata.obsm['{latent_key}']")

adata.obs[species_key] = adata.obs[species_key].astype(str)
adata.obs[label_key] = adata.obs[label_key].astype(str)

mouse_mask = adata.obs[species_key] == mouse_species
human_mask = adata.obs[species_key] == human_species

if mouse_mask.sum() == 0:
    raise ValueError(f"No cells found for species == '{mouse_species}'")

if human_mask.sum() == 0:
    raise ValueError(f"No cells found for species == '{human_species}'")

mouse = adata[mouse_mask].copy()
human = adata[human_mask].copy()

mouse_labels = mouse.obs[label_key].astype(str).values
human_labels = human.obs[label_key].astype(str).values

# -----------------------------
# 1. Human-to-mouse nearest-neighbor label agreement
# -----------------------------
# For each human cell, find its k nearest mouse neighbors in the
# scANVI latent space and calculate the fraction of mouse neighbors
# with the same final transferred label.

nn_mouse = NearestNeighbors(
    n_neighbors=k,
    metric="cosine"
)

nn_mouse.fit(mouse.obsm[latent_key])

mouse_neighbor_idx = nn_mouse.kneighbors(
    human.obsm[latent_key],
    return_distance=False
)

agreement_key = f"mouse_label_agreement_k{k}"

human.obs[agreement_key] = [
    np.mean(mouse_labels[idx] == human_labels[i])
    for i, idx in enumerate(mouse_neighbor_idx)
]

human.obs["dominant_mouse_neighbor_label"] = [
    pd.Series(mouse_labels[idx]).value_counts().idxmax()
    for idx in mouse_neighbor_idx
]

# Add human-only metrics back to the full object
adata.obs[agreement_key] = np.nan
adata.obs["dominant_mouse_neighbor_label"] = pd.NA

adata.obs.loc[human.obs_names, agreement_key] = human.obs[agreement_key]
adata.obs.loc[
    human.obs_names,
    "dominant_mouse_neighbor_label"
] = human.obs["dominant_mouse_neighbor_label"]

agreement_summary = (
    human.obs
    .groupby(label_key, observed=True)[agreement_key]
    .agg(["mean", "median", "std", "count"])
    .reset_index()
    .rename(
        columns={
            "mean": "mean_human_mouse_agreement",
            "median": "median_human_mouse_agreement",
            "std": "sd_human_mouse_agreement",
            "count": "n_human_cells"
        }
    )
    .sort_values("mean_human_mouse_agreement", ascending=False)
)

agreement_summary.to_csv(
    outdir / f"human_to_mouse_label_agreement_by_state_k{k}.csv",
    index=False
)

# -----------------------------
# 2. scANVI confidence summary for human cells
# -----------------------------

if confidence_key in human.obs.columns:
    confidence_summary = (
        human.obs
        .groupby(label_key, observed=True)[confidence_key]
        .agg(["mean", "median", "std", "count"])
        .reset_index()
        .rename(
            columns={
                "mean": "mean_scanvi_confidence",
                "median": "median_scanvi_confidence",
                "std": "sd_scanvi_confidence",
                "count": "n_human_cells_confidence"
            }
        )
    )

    confidence_summary.to_csv(
        outdir / "human_scanvi_confidence_by_mouse_reference_label.csv",
        index=False
    )
else:
    confidence_summary = None
    print(f"Warning: '{confidence_key}' not found; skipping confidence summary.")

# -----------------------------
# 3. Human assigned label vs dominant mouse-neighbor label
# -----------------------------

confusion_counts = pd.crosstab(
    human.obs[label_key],
    human.obs["dominant_mouse_neighbor_label"]
)

confusion_fraction = pd.crosstab(
    human.obs[label_key],
    human.obs["dominant_mouse_neighbor_label"],
    normalize="index"
)

confusion_counts.to_csv(
    outdir / f"human_label_vs_dominant_mouse_neighbor_label_counts_k{k}.csv"
)

confusion_fraction.to_csv(
    outdir / f"human_label_vs_dominant_mouse_neighbor_label_fraction_k{k}.csv"
)

# -----------------------------
# 4. Species composition per final label
# -----------------------------

species_composition_counts = pd.crosstab(
    adata.obs[label_key],
    adata.obs[species_key]
)

species_composition_fraction = pd.crosstab(
    adata.obs[label_key],
    adata.obs[species_key],
    normalize="index"
)

species_composition_counts.to_csv(
    outdir / "species_composition_counts_by_final_label.csv"
)

species_composition_fraction.to_csv(
    outdir / "species_composition_fraction_by_final_label.csv"
)

# -----------------------------
# 5. Species iLISI in scANVI latent space
# -----------------------------
# iLISI is calculated from each cell's k nearest neighbors.
# Values near 1 indicate mostly same-species neighborhoods;
# values approaching 2 indicate stronger mouse-human mixing.

species_values = adata.obs[species_key].astype(str).values

nn_all = NearestNeighbors(
    n_neighbors=k + 1,
    metric="cosine"
)

nn_all.fit(adata.obsm[latent_key])
idx_all = nn_all.kneighbors(return_distance=False)

ilisi_key = f"iLISI_species_k{k}"
ilisi = np.zeros(adata.n_obs)

for i in range(adata.n_obs):
    neighbors = idx_all[i, 1:]
    neighbor_species = species_values[neighbors]
    proportions = pd.Series(neighbor_species).value_counts(normalize=True)
    ilisi[i] = 1.0 / np.sum(proportions ** 2)

adata.obs[ilisi_key] = ilisi

ilisi_summary = (
    adata.obs
    .groupby(label_key, observed=True)[ilisi_key]
    .agg(["mean", "median", "std", "count"])
    .reset_index()
    .rename(
        columns={
            "mean": "mean_iLISI",
            "median": "median_iLISI",
            "std": "sd_iLISI",
            "count": "n_total_cells"
        }
    )
    .sort_values("mean_iLISI", ascending=False)
)

ilisi_summary.to_csv(
    outdir / f"species_iLISI_by_final_label_k{k}.csv",
    index=False
)

# -----------------------------
# 6. Combined manuscript summary table
# -----------------------------

mouse_counts = (
    mouse.obs
    .groupby(label_key, observed=True)
    .size()
    .rename("n_mouse_cells")
    .reset_index()
)

human_counts = (
    human.obs
    .groupby(label_key, observed=True)
    .size()
    .rename("n_human_cells")
    .reset_index()
)

summary_table = (
    mouse_counts
    .merge(human_counts, on=label_key, how="outer")
    .merge(
        agreement_summary[
            [
                label_key,
                "mean_human_mouse_agreement",
                "median_human_mouse_agreement",
                "sd_human_mouse_agreement"
            ]
        ],
        on=label_key,
        how="left"
    )
    .merge(
        ilisi_summary[
            [
                label_key,
                "mean_iLISI",
                "median_iLISI",
                "sd_iLISI"
            ]
        ],
        on=label_key,
        how="left"
    )
)

if confidence_summary is not None:
    summary_table = summary_table.merge(
        confidence_summary[
            [
                label_key,
                "mean_scanvi_confidence",
                "median_scanvi_confidence",
                "sd_scanvi_confidence"
            ]
        ],
        on=label_key,
        how="left"
    )

summary_table["n_mouse_cells"] = (
    summary_table["n_mouse_cells"].fillna(0).astype(int)
)

summary_table["n_human_cells"] = (
    summary_table["n_human_cells"].fillna(0).astype(int)
)

summary_table = summary_table.sort_values(label_key).round(3)

summary_table.to_csv(
    outdir / f"summary_scanvi_confidence_mouse_agreement_iLISI_k{k}.csv",
    index=False
)

# -----------------------------
# 7. Global metrics table
# -----------------------------

global_metrics = {
    "n_mouse_cells": int(mouse_mask.sum()),
    "n_human_cells": int(human_mask.sum()),
    "mean_human_to_mouse_label_agreement": human.obs[agreement_key].mean(),
    "median_human_to_mouse_label_agreement": human.obs[agreement_key].median(),
    "mean_species_iLISI": adata.obs[ilisi_key].mean(),
    "median_species_iLISI": adata.obs[ilisi_key].median()
}

if confidence_key in human.obs.columns:
    global_metrics["mean_human_scanvi_confidence"] = human.obs[confidence_key].mean()
    global_metrics["median_human_scanvi_confidence"] = human.obs[confidence_key].median()

global_metrics = pd.DataFrame(
    list(global_metrics.items()),
    columns=["metric", "value"]
)

global_metrics.to_csv(
    outdir / f"global_mouse_reference_label_transfer_metrics_k{k}.csv",
    index=False
)

# -----------------------------
# 8. Save updated metadata/object
# -----------------------------

adata.obs.to_csv(
    outdir / f"scanvi_metadata_with_mouse_reference_metrics_k{k}.csv"
)

adata.write_h5ad(
    outdir / f"tany_mouse_human_scanvi_final_with_metrics_k{k}.h5ad"
)

print("Saved summary tables to:", outdir)
print(summary_table)
print(global_metrics)
