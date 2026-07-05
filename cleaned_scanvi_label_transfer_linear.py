#!/usr/bin/env python
# coding: utf-8

"""
Cross-species tanycyte scVI/scANVI integration and label transfer.

Linear script version without helper functions.

Key behavior:
- Train scVI on mouse + human cells.
- Train scANVI using mouse reference labels and human cells labelled as "Unknown".
- Preserve curated mouse labels in the final output.
- Assign scANVI-predicted labels only to human cells.
"""

from pathlib import Path
import json
import gc

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from sklearn.neighbors import NearestNeighbors


# -----------------------------
# Configuration
# -----------------------------

PARAM_FILE = Path(
    "~/snRNA_seq_integ/cross_species_tanycytes_full/"
    "parameters_cross_species_tanycytes_scvi_full.json"
).expanduser()

LABEL_KEY = "transfer_label_2"
UNKNOWN_LABEL = "Unknown"
SPECIES_KEY = "species"
SAMPLE_KEY = "sample_id"

SCVI_LATENT_KEY = "X_scVI"
SCANVI_LATENT_KEY = "X_scANVI"

LEIDEN_RESOLUTIONS = [0.4, 0.6, 0.8]
N_NEIGHBORS = 30
USE_CUDA = True

print("scvi-tools version:", scvi.__version__)


# -----------------------------
# Load parameters
# -----------------------------

if not PARAM_FILE.exists():
    raise FileNotFoundError(f"Parameter file not found: {PARAM_FILE}")

with open(PARAM_FILE, "r") as f:
    params = json.load(f)

required_params = [
    "harmonization_folder_path",
    "data_filepath_full",
    "feature_set_file",
    "new_name_suffix",
    "global_seed",
    "batch_var",
    "n_layers",
    "n_latent",
    "n_hidden",
    "dropout_rate",
    "dispersion",
    "gene_likelihood",
    "max_epochs",
    "early_stopping",
    "k_param",
]

missing_params = [x for x in required_params if x not in params]
if missing_params:
    raise ValueError(f"Missing parameters in JSON: {missing_params}")

results_path = Path(params["harmonization_folder_path"]).expanduser()
results_path.mkdir(parents=True, exist_ok=True)

data_file = Path(params["data_filepath_full"]).expanduser()
feature_file = Path(params["feature_set_file"]).expanduser()
suffix = params["new_name_suffix"]
seed = int(params["global_seed"])

scvi.settings.seed = seed
np.random.seed(seed)

print("\nInput/output paths")
print("AnnData:", data_file)
print("Feature set:", feature_file)
print("Results:", results_path)


# -----------------------------
# Read AnnData and subset features
# -----------------------------

adata = sc.read_h5ad(data_file)

adata.obs_names = adata.obs_names.astype(str)
adata.var_names = adata.var_names.astype(str)
adata.var_names_make_unique()

print("\nLoaded AnnData:")
print(adata)
print("Metadata columns:", list(adata.obs.columns))

with open(feature_file, "r") as f:
    feature_set = json.load(f)

feature_set = [gene for gene in feature_set if gene in adata.var_names]
print("Feature-set genes found in AnnData:", len(feature_set))

if len(feature_set) == 0:
    raise ValueError("No feature-set genes found in adata.var_names.")

adata = adata[:, feature_set].copy()


# -----------------------------
# Prepare counts layer for scvi-tools
# -----------------------------

print("\nAvailable layers before counts setup:", list(adata.layers.keys()))
print("Has raw:", adata.raw is not None)

if "counts" in adata.layers:
    adata.layers["scvi_counts"] = adata.layers["counts"].copy()
elif adata.raw is not None:
    raw_adata = adata.raw.to_adata()
    raw_adata = raw_adata[:, adata.var_names].copy()
    adata.layers["scvi_counts"] = raw_adata.X.copy()
else:
    adata.layers["scvi_counts"] = adata.X.copy()

print("Available layers after counts setup:", list(adata.layers.keys()))


# -----------------------------
# Metadata checks
# -----------------------------

required_metadata = [SPECIES_KEY, SAMPLE_KEY, LABEL_KEY]
missing_metadata = [col for col in required_metadata if col not in adata.obs.columns]

if missing_metadata:
    raise ValueError(f"Missing metadata columns: {missing_metadata}")

for col in required_metadata:
    print(f"\n{col}:")
    print(adata.obs[col].value_counts(dropna=False))


# -----------------------------
# Train scVI
# -----------------------------

scvi.model.SCVI.setup_anndata(
    adata,
    layer="scvi_counts",
    batch_key=SPECIES_KEY,
    categorical_covariate_keys=[SAMPLE_KEY],
)

vae = scvi.model.SCVI(
    adata,
    n_layers=int(params["n_layers"]),
    n_latent=int(params["n_latent"]),
    n_hidden=int(params["n_hidden"]),
    dropout_rate=float(params["dropout_rate"]),
    dispersion=str(params["dispersion"]),
    gene_likelihood=str(params["gene_likelihood"]),
)

vae.train(
    max_epochs=int(params["max_epochs"]),
    early_stopping=bool(params["early_stopping"]),
    accelerator="gpu" if USE_CUDA else "cpu",
    devices=1,
    batch_size=1024,
)

adata.obsm[SCVI_LATENT_KEY] = vae.get_latent_representation()

scvi_latent = pd.DataFrame(
    adata.obsm[SCVI_LATENT_KEY],
    index=adata.obs_names,
    columns=[f"scVI_{i + 1}" for i in range(adata.obsm[SCVI_LATENT_KEY].shape[1])],
)

scvi_latent_file = results_path / f"{suffix}_scVI_reduction.txt"
scvi_latent.to_csv(scvi_latent_file, sep="\t")
print("Saved scVI latent space:", scvi_latent_file)


# -----------------------------
# scVI neighbors, UMAP, and Leiden clustering
# -----------------------------

sc.pp.neighbors(
    adata,
    use_rep=SCVI_LATENT_KEY,
    n_neighbors=int(params["k_param"]),
    metric=params.get("dist_type", "cosine"),
)

sc.tl.umap(adata, random_state=seed)

for res in LEIDEN_RESOLUTIONS:
    sc.tl.leiden(
        adata,
        resolution=res,
        key_added=f"scvi_leiden_{res}",
        flavor="igraph",
        directed=False,
        n_iterations=2,
    )

for res in LEIDEN_RESOLUTIONS:
    key = f"scvi_leiden_{res}"
    print(f"\nSpecies overlap for {key}")
    print(pd.crosstab(adata.obs[key], adata.obs[SPECIES_KEY], normalize="index"))


# -----------------------------
# Train scANVI
# -----------------------------

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    vae,
    labels_key=LABEL_KEY,
    unlabeled_category=UNKNOWN_LABEL,
)

scanvi_model.train(
    max_epochs=500,
    accelerator="gpu" if USE_CUDA else "cpu",
    devices=1,
)

adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation()


# -----------------------------
# scANVI label transfer while preserving mouse reference labels
# -----------------------------

adata.obs["reference_input_label"] = adata.obs[LABEL_KEY].astype(str)

mouse_mask = adata.obs[SPECIES_KEY].astype(str).eq("mouse")
human_mask = adata.obs[SPECIES_KEY].astype(str).eq("human")
unknown_mask = adata.obs[LABEL_KEY].astype(str).eq(UNKNOWN_LABEL)

# Predict all cells for diagnostics only.
adata.obs["scanvi_model_prediction"] = scanvi_model.predict(adata)

pred_probs = scanvi_model.predict(adata, soft=True)
adata.obs["scanvi_confidence"] = pred_probs.max(axis=1)

# Final label strategy:
# - Mouse cells retain original curated reference labels.
# - Human cells receive scANVI-predicted mT labels.
adata.obs["scanvi_final_label"] = adata.obs["reference_input_label"].copy()
adata.obs.loc[human_mask, "scanvi_final_label"] = scanvi_model.predict(
    adata[human_mask].copy()
)

# Alternative option if you prefer to label all Unknown cells rather than all human cells:
# adata.obs.loc[unknown_mask, "scanvi_final_label"] = scanvi_model.predict(
#     adata[unknown_mask].copy()
# )

adata.obs["scanvi_final_label"] = pd.Categorical(
    adata.obs["scanvi_final_label"],
    categories=[f"mT.{i}" for i in range(1, 9)],
    ordered=True,
)

print("\nMouse reference labels preserved in final labels:")
print(pd.crosstab(
    adata.obs.loc[mouse_mask, "reference_input_label"],
    adata.obs.loc[mouse_mask, "scanvi_final_label"],
    normalize="index",
))

print("\nMouse labels vs model prediction, diagnostic only:")
print(pd.crosstab(
    adata.obs.loc[mouse_mask, "reference_input_label"],
    adata.obs.loc[mouse_mask, "scanvi_model_prediction"],
    normalize="index",
))

print("\nHuman projected labels:")
print(adata.obs.loc[human_mask, "scanvi_final_label"].value_counts().sort_index())

print("\nHuman prediction confidence:")
print(
    adata.obs.loc[human_mask]
    .groupby("scanvi_final_label", observed=False)["scanvi_confidence"]
    .describe()
)


# -----------------------------
# scANVI neighbors, UMAP, and Leiden clustering
# -----------------------------

sc.pp.neighbors(
    adata,
    use_rep=SCANVI_LATENT_KEY,
    n_neighbors=N_NEIGHBORS,
    metric="cosine",
)

sc.tl.umap(adata, random_state=seed)

for res in LEIDEN_RESOLUTIONS:
    sc.tl.leiden(
        adata,
        resolution=res,
        key_added=f"scanvi_leiden_{res}",
        flavor="igraph",
        directed=False,
        n_iterations=2,
    )

sc.pl.umap(
    adata,
    color=[SPECIES_KEY, "scanvi_final_label", "scanvi_confidence"],
    frameon=False,
    legend_loc="right margin",
    save=f"_{suffix}_scanvi_final_labels.png",
)

sc.pl.umap(
    adata,
    color=[SAMPLE_KEY, "scanvi_leiden_0.6"],
    frameon=False,
    save=f"_{suffix}_scanvi_sample_leiden.png",
)

print("\nFinal labels by species:")
final_label_table = pd.crosstab(
    adata.obs["scanvi_final_label"],
    adata.obs[SPECIES_KEY],
    margins=True,
)
print(final_label_table)
final_label_table.to_csv(results_path / "scanvi_final_label_by_species.csv")


# -----------------------------
# Species neighbor mixing on scANVI latent space
# -----------------------------

nn = NearestNeighbors(n_neighbors=N_NEIGHBORS + 1, metric="cosine")
nn.fit(adata.obsm[SCANVI_LATENT_KEY])

idx = nn.kneighbors(return_distance=False)
species = adata.obs[SPECIES_KEY].values

species_mixing = []
for i in range(adata.n_obs):
    neighbors = idx[i, 1:]
    frac_other = np.mean(species[neighbors] != species[i])
    species_mixing.append(frac_other)

adata.obs["species_mixing"] = species_mixing

mixing_summary = (
    adata.obs
    .groupby(["scanvi_final_label", SPECIES_KEY], observed=False)["species_mixing"]
    .agg(["mean", "median", "count"])
)

print("\nSpecies neighbor mixing:")
print(mixing_summary)
mixing_summary.to_csv(results_path / "scanvi_species_neighbor_mixing_summary.csv")


# -----------------------------
# Marker testing: all cells
# -----------------------------

sc.tl.rank_genes_groups(
    adata,
    groupby="scanvi_final_label",
    method="wilcoxon",
    use_raw=False,
)

markers_all = sc.get.rank_genes_groups_df(adata, group=None)
markers_all.to_csv(
    results_path / "markers_scanvi_final_mT_labels_all_cells.csv",
    index=False,
)


# -----------------------------
# Marker testing: mouse cells only
# -----------------------------

mouse = adata[mouse_mask].copy()

sc.tl.rank_genes_groups(
    mouse,
    groupby="scanvi_final_label",
    method="wilcoxon",
    use_raw=False,
)

markers_mouse = sc.get.rank_genes_groups_df(mouse, group=None)
markers_mouse.to_csv(
    results_path / "mouse_markers_scanvi_final_mT_labels.csv",
    index=False,
)


# -----------------------------
# Marker testing: human cells only
# -----------------------------

human = adata[human_mask].copy()

sc.tl.rank_genes_groups(
    human,
    groupby="scanvi_final_label",
    method="wilcoxon",
    use_raw=False,
)

markers_human = sc.get.rank_genes_groups_df(human, group=None)
markers_human.to_csv(
    results_path / "human_markers_scanvi_final_mT_labels.csv",
    index=False,
)


# -----------------------------
# Export metadata, latent space, AnnData, and models
# -----------------------------

metadata_cols = [
    SPECIES_KEY,
    SAMPLE_KEY,
    "reference_input_label",
    "scanvi_model_prediction",
    "scanvi_final_label",
    "scanvi_confidence",
    "species_mixing",
]
metadata_cols = [col for col in metadata_cols if col in adata.obs.columns]

adata.obs[metadata_cols].to_csv(results_path / "scanvi_metadata_full.csv")

scanvi_latent = pd.DataFrame(
    adata.obsm[SCANVI_LATENT_KEY],
    index=adata.obs_names,
    columns=[f"scANVI_{i + 1}" for i in range(adata.obsm[SCANVI_LATENT_KEY].shape[1])],
)
scanvi_latent.to_csv(results_path / "scanvi_latent_full.txt", sep="\t")

adata.write_h5ad(results_path / "tany_mouse_human_scanvi_final.h5ad")

vae.save(results_path / "scvi_model", overwrite=True)
scanvi_model.save(results_path / "scanvi_model", overwrite=True)

gc.collect()
print("\nDone.")
