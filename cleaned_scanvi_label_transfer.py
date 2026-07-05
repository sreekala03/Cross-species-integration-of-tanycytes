#!/usr/bin/env python
# coding: utf-8

"""
Cross-species tanycyte scVI/scANVI integration and label transfer.

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


# -----------------------------
# Helper functions
# -----------------------------

def load_parameters(param_file: Path) -> dict:
    """Load analysis parameters from JSON."""
    if not param_file.exists():
        raise FileNotFoundError(f"Parameter file not found: {param_file}")

    with open(param_file, "r") as f:
        params = json.load(f)

    required = [
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

    missing = [x for x in required if x not in params]
    if missing:
        raise ValueError(f"Missing parameters in JSON: {missing}")

    return params


def prepare_output_paths(params: dict) -> tuple[Path, Path, Path, str, int]:
    """Create output directory and resolve input/output paths."""
    results_path = Path(params["harmonization_folder_path"]).expanduser()
    results_path.mkdir(parents=True, exist_ok=True)

    data_file = Path(params["data_filepath_full"]).expanduser()
    feature_file = Path(params["feature_set_file"]).expanduser()
    suffix = params["new_name_suffix"]
    seed = int(params["global_seed"])

    return results_path, data_file, feature_file, suffix, seed


def read_and_subset_adata(data_file: Path, feature_file: Path) -> sc.AnnData:
    """Read AnnData object and subset to selected features."""
    adata = sc.read_h5ad(data_file)

    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.var_names_make_unique()

    with open(feature_file, "r") as f:
        feature_set = json.load(f)

    feature_set = [gene for gene in feature_set if gene in adata.var_names]
    print(f"Feature-set genes found in AnnData: {len(feature_set)}")

    if not feature_set:
        raise ValueError("No feature-set genes found in adata.var_names.")

    return adata[:, feature_set].copy()


def add_scvi_counts_layer(adata: sc.AnnData, layer_name: str = "scvi_counts") -> sc.AnnData:
    """Create a counts layer for scvi-tools."""
    print("Available layers before counts setup:", list(adata.layers.keys()))
    print("Has raw:", adata.raw is not None)

    if "counts" in adata.layers:
        adata.layers[layer_name] = adata.layers["counts"].copy()
    elif adata.raw is not None:
        raw_adata = adata.raw.to_adata()
        raw_adata = raw_adata[:, adata.var_names].copy()
        adata.layers[layer_name] = raw_adata.X.copy()
    else:
        adata.layers[layer_name] = adata.X.copy()

    print("Available layers after counts setup:", list(adata.layers.keys()))
    return adata


def check_metadata(adata: sc.AnnData, required_columns: list[str]) -> None:
    """Confirm required metadata columns are present."""
    missing = [col for col in required_columns if col not in adata.obs.columns]
    if missing:
        raise ValueError(f"Missing metadata columns: {missing}")

    for col in required_columns:
        print(f"\n{col}:")
        print(adata.obs[col].value_counts(dropna=False))


def save_latent(adata: sc.AnnData, latent_key: str, outfile: Path, prefix: str) -> None:
    """Save latent embedding from adata.obsm as a tab-delimited file."""
    latent = pd.DataFrame(
        adata.obsm[latent_key],
        index=adata.obs_names,
        columns=[f"{prefix}_{i + 1}" for i in range(adata.obsm[latent_key].shape[1])],
    )
    latent.to_csv(outfile, sep="\t")
    print(f"Saved latent space: {outfile}")


def run_neighbors_umap_leiden(
    adata: sc.AnnData,
    latent_key: str,
    seed: int,
    leiden_prefix: str,
    n_neighbors: int = 30,
    metric: str = "cosine",
    resolutions: list[float] | None = None,
) -> None:
    """Run neighbors, UMAP, and Leiden clustering on a latent representation."""
    if resolutions is None:
        resolutions = LEIDEN_RESOLUTIONS

    sc.pp.neighbors(
        adata,
        use_rep=latent_key,
        n_neighbors=n_neighbors,
        metric=metric,
    )

    sc.tl.umap(adata, random_state=seed)

    for res in resolutions:
        sc.tl.leiden(
            adata,
            resolution=res,
            key_added=f"{leiden_prefix}_{res}",
            flavor="igraph",
            directed=False,
            n_iterations=2,
        )


def summarize_species_overlap(
    adata: sc.AnnData,
    cluster_keys: list[str],
    species_key: str = SPECIES_KEY,
) -> None:
    """Print species composition per cluster."""
    for key in cluster_keys:
        print(f"\nSpecies overlap for {key}")
        print(pd.crosstab(adata.obs[key], adata.obs[species_key], normalize="index"))


def calculate_species_neighbor_mixing(
    adata: sc.AnnData,
    latent_key: str = SCANVI_LATENT_KEY,
    species_key: str = SPECIES_KEY,
    label_key: str = "scanvi_final_label",
    k: int = 30,
) -> pd.DataFrame:
    """Calculate fraction of nearest neighbors from the other species."""
    nn = NearestNeighbors(n_neighbors=k + 1, metric="cosine")
    nn.fit(adata.obsm[latent_key])

    idx = nn.kneighbors(return_distance=False)
    species = adata.obs[species_key].values

    adata.obs["species_mixing"] = [
        np.mean(species[neighbors[1:]] != species[i])
        for i, neighbors in enumerate(idx)
    ]

    summary = (
        adata.obs
        .groupby([label_key, species_key], observed=False)["species_mixing"]
        .agg(["mean", "median", "count"])
    )

    return summary


def save_markers(
    adata: sc.AnnData,
    groupby: str,
    outfile: Path,
    subset_query: str | None = None,
) -> None:
    """Run Wilcoxon marker testing and save ranked genes."""
    test_adata = adata.copy()
    if subset_query is not None:
        test_adata = test_adata[test_adata.obs.eval(subset_query)].copy()

    sc.tl.rank_genes_groups(
        test_adata,
        groupby=groupby,
        method="wilcoxon",
        use_raw=False,
    )

    markers = sc.get.rank_genes_groups_df(test_adata, group=None)
    markers.to_csv(outfile, index=False)
    print(f"Saved markers: {outfile}")


# -----------------------------
# Main workflow
# -----------------------------

def main() -> None:
    print("scvi-tools version:", scvi.__version__)

    params = load_parameters(PARAM_FILE)
    results_path, data_file, feature_file, suffix, seed = prepare_output_paths(params)

    scvi.settings.seed = seed
    np.random.seed(seed)

    print("\nInput/output paths")
    print("AnnData:", data_file)
    print("Feature set:", feature_file)
    print("Results:", results_path)

    adata = read_and_subset_adata(data_file, feature_file)
    adata = add_scvi_counts_layer(adata)

    check_metadata(
        adata,
        required_columns=[
            SPECIES_KEY,
            SAMPLE_KEY,
            LABEL_KEY,
        ],
    )

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

    save_latent(
        adata,
        latent_key=SCVI_LATENT_KEY,
        outfile=results_path / f"{suffix}_scVI_reduction.txt",
        prefix="scVI",
    )

    run_neighbors_umap_leiden(
        adata,
        latent_key=SCVI_LATENT_KEY,
        seed=seed,
        leiden_prefix="scvi_leiden",
        n_neighbors=int(params["k_param"]),
        metric=params.get("dist_type", "cosine"),
    )

    summarize_species_overlap(
        adata,
        cluster_keys=[f"scvi_leiden_{res}" for res in LEIDEN_RESOLUTIONS],
    )

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

    # Predict all cells for diagnostics only.
    adata.obs["reference_input_label"] = adata.obs[LABEL_KEY].astype(str)
    adata.obs["scanvi_model_prediction"] = scanvi_model.predict(adata)

    pred_probs = scanvi_model.predict(adata, soft=True)
    adata.obs["scanvi_confidence"] = pred_probs.max(axis=1)

    # -----------------------------
    # scANVI label transfer while preserving mouse reference labels
    # -----------------------------

    mouse_mask = adata.obs[SPECIES_KEY].astype(str).eq("mouse")
    human_mask = adata.obs[SPECIES_KEY].astype(str).eq("human")

    adata.obs["scanvi_final_label"] = adata.obs["reference_input_label"].copy()

    # Preserve curated mouse labels; assign predictions only to human cells.
    adata.obs.loc[human_mask, "scanvi_final_label"] = scanvi_model.predict(
        adata[human_mask].copy()
    )

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
    # scANVI embedding, clustering, plots, and QC
    # -----------------------------

    run_neighbors_umap_leiden(
        adata,
        latent_key=SCANVI_LATENT_KEY,
        seed=seed,
        leiden_prefix="scanvi_leiden",
        n_neighbors=N_NEIGHBORS,
        metric="cosine",
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

    mixing_summary = calculate_species_neighbor_mixing(
        adata,
        latent_key=SCANVI_LATENT_KEY,
        species_key=SPECIES_KEY,
        label_key="scanvi_final_label",
        k=N_NEIGHBORS,
    )
    print("\nSpecies neighbor mixing:")
    print(mixing_summary)
    mixing_summary.to_csv(results_path / "scanvi_species_neighbor_mixing_summary.csv")

    # -----------------------------
    # Marker testing
    # -----------------------------

    save_markers(
        adata,
        groupby="scanvi_final_label",
        outfile=results_path / "markers_scanvi_final_mT_labels_all_cells.csv",
    )

    save_markers(
        adata,
        groupby="scanvi_final_label",
        outfile=results_path / "mouse_markers_scanvi_final_mT_labels.csv",
        subset_query=f"{SPECIES_KEY} == 'mouse'",
    )

    save_markers(
        adata,
        groupby="scanvi_final_label",
        outfile=results_path / "human_markers_scanvi_final_mT_labels.csv",
        subset_query=f"{SPECIES_KEY} == 'human'",
    )

    # -----------------------------
    # Export final objects
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

    save_latent(
        adata,
        latent_key=SCANVI_LATENT_KEY,
        outfile=results_path / "scanvi_latent_full.txt",
        prefix="scANVI",
    )

    adata.write_h5ad(results_path / "tany_mouse_human_scanvi_final.h5ad")

    vae.save(results_path / "scvi_model", overwrite=True)
    scanvi_model.save(results_path / "scanvi_model", overwrite=True)

    gc.collect()
    print("\nDone.")


if __name__ == "__main__":
    main()
