import numpy as np
import pandas as pd
import scanpy as sc
import scvi

from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    classification_report,
    confusion_matrix
)

# -----------------------------
# Settings
# -----------------------------

seed = 1234

input_h5ad = Path(
    "~/snRNA_seq_integ/cross_species_tanycytes_full/Tanycytes_mouse_RNA_counts_only_scanvi_latent.h5ad"
).expanduser()

true_label_col = "Tany.cell.state"
batch_col = "Sample.name"   # set to None if no batch correction is needed

out_dir = Path(
    "~/snRNA_seq_integ/cross_species_tanycytes_full/mouse_scanvi_holdout_validation"
).expanduser()

out_dir.mkdir(parents=True, exist_ok=True)

# -----------------------------
# Load data
# -----------------------------

adata = sc.read_h5ad(input_h5ad)

if true_label_col not in adata.obs.columns:
    raise ValueError(f"Missing label column: {true_label_col}")

if batch_col is not None and batch_col not in adata.obs.columns:
    raise ValueError(f"Missing batch column: {batch_col}")

adata.obs["true_label"] = adata.obs[true_label_col].astype(str)

print(adata)
print(adata.obs["true_label"].value_counts())

# -----------------------------
# Create stratified train/test split
# -----------------------------

train_cells, test_cells = train_test_split(
    adata.obs_names,
    test_size=0.2,
    random_state=seed,
    stratify=adata.obs["true_label"]
)

adata.obs["scanvi_train_label"] = adata.obs["true_label"]
adata.obs.loc[test_cells, "scanvi_train_label"] = "Unknown"

print(adata.obs["scanvi_train_label"].value_counts())

# -----------------------------
# Train scVI
# -----------------------------

setup_kwargs = {
    "labels_key": "scanvi_train_label"
}

if batch_col is not None:
    setup_kwargs["batch_key"] = batch_col

scvi.model.SCVI.setup_anndata(
    adata,
    **setup_kwargs
)

vae = scvi.model.SCVI(
    adata,
    n_latent=30,
    gene_likelihood="nb"
)

vae.train(
    max_epochs=150,
    early_stopping=True
)

# -----------------------------
# Train scANVI
# -----------------------------

scanvi = scvi.model.SCANVI.from_scvi_model(
    vae,
    unlabeled_category="Unknown"
)

scanvi.train(
    max_epochs=100,
    early_stopping=True
)

# -----------------------------
# Predict held-out labels
# -----------------------------

adata.obs["scanvi_predicted_label"] = scanvi.predict(adata)

test_df = adata.obs.loc[test_cells].copy()

y_true = test_df["true_label"]
y_pred = test_df["scanvi_predicted_label"]

# -----------------------------
# Accuracy metrics
# -----------------------------

overall_accuracy = accuracy_score(y_true, y_pred)
balanced_accuracy = balanced_accuracy_score(y_true, y_pred)

metrics_summary = pd.DataFrame(
    {
        "metric": ["overall_accuracy", "balanced_accuracy"],
        "value": [overall_accuracy, balanced_accuracy]
    }
)

metrics_summary.to_csv(
    out_dir / "mouse_scanvi_holdout_global_accuracy.csv",
    index=False
)

print(metrics_summary)

report = classification_report(
    y_true,
    y_pred,
    digits=3,
    output_dict=True
)

report_df = pd.DataFrame(report).transpose()

report_df.to_csv(
    out_dir / "mouse_scanvi_holdout_classification_report.csv"
)

print(classification_report(y_true, y_pred, digits=3))

# -----------------------------
# Confusion matrix
# -----------------------------

labels_order = sorted(y_true.unique())

conf_mat = pd.DataFrame(
    confusion_matrix(
        y_true,
        y_pred,
        labels=labels_order
    ),
    index=labels_order,
    columns=labels_order
)

conf_mat.to_csv(
    out_dir / "mouse_scanvi_holdout_confusion_matrix.csv"
)

print(conf_mat)

# -----------------------------
# Per-label accuracy
# -----------------------------

accuracy_by_label = (
    test_df
    .assign(correct=y_true.values == y_pred.values)
    .groupby("true_label")
    .agg(
        n_cells=("correct", "size"),
        correct=("correct", "sum"),
        accuracy=("correct", "mean")
    )
    .reset_index()
    .sort_values("accuracy", ascending=False)
)

accuracy_by_label.to_csv(
    out_dir / "mouse_scanvi_holdout_accuracy_by_label.csv",
    index=False
)

print(accuracy_by_label)

# -----------------------------
# Optional UMAP using scANVI latent space
# -----------------------------

adata.obsm["X_scanvi"] = scanvi.get_latent_representation(adata)

sc.pp.neighbors(
    adata,
    use_rep="X_scanvi",
    n_neighbors=30
)

sc.tl.umap(
    adata,
    random_state=seed
)

# -----------------------------
# Save AnnData
# -----------------------------

adata.write_h5ad(
    out_dir / "mouse_scanvi_holdout_validation_predictions.h5ad"
)

print("Done.")
print(f"Results saved to: {out_dir}")
