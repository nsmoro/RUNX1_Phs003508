#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import types
import platform
import pickle

import numpy as np
import pandas as pd

# -----------------------------------------------------------------------------
# Shims so a model pickled with scikit-learn 0.21.x can be loaded on 1.x
# -----------------------------------------------------------------------------
# Map 'sklearn.preprocessing.data' -> 'sklearn.preprocessing'
import sklearn.preprocessing as _sk_pre
_pre_shim = types.ModuleType("sklearn.preprocessing.data")
_pre_shim.__dict__.update(_sk_pre.__dict__)
sys.modules["sklearn.preprocessing.data"] = _pre_shim

# Map 'sklearn.neural_network.multilayer_perceptron' -> current MLP module
try:
    import sklearn.neural_network._multilayer_perceptron as _mlp_mod
except Exception:
    import sklearn.neural_network as _mlp_mod
_mlp_shim = types.ModuleType("sklearn.neural_network.multilayer_perceptron")
_mlp_shim.__dict__.update(_mlp_mod.__dict__)
sys.modules["sklearn.neural_network.multilayer_perceptron"] = _mlp_shim

# Some very old pickles referenced sklearn.externals.joblib
try:
    import joblib as _jl
    sys.modules["sklearn.externals.joblib"] = _jl
except Exception:
    pass
# -----------------------------------------------------------------------------

import sklearn
from sklearn.decomposition import PCA

import matplotlib
matplotlib.use("Agg")  # headless (HPC)
import matplotlib.pyplot as plt


def total_count_normalise(count_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Normalise count matrix for input into hscScore model.
    Expects cells x genes (rows=cells, cols=genes). Returns log1p-normalised DataFrame.
    """
    if not isinstance(count_matrix, pd.DataFrame):
        raise ValueError("count_matrix must be a pandas DataFrame")

    if count_matrix.isnull().values.any() or np.isinf(count_matrix.values).any():
        raise ValueError("count_matrix contains NaN or Inf values")

    wilson_molo_genes_median_counts = 18704.5

    # Convert to float numpy array
    X = count_matrix.astype(float).to_numpy(copy=False)

    # Per-cell depth
    counts_per_cell = X.sum(axis=1)
    counts_per_cell = counts_per_cell + (counts_per_cell == 0)  # avoid /0
    counts_per_cell = counts_per_cell / wilson_molo_genes_median_counts

    # Depth normalise
    norm = X / counts_per_cell[:, None]

    # Back to DataFrame and log1p
    norm_df = pd.DataFrame(norm, index=count_matrix.index, columns=count_matrix.columns)
    norm_df = np.log1p(norm_df)

    return norm_df


def main():
    print("Python:", platform.python_version())
    print("scikit-learn:", sklearn.__version__)

    # ---- Paths ----
    base_path = "/home/"
    tsv_path = os.path.join(base_path, "all_GSE168807.tsv")
    pkl_path = os.path.join(base_path, "hscScore_model.pkl")
    model_genes_path = os.path.join(base_path, "model_molo_genes.txt")

    score_csv_out = "/home/score.csv"
    plot_out = os.path.join(os.getcwd(), "score_plot.png")

    # ---- Load counts (genes as columns, cells as rows expected by pipeline) ----
    print(f"Reading counts from: {tsv_path}")
    df = pd.read_csv(tsv_path, sep="\t", index_col=0)
    print("Raw df (first 5 rows, 3 cols):")
    print(df.iloc[:5, :3])

    # If df is genes x cells, transpose to cells x genes
    # Your previous workflow used transpose, so we keep it explicit.
    count_matrix = df.transpose()
    print("count_matrix head (first 5 rows, 3 cols):")
    print(count_matrix.iloc[:5, :3])

    # ---- Load model safely (shims already set above) ----
    print(f"Loading model from: {pkl_path}")
    with open(pkl_path, "rb") as f:
        hsc_score = pickle.load(f)
    print("Model loaded.")

    # ---- Ensure gene order matches model expectations, fill missing with 0 ----
    print(f"Loading model gene list from: {model_genes_path}")
    model_genes = np.genfromtxt(model_genes_path, dtype=str)
    model_genes = list(model_genes)

    missing = [g for g in model_genes if g not in count_matrix.columns]
    if missing:
        print(f"[Info] {len(missing)} model genes missing in matrix; filling with 0.")

    count_data_molo = count_matrix.reindex(columns=model_genes, fill_value=0.0)

    # ---- Normalise ----
    print("Normalising data...")
    normalised_data_molo = total_count_normalise(count_data_molo)
    print("normalised_data_molo head:")
    print(normalised_data_molo.iloc[:5, :5])

    # ---- Predict ----
    print("Predicting HSC scores...")
    predicted_hsc_scores = hsc_score.predict(normalised_data_molo.to_numpy())
    predicted_hsc_scores_df = pd.DataFrame(
        predicted_hsc_scores, index=normalised_data_molo.index, columns=["HSC_Score"]
    )

    # ---- Save scores ----
    predicted_hsc_scores_df.to_csv(score_csv_out, index=True)  # keep cell IDs
    print(f"Saved scores -> {score_csv_out}")

    # ---- PCA + plot (color by HSC score) ----
    print("Running PCA and saving plot...")
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(normalised_data_molo.values)
    pca_df = pd.DataFrame(X_pca, index=normalised_data_molo.index, columns=["PC1", "PC2"])
    pca_df["HSC_Score"] = predicted_hsc_scores_df["HSC_Score"].values

    plt.figure()
    sc = plt.scatter(pca_df["PC1"], pca_df["PC2"], c=pca_df["HSC_Score"], cmap="plasma")
    plt.title("HSC-score")
    plt.colorbar(sc, label="HSC Score")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.savefig(plot_out, dpi=200, bbox_inches="tight")
    print(f"Saved plot -> {plot_out}")

    print("Done.")


if __name__ == "__main__":
    main()

