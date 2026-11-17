"""Real reference-based cell annotation: builds a reference centroid
profile per cell state from a labeled reference subset, then assigns
every other cell to its nearest centroid (correlation distance) --
the real nearest-centroid logic behind tools like scmap/SingleR,
implemented directly rather than via full scanpy clustering."""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad


def build_reference_centroids(adata: ad.AnnData, reference_frac: float = 0.2, seed: int = 839):
    rng = np.random.default_rng(seed)
    ref_idx = rng.choice(adata.n_obs, int(adata.n_obs * reference_frac), replace=False)
    ref_mask = np.zeros(adata.n_obs, dtype=bool)
    ref_mask[ref_idx] = True

    ref_adata = adata[ref_mask].copy()
    centroids = ref_adata.to_df().groupby(ref_adata.obs["true_state"].values).mean()
    return centroids, ~ref_mask


def _corr_to_centroids(query: np.ndarray, centroids: np.ndarray) -> np.ndarray:
    """Correlation between each query row and each centroid row,
    computed directly as a (n_query x n_centroids) matrix instead of
    the much larger (n_query+n_centroids)^2 full correlation matrix."""
    q = query - query.mean(axis=1, keepdims=True)
    c = centroids - centroids.mean(axis=1, keepdims=True)
    q_norm = q / np.linalg.norm(q, axis=1, keepdims=True)
    c_norm = c / np.linalg.norm(c, axis=1, keepdims=True)
    return q_norm @ c_norm.T


def annotate_by_nearest_centroid(adata: ad.AnnData, centroids: pd.DataFrame, query_mask) -> pd.Series:
    query_df = adata[query_mask].to_df()
    correlations = _corr_to_centroids(query_df.values, centroids.values)
    assigned_idx = np.argmax(correlations, axis=1)
    assigned_states = centroids.index[assigned_idx]
    return pd.Series(assigned_states, index=adata.obs_names[query_mask])


def run_annotation_pipeline(adata: ad.AnnData) -> dict:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    centroids, query_mask = build_reference_centroids(adata)
    predicted = annotate_by_nearest_centroid(adata, centroids, query_mask)

    true_labels = adata.obs.loc[query_mask, "true_state"]
    accuracy = (predicted.values == true_labels.values).mean()

    adata.obs["predicted_state"] = adata.obs["true_state"]
    adata.obs.loc[query_mask, "predicted_state"] = predicted.values

    return {"accuracy_pct": round(float(accuracy) * 100, 3), "n_reference_cells": int((~query_mask).sum()),
            "n_query_cells": int(query_mask.sum())}
