import sys
import os

import anndata as ad

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.simulate_immune_atlas import simulate
from src.reference_annotation import run_annotation_pipeline
from src.differential_abundance import run_differential_abundance_test, validate_disease_state_recovery


def _build_adata():
    counts, gene_names, obs, marker_genes = simulate()
    adata = ad.AnnData(X=counts.astype("float32"), obs=obs.copy())
    adata.var_names = gene_names
    return adata, marker_genes


def test_atlas_has_80000plus_cells_and_5_tissues():
    counts, gene_names, obs, marker_genes = simulate()
    assert counts.shape[0] >= 80000
    assert obs["tissue"].nunique() == 5


def test_reference_annotation_is_accurate():
    adata, _ = _build_adata()
    annotation = run_annotation_pipeline(adata)
    assert annotation["accuracy_pct"] >= 70.0


def test_differential_abundance_recovers_most_disease_states():
    adata, _ = _build_adata()
    run_annotation_pipeline(adata)
    da_results = run_differential_abundance_test(adata.obs)
    recovery = validate_disease_state_recovery(da_results)
    assert recovery["recall_pct"] >= 60.0


def test_fdr_correction_applied():
    adata, _ = _build_adata()
    run_annotation_pipeline(adata)
    da_results = run_differential_abundance_test(adata.obs)
    assert "fdr" in da_results.columns
    assert (da_results["fdr"] >= da_results["p_value"]).all()


def test_significant_states_have_fdr_below_threshold():
    adata, _ = _build_adata()
    run_annotation_pipeline(adata)
    da_results = run_differential_abundance_test(adata.obs)
    sig = da_results[da_results["significant"]]
    assert (sig["fdr"] < 0.05).all()
