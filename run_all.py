import json
import os

import anndata as ad

from src.simulate_immune_atlas import simulate
from src.reference_annotation import run_annotation_pipeline
from src.differential_abundance import run_differential_abundance_test, validate_disease_state_recovery
from src.pathway_enrichment import build_gene_sets, enrich_state_markers


def main():
    counts, gene_names, obs, marker_genes = simulate()
    adata = ad.AnnData(X=counts.astype("float32"), obs=obs.copy())
    adata.var_names = gene_names

    annotation = run_annotation_pipeline(adata)

    da_results = run_differential_abundance_test(adata.obs)
    da_recovery = validate_disease_state_recovery(da_results)

    disease_states = da_results.loc[da_results["significant"], "state"].tolist()
    disease_markers = {s: marker_genes[s] for s in disease_states if s in marker_genes}
    gene_sets = build_gene_sets(gene_names)
    enrichment = enrich_state_markers(disease_markers, gene_names, gene_sets)
    n_enriched = int(enrichment["significant"].sum())

    summary = {"n_cells": counts.shape[0], "n_tissues": obs["tissue"].nunique(),
               "annotation": annotation, "differential_abundance_recovery": da_recovery,
               "n_significant_pathway_hits": n_enriched}

    os.makedirs("results", exist_ok=True)
    with open("results/summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
