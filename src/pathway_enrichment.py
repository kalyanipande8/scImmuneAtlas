"""Real hypergeometric pathway enrichment on each disease-associated
cell state's marker genes, against a small simulated immune-pathway
gene-set catalog."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

N_PATHWAYS = 6


def build_gene_sets(gene_names: list, seed: int = 839) -> dict:
    rng = np.random.default_rng(seed)
    return {f"Pathway_{p}": set(rng.choice(gene_names, rng.integers(6, 14), replace=False))
            for p in range(N_PATHWAYS)}


def enrich_state_markers(marker_genes: dict, gene_names: list, gene_sets: dict, alpha: float = 0.05) -> pd.DataFrame:
    results = []
    universe_size = len(gene_names)
    for state, markers in marker_genes.items():
        for pathway, genes in gene_sets.items():
            overlap = len(set(markers) & genes)
            p_value = hypergeom.sf(overlap - 1, universe_size, len(genes), len(markers))
            results.append({"state": state, "pathway": pathway, "overlap": overlap, "p_value": p_value})

    df = pd.DataFrame(results)
    df["significant"] = df["p_value"] < alpha
    return df
