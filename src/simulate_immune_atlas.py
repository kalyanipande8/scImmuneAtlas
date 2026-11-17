"""Simulate an 82,000-cell immune atlas across 5 tissue sources and 2
conditions (disease/healthy), with 20 known true cell states, 14 of
which are truly disease-associated (real differential abundance in
disease vs healthy) -- the ground truth a reference-annotation +
differential-abundance pipeline must recover."""
from __future__ import annotations

import numpy as np
import pandas as pd

TISSUES = ["blood", "lymph_node", "spleen", "bone_marrow", "tumor"]
N_CELLS_PER_TISSUE = 16400  # 5 x 16400 = 82,000
N_GENES = 80
N_CELL_STATES = 20
N_MARKERS_PER_STATE = 4
SEED = 839

DISEASE_ASSOCIATED_STATES = [f"State_{i:02d}" for i in range(14)]  # 14 of 20 truly disease-associated


def simulate(n_cells_per_tissue: int = N_CELLS_PER_TISSUE, seed: int = SEED):
    rng = np.random.default_rng(seed)
    gene_names = [f"GENE_{i:03d}" for i in range(N_GENES)]
    state_names = [f"State_{i:02d}" for i in range(N_CELL_STATES)]

    marker_genes = {state: gene_names[i * N_MARKERS_PER_STATE:(i + 1) * N_MARKERS_PER_STATE]
                     for i, state in enumerate(state_names)}

    # baseline proportions (healthy), and a real shift for disease-associated states
    baseline_props = rng.dirichlet(np.ones(N_CELL_STATES) * 3)
    disease_props = baseline_props.copy()
    for i, state in enumerate(state_names):
        if state in DISEASE_ASSOCIATED_STATES:
            disease_props[i] *= rng.choice([0.3, 3.0])
    disease_props /= disease_props.sum()

    all_counts, all_states, all_tissues, all_conditions = [], [], [], []

    for tissue in TISSUES:
        for condition, props in [("healthy", baseline_props), ("disease", disease_props)]:
            n_this = n_cells_per_tissue // 2
            states = rng.choice(state_names, n_this, p=props)

            base_expr = rng.gamma(2.0, 2.0, N_GENES)
            counts = rng.poisson(base_expr, size=(n_this, N_GENES)).astype(float)

            for state in state_names:
                mask = states == state
                for g in marker_genes[state]:
                    counts[mask, gene_names.index(g)] += rng.poisson(14, mask.sum())

            all_counts.append(counts)
            all_states.extend(states)
            all_tissues.extend([tissue] * n_this)
            all_conditions.extend([condition] * n_this)

    counts_matrix = np.vstack(all_counts)
    obs = pd.DataFrame({"true_state": all_states, "tissue": all_tissues, "condition": all_conditions})
    return counts_matrix, gene_names, obs, marker_genes


if __name__ == "__main__":
    counts, genes, obs, markers = simulate()
    print(counts.shape, obs["condition"].value_counts().to_dict())
