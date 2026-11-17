"""Real differential-abundance testing: a chi-square test of cell-
state proportion between disease and healthy conditions, per state,
with FDR (Benjamini-Hochberg) correction across all 20 states --
validated against the simulation's known true disease-associated
states."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests

from .simulate_immune_atlas import DISEASE_ASSOCIATED_STATES


def run_differential_abundance_test(adata_obs: pd.DataFrame, state_col: str = "predicted_state") -> pd.DataFrame:
    results = []
    all_states = adata_obs[state_col].unique()
    total_disease = (adata_obs["condition"] == "disease").sum()
    total_healthy = (adata_obs["condition"] == "healthy").sum()

    for state in all_states:
        n_disease_state = ((adata_obs[state_col] == state) & (adata_obs["condition"] == "disease")).sum()
        n_healthy_state = ((adata_obs[state_col] == state) & (adata_obs["condition"] == "healthy")).sum()

        table = [[n_disease_state, total_disease - n_disease_state],
                  [n_healthy_state, total_healthy - n_healthy_state]]
        chi2, p_value, _, _ = chi2_contingency(table)
        results.append({"state": state, "n_disease": n_disease_state, "n_healthy": n_healthy_state,
                         "p_value": p_value})

    df = pd.DataFrame(results)
    df["fdr"] = multipletests(df["p_value"], method="fdr_bh")[1]
    df["significant"] = df["fdr"] < 0.05
    return df


def validate_disease_state_recovery(da_results: pd.DataFrame) -> dict:
    detected = set(da_results.loc[da_results["significant"], "state"])
    true_set = set(DISEASE_ASSOCIATED_STATES)
    tp = len(detected & true_set)
    return {"n_true_disease_states": len(true_set), "n_detected": len(detected), "n_true_positives": tp,
            "recall_pct": round(100 * tp / len(true_set), 2),
            "precision_pct": round(100 * tp / len(detected), 2) if detected else 0.0}
