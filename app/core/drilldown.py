"""Per-cell-line vectors for the drill-down scatter (candidate vs seed).

Kept separate from the agent pipeline so the UI can fetch two gene vectors
cheaply without recomputing correlations. Reads harmonized parquets only.
"""
from __future__ import annotations

import pandas as pd

from agents.data_agent import _load_all
from core.cohort import CohortSpec, select_models


def cohort_gene_vectors(
    spec: CohortSpec, gene_a: str, gene_b: str, modality: str
) -> pd.DataFrame | None:
    """Aligned (x=gene_a, y=gene_b) across the cohort's cell lines.

    modality: "expr" or "effect". Returns a DataFrame with columns
    [x, y, cell_line, disease] indexed by ModelID, or None if a gene is absent.
    """
    expr, effect, meta = _load_all()
    mat = expr if modality == "expr" else effect
    if gene_a not in mat.columns or gene_b not in mat.columns:
        return None
    models = select_models(meta, spec)
    both = [m for m in models if m in mat.index]
    sub = mat.loc[both, [gene_a, gene_b]].dropna()
    if sub.empty:
        return None
    labels = meta.loc[sub.index, "CellLineName"].fillna(sub.index.to_series())
    disease = meta.loc[sub.index, "OncotreePrimaryDisease"].fillna("")
    return pd.DataFrame(
        {"x": sub[gene_a].to_numpy(), "y": sub[gene_b].to_numpy(),
         "cell_line": labels.to_numpy(), "disease": disease.to_numpy()},
        index=sub.index,
    )
