"""Fold STRING physical/functional scores into the ranking (top-N re-rank).

STRING is one API call per gene, so we only score the top-N genes (by selectivity
score) and re-rank that subset:

    string_score = selectivity_score * (1 + w_phys*phys_max + w_func*func_max)

with separate physical and functional weights. Also builds candidate x seed score
matrices (for the heatmap) and an edge list (for the network graph).
Cached; graceful-empty offline.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from agents import string_db


def string_rerank(ranked: pd.DataFrame, seeds: tuple[str, ...], top_n: int = 25,
                  w_sel: float = 0.5, w_phys: float = 0.35, w_func: float = 0.15,
                  cutoff: float = 0.4
                  ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict]]:
    """Return (rerank_df, phys_matrix, func_matrix, edges).

    Normalized weighted-sum score (all terms in [0,1]):
        string_score = w_sel*sel_norm + w_phys*physical_max + w_func*functional_max
    where sel_norm is selectivity_score min-max scaled across the top-N, and the
    three weights are auto-normalized to sum to 1 so string_score stays in [0,1].

    rerank_df: top-N genes with sel_norm, string_phys, string_func, string_score,
      sorted by string_score desc. phys_matrix/func_matrix: candidate x seed STRING
      scores. edges: [{source, target, physical, functional}] for the network graph.
    """
    if ranked is None or ranked.empty:
        empty = pd.DataFrame()
        return empty, empty, empty, []
    top = ranked.head(top_n).copy()
    genes = list(top.index)

    phys_rows, func_rows, edges = {}, {}, []
    smax_p, smax_f, best_p, best_f = [], [], [], []
    sc_p, sc_f = [], []   # cutoff-thresholded maxima used for scoring
    for g in genes:
        sd = string_db.string_scores(g, seeds, cutoff)
        pl = sd.get("phys_links", {})
        fl = sd.get("func_links", {})
        phys_rows[g] = {s: pl.get(s, 0.0) for s in seeds}
        func_rows[g] = {s: fl.get(s, 0.0) for s in seeds}
        smax_p.append(sd["physical"]["max"])
        smax_f.append(sd["functional"]["max"])
        sc_p.append(sd["physical"]["max_thresh"])
        sc_f.append(sd["functional"]["max_thresh"])
        best_p.append(sd["physical"]["best_seed"])
        best_f.append(sd["functional"]["best_seed"])
        for s in seeds:
            pv, fv = pl.get(s, 0.0), fl.get(s, 0.0)
            if pv > 0 or fv > 0:
                edges.append({"source": g, "target": s,
                              "physical": round(pv, 3), "functional": round(fv, 3)})

    # min-max normalize selectivity across the top-N
    sel = top["selectivity_score"].to_numpy(dtype=float)
    lo, hi = float(sel.min()), float(sel.max())
    sel_norm = np.ones_like(sel) if hi <= lo else (sel - lo) / (hi - lo)

    # auto-normalize weights to sum to 1 (fallback to selectivity-only if all zero)
    wsum = w_sel + w_phys + w_func
    if wsum <= 0:
        w_sel, w_phys, w_func, wsum = 1.0, 0.0, 0.0, 1.0
    n_sel, n_phys, n_func = w_sel / wsum, w_phys / wsum, w_func / wsum

    # scoring uses cutoff-thresholded STRING maxima (weak links contribute 0)
    sc_p_a, sc_f_a = np.array(sc_p), np.array(sc_f)
    score = n_sel * sel_norm + n_phys * sc_p_a + n_func * sc_f_a

    top["sel_norm"] = sel_norm.round(3)
    top["string_phys"] = smax_p
    top["string_func"] = smax_f
    top["string_best_phys_seed"] = best_p
    top["string_best_func_seed"] = best_f
    top["string_score"] = score.round(4)
    top.attrs["weights"] = {"selectivity": round(n_sel, 3), "physical": round(n_phys, 3),
                            "functional": round(n_func, 3)}
    top = top.sort_values("string_score", ascending=False)

    phys_matrix = pd.DataFrame(phys_rows).T.reindex(index=top.index)
    func_matrix = pd.DataFrame(func_rows).T.reindex(index=top.index)
    return top, phys_matrix, func_matrix, edges
