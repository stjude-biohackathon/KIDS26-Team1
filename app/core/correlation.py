"""Vectorized Pearson correlation of a candidate matrix against a seed matrix,
with per-pair p-values and Benjamini-Hochberg FDR.

All correlations are computed *across cell lines* (rows). NaNs are handled
pairwise where cheap; for the fast path we require complete columns (callers
impute or drop first). Designed for m ~ 18k candidates x k ~ 23 seeds.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy import stats


def _zscore(mat: np.ndarray) -> np.ndarray:
    """Column-wise z-score; constant columns become all-zero (r -> 0)."""
    mu = np.nanmean(mat, axis=0)
    sd = np.nanstd(mat, axis=0, ddof=1)
    sd_safe = np.where(sd == 0, np.nan, sd)
    z = (mat - mu) / sd_safe
    return z


def _rank_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Average-rank each column (Spearman = Pearson on ranks). Matches scipy."""
    return df.rank(axis=0, method="average")


def pearson_matrix(
    cand: pd.DataFrame, seeds: pd.DataFrame, method: str = "pearson"
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Correlate every candidate column against every seed column across rows.

    method: "pearson" (default) or "spearman" (rank-transform then Pearson).
    cand, seeds share the same row index (cell lines), already aligned.
    Returns (r, p): DataFrames of shape (n_candidates, n_seeds).
    Assumes no NaNs (caller cleans). n = number of rows.
    """
    idx = cand.index
    assert seeds.index.equals(idx), "cand and seeds must be row-aligned"
    n = len(idx)
    if n < 3:
        raise ValueError(f"need >=3 cell lines, got {n}")

    if method == "spearman":
        cand = _rank_columns(cand)
        seeds = _rank_columns(seeds)
    elif method != "pearson":
        raise ValueError(f"unknown method {method!r}")

    Zc = _zscore(cand.to_numpy(dtype=np.float64))  # n x m
    Zs = _zscore(seeds.to_numpy(dtype=np.float64))  # n x k

    # r = (Zc^T Zs) / (n-1); constant columns are NaN -> r NaN
    r = (Zc.T @ np.nan_to_num(Zs)) / (n - 1)
    # zero out where a seed column was constant
    seed_bad = np.isnan(Zs).all(axis=0)
    cand_bad = np.isnan(Zc).all(axis=0)
    r[cand_bad, :] = np.nan
    r[:, seed_bad] = np.nan
    r = np.clip(r, -1.0, 1.0)

    # two-sided p from t = r sqrt((n-2)/(1-r^2))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = r * np.sqrt((n - 2) / (1.0 - r**2))
    p = 2 * stats.t.sf(np.abs(t), df=n - 2)
    p = np.where(np.isnan(r), np.nan, p)

    r_df = pd.DataFrame(r, index=cand.columns, columns=seeds.columns)
    p_df = pd.DataFrame(p, index=cand.columns, columns=seeds.columns)
    return r_df, p_df


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values over all non-NaN entries (flattened)."""
    flat = pvals.ravel()
    mask = ~np.isnan(flat)
    q = np.full_like(flat, np.nan, dtype=np.float64)
    p = flat[mask]
    if p.size == 0:
        return q.reshape(pvals.shape)
    order = np.argsort(p)
    ranked = p[order]
    m = p.size
    q_sorted = ranked * m / (np.arange(1, m + 1))
    # enforce monotonicity from the largest p downward
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0, 1)
    q_vals = np.empty(m)
    q_vals[order] = q_sorted
    q[mask] = q_vals
    return q.reshape(pvals.shape)


def aggregate_per_candidate(
    r_df: pd.DataFrame, q_df: pd.DataFrame, q_thresh: float = 0.05
) -> pd.DataFrame:
    """Collapse the candidate x seed matrix to one row per candidate."""
    r = r_df.to_numpy()
    q = q_df.to_numpy()
    absr = np.abs(r)
    all_nan_rows = np.isnan(absr).all(axis=1)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        best_idx = np.nanargmax(np.where(np.isnan(absr), -1, absr), axis=1)
        mean_r = np.nanmean(r, axis=1)
        mean_abs_r = np.nanmean(absr, axis=1)
        max_abs_r = np.nanmax(np.where(np.isnan(absr), -np.inf, absr), axis=1)
        min_q = np.nanmin(np.where(np.isnan(q), np.inf, q), axis=1)
    max_abs_r = np.where(all_nan_rows, np.nan, max_abs_r)
    min_q = np.where(np.isinf(min_q), np.nan, min_q)
    seeds = np.array(r_df.columns)
    out = pd.DataFrame(
        {
            "gene": r_df.index,
            "mean_r": mean_r,
            "mean_abs_r": mean_abs_r,
            "max_abs_r": max_abs_r,
            "best_seed": seeds[best_idx],
            "best_r": r[np.arange(r.shape[0]), best_idx],
            "min_q": min_q,
            "n_seeds_signif": np.nansum(q < q_thresh, axis=1).astype(int),
            "n_seeds": np.sum(~np.isnan(r), axis=1).astype(int),
        }
    ).set_index("gene")
    return out
