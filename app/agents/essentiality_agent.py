"""Per-cell-line essentiality metrics for co-expressed candidate genes.

For the SELECTED cancer cohort (age + disease), and for each candidate gene that
is co-expressed with the INO80/SRCAP seed proteins, compute how often across the
cohort's cell lines the gene is:
  - co-expressed  (candidate expressed AND the seed set is "on" in that line)
  - essential     (Chronos gene-effect <= ess_thresh)
  - both (joint)  -> the primary ranking metric

"Expressed" uses the log2(TPM+1) matrix: log2(TPM+1) > log2(expr_tpm + 1).
The seed set is "on" in a line if >= seed_on_frac of its genes are expressed.
Metrics are always cohort-scoped, even when candidate selection used pan-cancer.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from agents.base import AnalysisContext, Params
from agents.data_agent import _load_all
from core.cohort import select_models


def _expr_cut(expr_tpm: float) -> float:
    return float(np.log2(expr_tpm + 1.0))


def compute_essentiality(
    ctx: AnalysisContext, candidates: list[str]
) -> tuple[pd.DataFrame, int]:
    """Return (metrics_df indexed by gene, denom = cohort lines with both data)."""
    expr, effect, meta = _load_all()
    p: Params = ctx.params

    # cohort lines (age + disease) present in BOTH matrices — always cohort-scoped
    cohort = select_models(meta, ctx.spec)
    lines = [m for m in cohort if m in expr.index and m in effect.index]
    denom = len(lines)
    if denom == 0:
        return pd.DataFrame(), 0

    cut = _expr_cut(p.expr_tpm)
    if not candidates:
        return pd.DataFrame(), denom

    # seed set "on" per line: fraction of found seed genes expressed
    seeds_expr = [g for g in ctx.seeds.genes if g in expr.columns]
    if seeds_expr:
        seed_block = expr.loc[lines, seeds_expr].to_numpy() > cut   # lines x k
        seed_on = seed_block.mean(axis=1) >= p.seed_on_frac         # lines,
    else:
        seed_on = np.zeros(denom, dtype=bool)

    # candidate genes present in both matrices
    cand = [g for g in candidates if g in expr.columns and g in effect.columns]
    if not cand:
        return pd.DataFrame(), denom

    cand_expr = expr.loc[lines, cand].to_numpy() > cut               # lines x m
    cand_ess = effect.loc[lines, cand].to_numpy() <= p.ess_thresh    # lines x m

    coexpr_line = cand_expr & seed_on[:, None]                       # lines x m
    joint = coexpr_line & cand_ess

    n_coexpr = coexpr_line.sum(axis=0)
    n_ess = cand_ess.sum(axis=0)
    n_joint = joint.sum(axis=0)

    # Pan-cancer essentiality baseline: fraction of analyzable lines (both
    # matrices) where the gene is essential -> selectivity delta. Optionally
    # exclude the selected cohort lines for a purer in- vs out-of-cohort contrast.
    pan_lines = [m for m in effect.index if m in expr.index]
    if p.pan_exclude_cohort:
        cohort_set = set(lines)
        pan_lines = [m for m in pan_lines if m not in cohort_set]
    if not pan_lines:  # degenerate: cohort == whole panel
        pan_lines = [m for m in effect.index if m in expr.index]
    pan_ess = effect.loc[pan_lines, cand].to_numpy() <= p.ess_thresh
    pct_ess_pan = 100.0 * pan_ess.sum(axis=0) / len(pan_lines)
    pct_ess_cohort = 100.0 * n_ess / denom

    with np.errstate(divide="ignore", invalid="ignore"):
        pct_cond = np.where(n_coexpr > 0, 100.0 * n_joint / n_coexpr, np.nan)

    df = pd.DataFrame(
        {
            "gene": cand,
            "pct_coexpressed": 100.0 * n_coexpr / denom,
            "pct_essential": pct_ess_cohort,
            "pct_essential_pan": pct_ess_pan,
            "selectivity_delta": pct_ess_cohort - pct_ess_pan,
            "pct_joint": 100.0 * n_joint / denom,
            "pct_essential_given_coexpr": pct_cond,
            "n_coexpressed": n_coexpr.astype(int),
            "n_essential": n_ess.astype(int),
            "n_joint": n_joint.astype(int),
            "n_lines": denom,
        }
    ).set_index("gene")
    return df, denom


class EssentialityAgent:
    """Compute per-cell-line essentiality metrics for co-expressed candidates,
    then rank the co-expressed hits by % of lines co-expressed AND essential."""

    name = "EssentialityAgent"

    def run(self, ctx: AnalysisContext) -> AnalysisContext:
        cox = ctx.coexpr
        if cox is None or cox.empty:
            ctx.essentiality = pd.DataFrame()
            ctx.ranked = pd.DataFrame()
            return ctx

        # co-expression hits: pass the correlation thresholds
        p = ctx.params
        score = cox["mean_r"] if p.positive_only else cox["mean_abs_r"]
        hits = cox[(score.abs() >= p.r_min) & (cox["min_q"] <= p.q_max)
                   & (cox["n_seeds_signif"] >= 1)]
        candidates = list(hits.index)

        metrics, denom = compute_essentiality(ctx, candidates)
        ctx.essentiality = metrics
        ctx.ess_denom = denom

        if metrics.empty:
            ctx.ranked = pd.DataFrame()
            return ctx

        # join co-expression summary with per-line essentiality metrics.
        # Final selectivity scoring + sort happens in SynthesisAgent (it owns the
        # common_essential flag). Provisional sort by joint % here.
        ranked = hits.join(metrics, how="inner")
        ranked = ranked.sort_values("pct_joint", ascending=False)
        ctx.ranked = ranked
        return ctx
