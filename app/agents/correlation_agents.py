"""Coexpression / Coessentiality agents: candidate x seed Pearson + FDR + aggregate.

Both share the same logic on different matrices, so a base class holds it.
Seed genes are removed from the candidate pool (a seed vs itself is r=1).
"""
from __future__ import annotations

import pandas as pd

from agents.base import AnalysisContext
from core.correlation import pearson_matrix, bh_fdr, aggregate_per_candidate


class _CorrelationAgent:
    name = "CorrelationAgent"
    modality = ""  # "expr" | "effect"

    def _matrix(self, ctx: AnalysisContext) -> pd.DataFrame:
        raise NotImplementedError

    def _seeds_found(self, ctx: AnalysisContext) -> list[str]:
        raise NotImplementedError

    def run(self, ctx: AnalysisContext) -> AnalysisContext:
        mat = self._matrix(ctx)
        seeds_found = self._seeds_found(ctx)
        if len(seeds_found) == 0 or mat.shape[0] < 3:
            empty = pd.DataFrame()
            self._store(ctx, empty, empty)
            ctx.warnings.append(f"[{self.name}] no seeds found or cohort too small.")
            return ctx

        seed_mat = mat[seeds_found]
        candidates = [g for g in mat.columns if g not in set(seeds_found)]
        cand_mat = mat[candidates]

        r_df, p_df = pearson_matrix(cand_mat, seed_mat, method=ctx.params.method)
        q = bh_fdr(p_df.to_numpy())
        q_df = pd.DataFrame(q, index=r_df.index, columns=r_df.columns)

        agg = aggregate_per_candidate(r_df, q_df, ctx.params.q_max)
        agg = agg.sort_values("mean_abs_r", ascending=False)
        self._store(ctx, agg, r_df, q_df)
        return ctx

    def _store(self, ctx, agg, r_df, q_df=None):
        raise NotImplementedError


class CoexpressionAgent(_CorrelationAgent):
    name = "CoexpressionAgent"
    modality = "expr"

    def _matrix(self, ctx):
        return ctx._cohort.expr  # type: ignore[attr-defined]

    def _seeds_found(self, ctx):
        return ctx.seeds_found_expr

    def _store(self, ctx, agg, r_df, q_df=None):
        ctx.coexpr = agg
        ctx.coexpr_matrix = r_df


class CoessentialityAgent(_CorrelationAgent):
    name = "CoessentialityAgent"
    modality = "effect"

    def _matrix(self, ctx):
        return ctx._cohort.effect  # type: ignore[attr-defined]

    def _seeds_found(self, ctx):
        return ctx.seeds_found_effect

    def _store(self, ctx, agg, r_df, q_df=None):
        ctx.coess = agg
        ctx.coess_matrix = r_df
