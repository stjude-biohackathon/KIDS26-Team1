"""DataAgent: load harmonized parquets, slice the cohort, variance-filter."""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd

from agents.base import AnalysisContext
from core.cohort import CohortData, select_models, variance_filter

PROCESSED = Path(__file__).resolve().parent.parent.parent / "data" / "processed"


@lru_cache(maxsize=1)
def _load_all() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expr = pd.read_parquet(PROCESSED / "expr_by_model.parquet")
    effect = pd.read_parquet(PROCESSED / "effect_by_model.parquet")
    meta = pd.read_parquet(PROCESSED / "model_meta.parquet")
    return expr, effect, meta


class DataAgent:
    name = "DataAgent"

    def build_cohort(self, ctx: AnalysisContext) -> CohortData:
        expr, effect, meta = _load_all()
        cohort_models = select_models(meta, ctx.spec)
        cohort_both = [m for m in cohort_models if m in expr.index and m in effect.index]

        # Correlation population: cohort (age+disease) or pan-cancer (all analyzable).
        if ctx.params.population == "pan_cancer":
            pop = [m for m in meta.index if m in expr.index and m in effect.index]
        else:
            pop = cohort_both

        expr_c = expr.loc[pop]
        effect_c = effect.loc[pop]
        # Variance-filter the CANDIDATE pool, but never drop seed genes we need.
        keep_expr = set(ctx.seeds.genes) & set(expr_c.columns)
        keep_effect = set(ctx.seeds.genes) & set(effect_c.columns)
        expr_c = variance_filter(expr_c, ctx.params.var_pctile, keep=keep_expr)
        effect_c = variance_filter(effect_c, ctx.params.var_pctile, keep=keep_effect)
        cd = CohortData(ctx.spec, pop, expr_c, effect_c)
        cd.cohort_models = cohort_both  # the age+disease cohort (for the card/warning)
        return cd

    def run(self, ctx: AnalysisContext) -> AnalysisContext:
        cohort = self.build_cohort(ctx)
        # cohort_n reflects the age+disease selection, not the correlation population
        ctx.models = cohort.cohort_models
        ctx.cohort_n = len(cohort.cohort_models)
        ctx.corr_n = cohort.n  # number of lines actually correlated
        ctx._cohort = cohort  # type: ignore[attr-defined]
        if ctx.cohort_n < ctx.params.min_cohort_n:
            ctx.warnings.append(
                f"Cohort has only {ctx.cohort_n} cell lines "
                f"(< min {ctx.params.min_cohort_n}); correlations are unstable."
            )
        if ctx.params.population == "pan_cancer":
            ctx.warnings.append(
                f"Pan-cancer mode: correlations use all {cohort.n} analyzable lines "
                "(disease selection affects only the drill-down, not the ranking)."
            )
        # seed presence
        seeds = ctx.seeds.genes
        ctx.seeds_found_expr = [g for g in seeds if g in cohort.expr.columns]
        ctx.seeds_found_effect = [g for g in seeds if g in cohort.effect.columns]
        found_any = set(ctx.seeds_found_expr) | set(ctx.seeds_found_effect)
        ctx.seeds_missing = [g for g in seeds if g not in found_any]
        return ctx
