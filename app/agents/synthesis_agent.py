"""SynthesisAgent: annotate the ranked co-expressed/essential candidates with
common-essential flags and produce the interpretive narrative.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd

from agents.base import AnalysisContext

RAW = Path(__file__).resolve().parent.parent.parent / "user_data"


@lru_cache(maxsize=1)
def _common_essentials() -> set[str]:
    fp = RAW / "CRISPRInferredCommonEssentials.csv"
    if not fp.exists():
        return set()
    s = pd.read_csv(fp)
    col = s.columns[0]
    genes = s[col].astype(str).str.replace(r"\s*\(\d+\)$", "", regex=True)
    return set(genes)


class SynthesisAgent:
    name = "SynthesisAgent"

    def run(self, ctx: AnalysisContext) -> AnalysisContext:
        ce = _common_essentials()
        ranked = ctx.ranked
        if ranked is not None and not ranked.empty:
            ranked = ranked.copy()
            ranked["common_essential"] = ranked.index.isin(ce)
            # Selectivity-weighted score:
            #   score = pct_joint * (1 + selectivity_delta/100) * ce_mult
            # delta = cohort essential% - pan-cancer essential%  (in [-100,100])
            #   (1 + delta/100) in [0,2]: boosts genes more essential HERE, suppresses
            #   anti-selective ones; housekeeping (delta~0) stays ~1x then hit by penalty.
            #   ce_mult = ce_penalty for DepMap common-essential genes, else 1.
            p = ctx.params
            delta_factor = (1.0 + ranked["selectivity_delta"] / 100.0).clip(lower=0.0)
            ce_mult = ranked["common_essential"].map(lambda x: p.ce_penalty if x else 1.0)
            ranked["selectivity_score"] = ranked["pct_joint"] * delta_factor * ce_mult
            ranked = ranked.sort_values("selectivity_score", ascending=False)
            ctx.ranked = ranked
        from agents.narrative import deterministic_narrative
        ctx.narrative = deterministic_narrative(ctx)
        return ctx
