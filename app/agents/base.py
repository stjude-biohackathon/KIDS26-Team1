"""Agent contract: a single AnalysisContext threaded through role-based agents.

Agents are pure-ish services with `run(ctx) -> ctx`. The Streamlit layer is the
orchestrator; heavy stages are cached upstream by cohort/seed/params keys.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol

import pandas as pd

from core.cohort import CohortSpec


@dataclass
class SeedSelection:
    source: str  # "INO80" | "SRCAP" | "Both" | "Custom"
    genes: list[str]  # resolved, de-duplicated

    def seed_hash(self) -> str:
        return f"{self.source}:{','.join(sorted(self.genes))}"


@dataclass
class Params:
    r_min: float = 0.3
    q_max: float = 0.05
    min_cohort_n: int = 15
    var_pctile: float = 25.0
    positive_only: bool = False
    method: str = "pearson"        # "pearson" | "spearman"
    population: str = "cohort"     # "cohort" | "pan_cancer" (candidate selection only)
    expr_tpm: float = 3.0          # "expressed" if TPM > this (matrix is log2(TPM+1))
    ess_thresh: float = -0.5       # "essential" if Chronos gene-effect <= this
    seed_on_frac: float = 0.5      # seed set "on" in a line if >= this frac of seeds expressed
    ce_penalty: float = 0.25       # multiplier on common-essential genes in ranking (0=exclude,1=none)
    pan_exclude_cohort: bool = False  # pan-cancer baseline excludes the selected cohort lines


@dataclass
class AnalysisContext:
    spec: CohortSpec
    seeds: SeedSelection
    params: Params
    # populated by agents
    models: list[str] = field(default_factory=list)
    cohort_n: int = 0
    corr_n: int = 0
    seeds_found_expr: list[str] = field(default_factory=list)
    seeds_found_effect: list[str] = field(default_factory=list)
    seeds_missing: list[str] = field(default_factory=list)
    coexpr: pd.DataFrame | None = None      # per-candidate aggregate
    coexpr_matrix: pd.DataFrame | None = None  # candidate x seed r
    essentiality: pd.DataFrame | None = None   # per-candidate % metrics (cohort lines)
    ranked: pd.DataFrame | None = None         # co-expressed candidates ranked by joint %
    ess_denom: int = 0                          # cohort lines with both expr+effect data
    coess: pd.DataFrame | None = None
    coess_matrix: pd.DataFrame | None = None
    consensus: pd.DataFrame | None = None
    warnings: list[str] = field(default_factory=list)
    narrative: str | None = None


class Agent(Protocol):
    name: str

    def run(self, ctx: AnalysisContext) -> AnalysisContext: ...
