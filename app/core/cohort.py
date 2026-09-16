"""Cohort selection from Model.csv metadata + candidate-gene variance filtering."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

PEDIATRIC_CATS = ("Pediatric",)
ADULT_CATS = ("Adult",)


@dataclass
class CohortSpec:
    age_group: str  # "Pediatric" | "Adult"
    disease: str  # OncotreePrimaryDisease value, or "All"
    include_fetus: bool = False

    def age_cats(self) -> tuple[str, ...]:
        if self.age_group == "Pediatric":
            return PEDIATRIC_CATS + (("Fetus",) if self.include_fetus else ())
        return ADULT_CATS

    def cohort_id(self) -> str:
        return "|".join(
            [
                self.age_group,
                "fetus" if self.include_fetus else "nofetus",
                self.disease,
            ]
        )


def select_models(meta: pd.DataFrame, spec: CohortSpec) -> list[str]:
    """ModelIDs matching the age + primary-disease filter."""
    m = meta
    sel = m["AgeCategory"].isin(spec.age_cats())
    if spec.disease and spec.disease != "All":
        sel &= m["OncotreePrimaryDisease"] == spec.disease
    return list(m.index[sel])


def disease_options(meta: pd.DataFrame, spec_age_cats: tuple[str, ...]) -> pd.Series:
    """Disease -> cohort size for the current age group (for the UI selectbox)."""
    sub = meta[meta["AgeCategory"].isin(spec_age_cats)]
    return sub["OncotreePrimaryDisease"].value_counts()


def variance_filter(
    mat: pd.DataFrame, min_pctile: float = 25.0, min_var: float = 1e-8,
    keep: set[str] | None = None,
) -> pd.DataFrame:
    """Drop near-constant / low-variance gene columns within the cohort.

    Keeps genes whose variance is > the `min_pctile` percentile of the nonzero
    variances (and strictly > min_var). Columns are genes. Any gene in `keep`
    (e.g. seed genes) is retained regardless of variance, as long as it is not
    perfectly constant (which would make correlation undefined).
    """
    v = mat.var(axis=0, ddof=1)
    keep = set(keep or ())
    nonzero = v[v > min_var]
    if nonzero.empty:
        forced = [g for g in mat.columns if g in keep and v.get(g, 0) > min_var]
        return mat[forced]
    thresh = float(np.percentile(nonzero.to_numpy(), min_pctile))
    passes = (v > min_var) & (v >= thresh)
    if keep:
        forced = mat.columns.isin(keep) & (v.to_numpy() > min_var)
        passes = passes.to_numpy() | forced
    else:
        passes = passes.to_numpy()
    # boolean mask over columns, original order preserved
    return mat.loc[:, passes]


@dataclass
class CohortData:
    spec: CohortSpec
    models: list[str]
    expr: pd.DataFrame  # models x genes (cohort, variance-filtered)
    effect: pd.DataFrame  # models x genes (cohort, variance-filtered)
    n: int = field(init=False)
    cohort_models: list[str] = field(default_factory=list)

    def __post_init__(self):
        self.n = len(self.models)
        if not self.cohort_models:
            self.cohort_models = self.models
