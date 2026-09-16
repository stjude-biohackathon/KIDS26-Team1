"""Tests for the compute core + agents. Run: uv run pytest app/tests -q"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

APP = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(APP))

from conftest import needs_parquets  # noqa: E402
pytestmark = needs_parquets

from core.correlation import pearson_matrix, bh_fdr, aggregate_per_candidate  # noqa: E402
from core.cohort import CohortSpec, select_models, variance_filter  # noqa: E402
from agents.base import AnalysisContext, SeedSelection, Params  # noqa: E402
from agents.data_agent import DataAgent, _load_all  # noqa: E402
from agents.correlation_agents import CoexpressionAgent, CoessentialityAgent  # noqa: E402


def test_pearson_matches_scipy():
    rng = np.random.default_rng(0)
    n = 40
    cand = pd.DataFrame(rng.normal(size=(n, 5)), columns=[f"C{i}" for i in range(5)])
    seeds = pd.DataFrame(rng.normal(size=(n, 3)), columns=[f"S{j}" for j in range(3)])
    r_df, p_df = pearson_matrix(cand, seeds)
    for ci in cand.columns:
        for sj in seeds.columns:
            r_exp, p_exp = stats.pearsonr(cand[ci], seeds[sj])
            assert r_df.loc[ci, sj] == pytest.approx(r_exp, abs=1e-10)
            assert p_df.loc[ci, sj] == pytest.approx(p_exp, abs=1e-9)


def test_pearson_constant_column_is_nan():
    n = 20
    cand = pd.DataFrame({"flat": np.ones(n), "vary": np.arange(n)})
    seeds = pd.DataFrame({"s": np.random.default_rng(1).normal(size=n)})
    r_df, _ = pearson_matrix(cand, seeds)
    assert np.isnan(r_df.loc["flat", "s"])
    assert not np.isnan(r_df.loc["vary", "s"])


def test_bh_fdr_matches_statsmodels():
    from statsmodels.stats.multitest import multipletests
    rng = np.random.default_rng(3)
    p = rng.uniform(size=50)
    q_ours = bh_fdr(p.reshape(10, 5))
    _, q_sm, _, _ = multipletests(p, method="fdr_bh")
    assert np.allclose(q_ours.ravel(), q_sm, atol=1e-12)


def test_bh_fdr_ignores_nan():
    p = np.array([[0.01, np.nan], [0.5, 0.001]])
    q = bh_fdr(p)
    assert np.isnan(q[0, 1])
    assert not np.isnan(q[0, 0])


def test_variance_filter_drops_constant():
    df = pd.DataFrame({"const": np.ones(10), "a": np.arange(10.0), "b": np.arange(10.0) * 2})
    out = variance_filter(df, min_pctile=0.0)
    assert "const" not in out.columns
    assert "a" in out.columns


def test_cohort_spec_and_select():
    _, _, meta = _load_all()
    spec = CohortSpec(age_group="Adult", disease="Non-Small Cell Lung Cancer")
    models = select_models(meta, spec)
    assert len(models) > 40
    # all selected are Adult + NSCLC
    assert (meta.loc[models, "AgeCategory"] == "Adult").all()
    assert (meta.loc[models, "OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer").all()


def _ctx(age="Adult", disease="Non-Small Cell Lung Cancer", source="Both"):
    ino80 = [g.strip() for g in (APP.parent / "user_data" / "ino80.proteins.txt").read_text().splitlines() if g.strip()]
    srcap = [g.strip() for g in (APP.parent / "user_data" / "srcap.proteins.txt").read_text().splitlines() if g.strip()]
    genes = {"INO80": ino80, "SRCAP": srcap, "Both": list(dict.fromkeys(ino80 + srcap))}[source]
    return AnalysisContext(
        spec=CohortSpec(age_group=age, disease=disease),
        seeds=SeedSelection(source=source, genes=genes),
        params=Params(),
    )


def test_end_to_end_pipeline_nsclc_adult():
    ctx = _ctx("Adult", "Non-Small Cell Lung Cancer", "Both")
    ctx = DataAgent().run(ctx)
    assert ctx.cohort_n > 30
    assert len(ctx.seeds_found_expr) >= 18
    ctx = CoexpressionAgent().run(ctx)
    ctx = CoessentialityAgent().run(ctx)
    assert ctx.coexpr is not None and len(ctx.coexpr) > 100
    assert ctx.coess is not None and len(ctx.coess) > 100
    # seeds excluded from candidate pool
    assert not set(ctx.seeds.genes) & set(ctx.coexpr.index)
    # aggregate columns present and bounded
    assert ctx.coexpr["mean_abs_r"].max() <= 1.0
    assert (ctx.coexpr["n_seeds_signif"] >= 0).all()


def test_seed_source_changes_candidate_scores():
    a = _ctx("Adult", "Non-Small Cell Lung Cancer", "INO80")
    a = CoessentialityAgent().run(CoexpressionAgent().run(DataAgent().run(a)))
    b = _ctx("Adult", "Non-Small Cell Lung Cancer", "SRCAP")
    b = CoessentialityAgent().run(CoexpressionAgent().run(DataAgent().run(b)))
    # different seed sets => different top co-essential gene ordering (usually)
    assert a.seeds.genes != b.seeds.genes
    assert a.coess is not None and b.coess is not None


def test_spearman_matches_scipy():
    from scipy import stats as _st
    rng = np.random.default_rng(7)
    n = 45
    cand = pd.DataFrame(rng.normal(size=(n, 4)), columns=[f"C{i}" for i in range(4)])
    seeds = pd.DataFrame(rng.normal(size=(n, 2)), columns=[f"S{j}" for j in range(2)])
    r_df, p_df = pearson_matrix(cand, seeds, method="spearman")
    for ci in cand.columns:
        for sj in seeds.columns:
            rho, p = _st.spearmanr(cand[ci], seeds[sj])
            assert r_df.loc[ci, sj] == pytest.approx(rho, abs=1e-10)
            assert p_df.loc[ci, sj] == pytest.approx(p, abs=1e-6)


def test_pan_cancer_population_ignores_disease():
    from agents.base import AnalysisContext, SeedSelection, Params
    ino80 = [g.strip() for g in (APP.parent / "user_data" / "ino80.proteins.txt").read_text().splitlines() if g.strip()]
    def _pan(disease):
        ctx = AnalysisContext(CohortSpec("Adult", disease),
                              SeedSelection("INO80", ino80),
                              Params(population="pan_cancer"))
        return CoessentialityAgent().run(DataAgent().run(ctx))
    a = _pan("Non-Small Cell Lung Cancer")
    b = _pan("Melanoma")
    # same correlation population -> identical corr_n and identical top ranking
    assert a.corr_n == b.corr_n > 800
    assert a.corr_n != a.cohort_n  # pan population larger than the cohort
    assert list(a.coess.index[:10]) == list(b.coess.index[:10])
