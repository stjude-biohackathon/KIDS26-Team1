"""Tests for SynthesisAgent + essentiality metrics + viz. Run: uv run pytest app/tests -q"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

APP = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(APP))

from conftest import needs_parquets  # noqa: E402
pytestmark = needs_parquets

from agents.base import AnalysisContext, SeedSelection, Params  # noqa: E402
from core.cohort import CohortSpec  # noqa: E402
from agents.data_agent import DataAgent  # noqa: E402
from agents.correlation_agents import CoexpressionAgent  # noqa: E402
from agents.essentiality_agent import EssentialityAgent, compute_essentiality  # noqa: E402
from agents.synthesis_agent import SynthesisAgent, _common_essentials  # noqa: E402
from viz import plots  # noqa: E402


def _run(age="Adult", disease="Non-Small Cell Lung Cancer", source="Both", **pk):
    ino80 = [g.strip() for g in (APP.parent / "user_data" / "ino80.proteins.txt").read_text().splitlines() if g.strip()]
    srcap = [g.strip() for g in (APP.parent / "user_data" / "srcap.proteins.txt").read_text().splitlines() if g.strip()]
    genes = {"INO80": ino80, "SRCAP": srcap, "Both": list(dict.fromkeys(ino80 + srcap))}[source]
    ctx = AnalysisContext(CohortSpec(age_group=age, disease=disease),
                          SeedSelection(source, genes), Params(**pk))
    for a in (DataAgent(), CoexpressionAgent(), EssentialityAgent(), SynthesisAgent()):
        ctx = a.run(ctx)
    return ctx


def test_common_essentials_loaded():
    ce = _common_essentials()
    assert len(ce) > 500


def test_essentiality_metrics_bounds_and_denominator():
    ctx = _run(r_min=0.1, q_max=0.25)
    r = ctx.ranked
    assert r is not None and not r.empty
    assert ctx.ess_denom > 30
    # percentages in [0,100]; joint <= min(coexpr, essential)
    for col in ["pct_coexpressed", "pct_essential", "pct_joint"]:
        assert r[col].between(0, 100).all()
    assert (r["pct_joint"] <= r["pct_coexpressed"] + 1e-9).all()
    assert (r["pct_joint"] <= r["pct_essential"] + 1e-9).all()
    # counts consistent with denominator
    assert (r["n_lines"] == ctx.ess_denom).all()
    assert (r["n_joint"] <= r["n_coexpressed"]).all()


def test_ranked_sorted_by_selectivity():
    ctx = _run(r_min=0.1, q_max=0.25)
    assert ctx.ranked["selectivity_score"].is_monotonic_decreasing


def test_conditional_metric_definition():
    ctx = _run(r_min=0.1, q_max=0.25)
    r = ctx.ranked
    # pct_essential_given_coexpr == 100 * n_joint / n_coexpressed where n_coexpressed>0
    mask = r["n_coexpressed"] > 0
    expected = 100.0 * r.loc[mask, "n_joint"] / r.loc[mask, "n_coexpressed"]
    assert (abs(r.loc[mask, "pct_essential_given_coexpr"] - expected) < 1e-6).all()


def test_common_essential_flag_present():
    ctx = _run(r_min=0.1, q_max=0.25)
    assert "common_essential" in ctx.ranked.columns
    assert ctx.ranked["common_essential"].dtype == bool


def test_selectivity_score_and_ranking():
    ctx = _run(r_min=0.1, q_max=0.25, ce_penalty=0.25)
    r = ctx.ranked
    for col in ["pct_essential_pan", "selectivity_delta", "selectivity_score"]:
        assert col in r.columns
    # ranked by selectivity_score descending
    assert r["selectivity_score"].is_monotonic_decreasing
    # delta == cohort essential% - pan essential%
    assert (abs(r["selectivity_delta"] - (r["pct_essential"] - r["pct_essential_pan"])) < 1e-6).all()
    # score formula holds
    p = ctx.params
    ce_mult = r["common_essential"].map(lambda x: p.ce_penalty if x else 1.0)
    expected = r["pct_joint"] * (1 + r["selectivity_delta"] / 100.0).clip(lower=0.0) * ce_mult
    assert (abs(r["selectivity_score"] - expected) < 1e-6).all()


def test_ce_penalty_zero_pushes_common_essential_down():
    ctx = _run(r_min=0.1, q_max=0.25, ce_penalty=0.0)
    r = ctx.ranked
    # with penalty 0, common-essential genes get score 0 and sink to the bottom
    top = r.head(50)
    assert (~top["common_essential"]).sum() >= (top["common_essential"]).sum()


def test_empty_when_strict():
    ctx = _run(r_min=0.99, q_max=0.001)
    assert ctx.ranked is not None and ctx.ranked.empty


def test_plots_return_figures():
    ctx = _run(r_min=0.1, q_max=0.25)
    assert isinstance(plots.essentiality_scatter(ctx.ranked, "t"), go.Figure)
    assert isinstance(plots.ranked_bar(ctx.ranked, 20, "t"), go.Figure)
    assert isinstance(plots.essentiality_scatter(pd.DataFrame(), "t"), go.Figure)  # empty-safe


def test_gene_narrative_deterministic_grounded():
    from agents.narrative import gene_narrative_deterministic, gene_facts
    ctx = _run(r_min=0.1, q_max=0.25)
    top3 = list(ctx.ranked.index[:3])
    txt = gene_narrative_deterministic(ctx, top3 + ["NOT_A_GENE"])
    for g in top3:
        assert g in txt
        f = gene_facts(ctx, g)
        assert f["rank"] >= 1 and f["of_total"] == len(ctx.ranked)
    assert "not among" in txt  # missing gene reported


def test_gene_narrative_deterministic_highlights_selectivity_delta():
    """The headline selectivity_delta clause must be wrapped in this app's
    bounded ==highlight== syntax (not raw HTML/markdown color) -- rendered
    safely by ui.markdown_with_highlights() on-screen and by pdf_export.py's
    _split_inline() in the PDF."""
    from agents.narrative import gene_narrative_deterministic, gene_facts
    ctx = _run(r_min=0.1, q_max=0.25)
    top = ctx.ranked.index[0]
    txt = gene_narrative_deterministic(ctx, [top])
    f = gene_facts(ctx, top)
    assert f"==delta {f['selectivity_delta']:+.0f} pts," in txt
    assert txt.count("==") >= 2  # opening + closing marker present


def test_plot_highlight_genes_additive_and_backward_compatible():
    """volcano()/essentiality_scatter()/ranked_bar() gained an optional
    highlight_genes param for the Insights tab's cross-tab plot embedding --
    must be purely additive (omitting it reproduces the prior trace count)."""
    from viz import plots
    import plotly.graph_objects as go
    ctx = _run(r_min=0.1, q_max=0.25)
    top = list(ctx.ranked.index[:2])

    v_plain = plots.volcano(ctx.coexpr, 0.1, 0.25, "t")
    v_hi = plots.volcano(ctx.coexpr, 0.1, 0.25, "t", highlight_genes=top)
    assert isinstance(v_hi, go.Figure)
    assert len(v_hi.data) == len(v_plain.data) + 1  # one overlay trace added

    s_plain = plots.essentiality_scatter(ctx.ranked, "t")
    s_hi = plots.essentiality_scatter(ctx.ranked, "t", highlight_genes=top)
    assert len(s_hi.data) == len(s_plain.data) + 1

    b_plain = plots.ranked_bar(ctx.ranked, 10, "t")
    b_hi = plots.ranked_bar(ctx.ranked, 10, "t", highlight_genes=[top[0]])
    assert len(b_hi.data) == len(b_plain.data)  # bar recolors in place, no new trace
    assert list(b_hi.data[0].marker.color) != list(b_plain.data[0].marker.color)


def test_plot_highlight_genes_absent_gene_is_a_noop():
    from viz import plots
    ctx = _run(r_min=0.1, q_max=0.25)
    v_plain = plots.volcano(ctx.coexpr, 0.1, 0.25, "t")
    v_hi = plots.volcano(ctx.coexpr, 0.1, 0.25, "t", highlight_genes=["NOT_A_REAL_GENE"])
    assert len(v_hi.data) == len(v_plain.data)  # no overlay trace for a gene not present


def test_heatmap_order_param_selective_only():
    from viz import plots
    ctx = _run(r_min=0.1, q_max=0.25)
    selective = list(ctx.ranked[~ctx.ranked["common_essential"]].index)
    fig = plots.heatmap(ctx.coexpr_matrix, 20, ctx.coexpr, "t", order=selective)
    import plotly.graph_objects as go
    assert isinstance(fig, go.Figure)
    # y-axis genes must be a subset of the selective ordering (no common-essential)
    ce = set(ctx.ranked[ctx.ranked["common_essential"]].index)
    ydata = list(fig.data[0].y)
    assert not (set(ydata) & ce)


def test_pan_exclude_cohort_toggle():
    inc = _run(r_min=0.1, q_max=0.25, pan_exclude_cohort=False)
    exc = _run(r_min=0.1, q_max=0.25, pan_exclude_cohort=True)
    a = inc.ranked[["pct_essential_pan", "selectivity_delta"]]
    b = exc.ranked.reindex(a.index)[["pct_essential_pan", "selectivity_delta"]]
    # cohort essential% is unchanged; only the pan baseline (and thus delta) shifts
    assert (inc.ranked["pct_essential"].reindex(a.index)
            .equals(exc.ranked["pct_essential"].reindex(a.index)))
    # for at least some genes the pan baseline differs between the two modes
    assert not a["pct_essential_pan"].equals(b["pct_essential_pan"])
    # delta = cohort - pan still holds under exclusion
    r = exc.ranked
    assert (abs(r["selectivity_delta"] - (r["pct_essential"] - r["pct_essential_pan"])) < 1e-6).all()
