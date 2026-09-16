"""Tests for the drug-target scoring pipeline.

Offline: Scorer / CancerResolver logic on synthetic data, stage-to-phase mapping,
plot shapes. Network-gated: live smoke test against a known result (ALK +
Neuroblastoma -> CERITINIB Score 5, matching the collaborator's own sample CSV).
Run: uv run pytest app/tests -q
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

APP = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(APP))

from agents import drug_targets as dt  # noqa: E402


def _online() -> bool:
    try:
        import requests
        requests.get("https://mygene.info/v3/query", params={"q": "symbol:TP53", "size": 1},
                     timeout=10).raise_for_status()
        return True
    except Exception:
        return False


# ---------------- offline: stage-to-phase mapping ----------------
def test_stage_to_phase_mapping():
    ot = dt.OpenTargetsClient
    assert ot._stage_to_phase("APPROVAL") == 4.0
    assert ot._stage_to_phase("PHASE_3") == 3.0
    assert ot._stage_to_phase("PHASE_1_2") == 1.5
    assert ot._stage_to_phase("EARLY_PHASE_1") == 0.5
    assert ot._stage_to_phase(None) is None
    assert ot._stage_to_phase("bogus_stage") is None
    assert ot._stage_to_phase("phase_2") == 2.0  # case-insensitive


# ---------------- offline: 3-tier Scorer ----------------
def test_score_3tier_fda_plus_cancer():
    entry = {"max_phase": 4, "fda_approved": True, "cancer_indication": True}
    assert dt.Scorer.score(entry, gene_cancers=[]) == 3


def test_score_3tier_fda_only():
    entry = {"max_phase": 4, "fda_approved": True, "cancer_indication": False}
    assert dt.Scorer.score(entry, gene_cancers=[]) == 2


def test_score_3tier_experimental():
    entry = {"max_phase": 0, "fda_approved": False, "cancer_indication": False}
    assert dt.Scorer.score(entry, gene_cancers=[]) == 1


# ---------------- offline: 5-tier Scorer (using real CancerResolver, no network) ----------------
def test_score_5tier_approved_for_query_cancer():
    resolver = dt.CancerResolver()
    cancer_info = {"disease_name": "Neuroblastoma", "cancer_terms": {"neuroblastoma"}}
    entry = {
        "max_phase": 4, "fda_approved": True,
        "indications": [{"term": "neuroblastoma", "phase": "4", "mesh": ""}],
    }
    score, label, indicated, phase, other = dt.Scorer.score_cancer_specific(
        entry, resolver, cancer_info)
    assert score == 5
    assert indicated is True


def test_score_5tier_approved_other_cancer_only():
    resolver = dt.CancerResolver()
    cancer_info = {"disease_name": "Neuroblastoma", "cancer_terms": {"neuroblastoma"}}
    entry = {
        "max_phase": 4, "fda_approved": True,
        "indications": [{"term": "breast carcinoma", "phase": "4", "mesh": ""}],
    }
    score, label, indicated, phase, other = dt.Scorer.score_cancer_specific(
        entry, resolver, cancer_info)
    assert score == 3
    assert indicated is False
    assert "breast carcinoma" in other


def test_score_5tier_dgidb_only_no_indications():
    resolver = dt.CancerResolver()
    cancer_info = {"disease_name": "Neuroblastoma", "cancer_terms": {"neuroblastoma"}}
    entry = {"max_phase": -1, "fda_approved": True, "indications": []}
    score, label, indicated, phase, other = dt.Scorer.score_cancer_specific(
        entry, resolver, cancer_info)
    assert score == 2  # approved but cancer-specificity unknown -> capped at 2
    assert indicated is False


# ---------------- offline: plot shapes ----------------
def test_drug_target_plots_shapes():
    from viz import plots
    import plotly.graph_objects as go
    df = pd.DataFrame({
        "Gene": ["ALK", "ALK", "MYCN"],
        "Drug_Name": ["CERITINIB", "CRIZOTINIB", "X"],
        "Score": [5, 3, 1],
    })
    assert isinstance(plots.drug_target_heatmap(df, 15, 5, "t"), go.Figure)
    assert isinstance(plots.drug_target_score_bar(df, "t"), go.Figure)
    assert isinstance(plots.drug_target_heatmap(pd.DataFrame(), 15, 5, "t"), go.Figure)
    assert isinstance(plots.drug_target_score_bar(pd.DataFrame(), "t"), go.Figure)


# ---------------- network-gated: the actual bug-fix validation ----------------
def test_known_drugs_query_no_flaky_endpoint():
    """The Open Targets replacement for the old ChEMBL name-lookup: one call,
    keyed by Ensembl ID, must return CERITINIB with a Neuroblastoma indication."""
    if not _online():
        pytest.skip("offline")
    ot = dt.OpenTargetsClient()
    known = ot.query_known_drugs("ENSG00000171094")  # ALK
    assert "ceritinib" in known
    terms = [i["term"].lower() for i in known["ceritinib"]["indications"]]
    assert any("neuroblastoma" in t for t in terms)


def test_full_pipeline_reproduces_collaborator_validated_result():
    """Ground truth from user_data/drug_targets/drug_target_results.csv:
    ALK + Neuroblastoma -> CERITINIB, Score 5, FDA-approved, indicated."""
    if not _online():
        pytest.skip("offline")
    df = dt.run_pipeline(["ALK"], output_dir=str(APP.parent / "data" / "drug_targets_cache" / "_test"),
                         cancer_type="Neuroblastoma", min_score=2, skip_plots=True)
    assert df is not None and not df.empty
    row = df[df["Drug_Name"] == "CERITINIB"]
    assert not row.empty
    assert int(row.iloc[0]["Score"]) == 5
    assert row.iloc[0]["FDA_Approved"] == "Yes"
