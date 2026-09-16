"""Tests for the STRING-db module. Network-gated (skip offline).
Run: uv run pytest app/tests -q
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

APP = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(APP))

from agents import string_db  # noqa: E402

SEEDS = ("INO80", "RUVBL1", "RUVBL2", "ACTL6A", "ACTB", "YEATS4", "SRCAP",
         "VPS72", "DMAP1", "YY1")


def _online() -> bool:
    try:
        string_db._post("tsv-no-header/get_string_ids",
                        {"identifiers": "TP53", "species": 9606, "caller_identity": "scrap-ai-test"},
                        timeout=12)
        return True
    except Exception:
        return False


def test_subunit_scores_high():
    if not _online():
        pytest.skip("offline")
    r = string_db.string_scores("ACTL6A", SEEDS, 0.4)
    # ACTL6A is a real INO80/SRCAP subunit -> strong functional AND physical
    assert r["functional"]["max"] >= 0.9
    assert r["physical"]["max"] >= 0.9
    assert r["functional"]["n_seeds_ge_cutoff"] >= 5
    assert r["functional"]["best_seed"] in SEEDS


def test_noncomplex_gene_low_physical():
    if not _online():
        pytest.skip("offline")
    # LDB1 is co-expressed but not a physical subunit
    r = string_db.string_scores("LDB1", SEEDS, 0.4)
    assert r["physical"]["max"] <= r["functional"]["max"] + 1e-9
    assert 0.0 <= r["functional"]["max"] <= 1.0


def test_scores_bounded_and_shaped():
    if not _online():
        pytest.skip("offline")
    r = string_db.string_scores("YEATS4", SEEDS, 0.4)
    for k in ("functional", "physical"):
        a = r[k]
        assert 0.0 <= a["mean"] <= a["max"] <= 1.0
        assert isinstance(a["linked_seeds"], dict)
    assert r["string_url"].startswith("https://string-db.org")


def test_network_image_png_or_none():
    if not _online():
        pytest.skip("offline")
    img = string_db.network_image("ACTL6A", SEEDS[:5])
    assert img is None or img[:4] == b"\x89PNG"


def test_string_rerank_and_matrices():
    if not _online():
        pytest.skip("offline")
    import pandas as pd
    from agents.string_rank import string_rerank
    # small ranked frame: a subunit (high STRING) + a weak gene
    df = pd.DataFrame({"selectivity_score": [100.0, 80.0]}, index=["ACTL6A", "LDB1"])
    rr, phys_m, func_m, edges = string_rerank(df, SEEDS, top_n=2,
                                              w_sel=0.5, w_phys=0.35, w_func=0.15, cutoff=0.4)
    assert "string_score" in rr.columns and "sel_norm" in rr.columns
    # normalized weighted-sum stays in [0,1]
    assert (rr["string_score"] >= 0).all() and (rr["string_score"] <= 1.0 + 1e-9).all()
    # min-max: within top-N, sel_norm spans 0..1
    assert abs(rr["sel_norm"].max() - 1.0) < 1e-9 and abs(rr["sel_norm"].min()) < 1e-9
    # ACTL6A (subunit) should get a bigger physical component than LDB1
    assert rr.loc["ACTL6A", "string_phys"] >= rr.loc["LDB1", "string_phys"]
    assert phys_m.shape[1] == len(SEEDS) and func_m.shape[1] == len(SEEDS)
    assert rr["string_score"].is_monotonic_decreasing
    # weights recorded and normalized to ~1
    w = rr.attrs["weights"]
    assert abs(w["selectivity"] + w["physical"] + w["functional"] - 1.0) < 1e-6


def test_string_plots_shapes():
    import pandas as pd
    import plotly.graph_objects as go
    from viz import plots
    mat = pd.DataFrame([[0.9, 0.1], [0.0, 0.5]], index=["A", "B"], columns=["S1", "S2"])
    assert isinstance(plots.string_heatmap(mat, "t"), go.Figure)
    edges = [{"source": "A", "target": "S1", "physical": 0.9, "functional": 0.8}]
    assert isinstance(plots.string_network(edges, {"A": 100.0}, "physical", "t"), go.Figure)
    assert isinstance(plots.string_network([], {}, "physical", "t"), go.Figure)  # empty-safe


def test_cutoff_thresholds_score():
    if not _online():
        pytest.skip("offline")
    import pandas as pd
    from agents.string_rank import string_rerank
    df = pd.DataFrame({"selectivity_score": [100.0, 100.0]}, index=["LDB1", "SMARCA4"])
    # LDB1 functional best ~0.389: counts at cutoff 0.3, excluded at 0.5
    lo, _, _, _ = string_rerank(df, SEEDS, top_n=2, w_sel=0.0, w_phys=0.0, w_func=1.0, cutoff=0.3)
    hi, _, _, _ = string_rerank(df, SEEDS, top_n=2, w_sel=0.0, w_phys=0.0, w_func=1.0, cutoff=0.5)
    assert lo.loc["LDB1", "string_score"] > hi.loc["LDB1", "string_score"]
    assert hi.loc["LDB1", "string_score"] == 0.0  # weak link zeroed at higher cutoff
    # raw display value is preserved regardless of cutoff
    assert lo.loc["LDB1", "string_func"] == hi.loc["LDB1", "string_func"]
