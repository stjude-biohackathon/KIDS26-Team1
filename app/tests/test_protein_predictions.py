"""Tests for app/agents/protein_predictions.py (precomputed classifier CSV loader)."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from agents import protein_predictions as pp  # noqa: E402


def _write_csv(tmp_path, rows):
    p = tmp_path / "preds.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    return p


def test_load_predictions_missing_file_returns_empty(tmp_path):
    df = pp.load_predictions(tmp_path / "does_not_exist.csv")
    assert df.empty
    assert df.index.name == "Gene"


def test_load_predictions_malformed_file_returns_empty(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("not,the,right,columns\n1,2,3,4\n")
    df = pp.load_predictions(p)
    assert df.empty


def test_load_predictions_uppercases_and_dedupes(tmp_path):
    p = _write_csv(tmp_path, [
        {"Gene": "brd7", "UniProt_ID": "Q9NPI1", "Probability": 0.97,
         "Prediction": "INTERACTOR", "Confidence": "High"},
        {"Gene": "BRD7", "UniProt_ID": "Q9NPI1", "Probability": 0.50,
         "Prediction": "NON-INTERACTOR", "Confidence": "Low"},  # duplicate, should be dropped
    ])
    df = pp.load_predictions(p)
    assert list(df.index) == ["BRD7"]
    assert df.loc["BRD7", "Probability"] == 0.97  # first occurrence kept


def test_lookup_cached_and_not_cached(tmp_path):
    p = _write_csv(tmp_path, [
        {"Gene": "BRD7", "UniProt_ID": "Q9NPI1", "Probability": 0.97,
         "Prediction": "INTERACTOR", "Confidence": "High"},
    ])
    preds = pp.load_predictions(p)
    out = pp.lookup(["brd7", "NOTAGENE"], preds)
    assert out.loc["BRD7", "cached"] == True  # noqa: E712
    assert out.loc["BRD7", "Prediction"] == "INTERACTOR"
    assert out.loc["NOTAGENE", "cached"] == False  # noqa: E712
    assert out.loc["NOTAGENE", "Prediction"] == "Not generated/cached"
    assert pd.isna(out.loc["NOTAGENE", "Probability"])


def test_lookup_never_fabricates_a_probability_for_missing_gene(tmp_path):
    """The core requirement: missing genes must never get an invented score."""
    preds = pp.load_predictions(tmp_path / "missing.csv")  # empty
    out = pp.lookup(["ANYGENE"], preds)
    assert out.loc["ANYGENE", "cached"] == False  # noqa: E712
    assert pd.isna(out.loc["ANYGENE", "Probability"])


def test_lookup_preserves_requested_order():
    preds = pd.DataFrame(
        {"UniProt_ID": ["A1", "B1"], "Probability": [0.9, 0.1],
         "Prediction": ["INTERACTOR", "NON-INTERACTOR"], "Confidence": ["High", "Low"]},
        index=pd.Index(["GENEB", "GENEA"], name="Gene"))
    out = pp.lookup(["GENEA", "GENEB"], preds)
    assert list(out.index) == ["GENEA", "GENEB"]


def test_real_extended_predictions_file_loads_if_present():
    """If the user's actual uploaded file is present, sanity-check its shape
    (skips cleanly if the file isn't there -- e.g. a fresh checkout)."""
    real = (Path(__file__).resolve().parent.parent.parent / "user_data"
            / "top_neuroblastoma_predictions" / "neuroblastoma_extended_list_predictions.csv")
    if not real.exists():
        pytest.skip("neuroblastoma_extended_list_predictions.csv not present")
    df = pp.load_predictions(real)
    assert not df.empty
    assert set(df["Prediction"].unique()) <= {"INTERACTOR", "NON-INTERACTOR"}
    assert df["Probability"].between(0, 1).all()
