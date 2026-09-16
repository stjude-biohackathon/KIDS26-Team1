"""Tests for harmonization. Run: uv run pytest app/tests -q"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

APP = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(APP))

from core.harmonize import (  # noqa: E402
    strip_entrez,
    clean_gene_columns,
    load_expression,
    load_gene_effect,
    load_model_meta,
    analyzable_models,
)

RAW = APP.parent / "user_data"

_EXPR_CSV = RAW / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
_EFFECT_CSV = RAW / "CRISPRGeneEffect.csv"
_MODEL_CSV = RAW / "Model.csv"
# Raw DepMap CSVs are not bundled in a distribution (large, licensed). Skip the
# CSV-reading tests when they are absent; the app runs from prebuilt parquets.
needs_raw = pytest.mark.skipif(
    not (_EXPR_CSV.exists() and _EFFECT_CSV.exists() and _MODEL_CSV.exists()),
    reason="raw DepMap CSVs not present (distribution ships prebuilt parquets)",
)


def test_strip_entrez():
    assert strip_entrez("TSPAN6 (7105)") == "TSPAN6"
    assert strip_entrez("A1BG (1)") == "A1BG"
    assert strip_entrez("PLAIN") == "PLAIN"
    assert clean_gene_columns(["BAD (572)", "YY1 (7528)"]) == ["BAD", "YY1"]


@needs_raw
def test_load_expression_collapses_to_one_row_per_model():
    expr = load_expression(RAW / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
    # 1,719 unique ModelIDs with a default profile
    assert expr.index.is_unique
    assert expr.shape[0] == 1719
    assert expr.index.name == "ModelID"
    assert expr.values.dtype == np.float32
    assert not any("(" in c for c in expr.columns)  # symbols cleaned


@needs_raw
def test_load_gene_effect():
    eff = load_gene_effect(RAW / "CRISPRGeneEffect.csv")
    assert eff.index.name == "ModelID"
    assert eff.shape[0] == 1208
    assert eff.values.dtype == np.float32
    assert not any("(" in c for c in eff.columns)


@needs_raw
def test_model_meta_matches_eda():
    meta = load_model_meta(RAW / "Model.csv")
    assert meta.shape[0] == 2154
    lineage = meta["OncotreeLineage"].value_counts()
    assert lineage["Lung"] == 261
    assert lineage["Lymphoid"] == 266
    age = meta["AgeCategory"].value_counts()
    assert age["Pediatric"] == 291
    assert age["Adult"] == 1480


@needs_raw
def test_analyzable_intersection():
    expr = load_expression(RAW / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
    eff = load_gene_effect(RAW / "CRISPRGeneEffect.csv")
    meta = load_model_meta(RAW / "Model.csv")
    common = analyzable_models(expr, eff, meta)
    assert len(common) > 900  # substantial overlap expected
    assert len(common) == len(set(common))
    assert set(common) <= set(meta.index)


@needs_raw
def test_seed_files_present_in_matrices():
    expr = load_expression(RAW / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
    eff = load_gene_effect(RAW / "CRISPRGeneEffect.csv")
    for fn, n in [("ino80.proteins.txt", 16), ("srcap.proteins.txt", 11)]:
        seeds = [g.strip() for g in (RAW / fn).read_text().splitlines() if g.strip()]
        assert len(seeds) == n
        # most seeds should be found in both matrices
        assert sum(g in expr.columns for g in seeds) >= n - 2
        assert sum(g in eff.columns for g in seeds) >= n - 2
