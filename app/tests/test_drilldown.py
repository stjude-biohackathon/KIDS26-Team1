"""Tests for drill-down vectors. Run: uv run pytest app/tests -q"""
from __future__ import annotations

import sys
from pathlib import Path

APP = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(APP))

from conftest import needs_parquets  # noqa: E402
pytestmark = needs_parquets

from core.cohort import CohortSpec  # noqa: E402
from core.drilldown import cohort_gene_vectors  # noqa: E402


def test_drilldown_vectors_expr():
    spec = CohortSpec(age_group="Adult", disease="Non-Small Cell Lung Cancer")
    vec = cohort_gene_vectors(spec, "INO80B", "RUVBL1", "expr")
    assert vec is not None and not vec.empty
    assert list(vec.columns) == ["x", "y", "cell_line", "disease"]
    assert len(vec) > 30
    assert vec["cell_line"].notna().all()


def test_drilldown_vectors_effect():
    spec = CohortSpec(age_group="Adult", disease="Non-Small Cell Lung Cancer")
    vec = cohort_gene_vectors(spec, "INO80", "ACTL6A", "effect")
    assert vec is not None and not vec.empty


def test_drilldown_missing_gene_returns_none():
    spec = CohortSpec(age_group="Adult", disease="Non-Small Cell Lung Cancer")
    assert cohort_gene_vectors(spec, "NOT_A_GENE", "RUVBL1", "expr") is None
