"""Shared test fixtures / skip markers."""
from __future__ import annotations

from pathlib import Path

import pytest

_APP = Path(__file__).resolve().parent.parent
_PROC = _APP.parent / "data" / "processed"

# Many tests run the full pipeline, which reads the harmonized parquets. In a
# code-only distribution those are absent until `app/build_harmonized.py` runs,
# so skip (not fail) when they are missing.
needs_parquets = pytest.mark.skipif(
    not (_PROC / "expr_by_model.parquet").exists()
    or not (_PROC / "effect_by_model.parquet").exists()
    or not (_PROC / "model_meta.parquet").exists(),
    reason="harmonized parquets not built (run app/build_harmonized.py)",
)
