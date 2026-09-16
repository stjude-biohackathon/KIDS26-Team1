"""Harmonization helpers for DepMap matrices.

Turns the raw DepMap CSVs (read-only, under user_data/) into tidy, ModelID-indexed
float32 tables with clean gene-symbol columns. No file under user_data/ is ever
written.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

# Columns of the expression CSV that are metadata, not genes.
EXPR_META_COLS = [
    "SequencingID",
    "ModelConditionID",
    "ModelID",
    "IsDefaultEntryForMC",
    "IsDefaultEntryForModel",
]

_SYMBOL_RE = re.compile(r"^(.*?)\s*\(\d+\)\s*$")


def strip_entrez(col: str) -> str:
    """'TSPAN6 (7105)' -> 'TSPAN6'. Leaves already-clean symbols untouched."""
    m = _SYMBOL_RE.match(col)
    return m.group(1) if m else col


def clean_gene_columns(cols) -> list[str]:
    return [strip_entrez(c) for c in cols]


def dedupe_columns(df: pd.DataFrame) -> pd.DataFrame:
    """If symbol stripping produced duplicate gene columns, keep the first."""
    if df.columns.duplicated().any():
        df = df.loc[:, ~df.columns.duplicated()]
    return df


def load_expression(path: str | Path) -> pd.DataFrame:
    """Load expression, collapse to one row per ModelID, clean gene symbols.

    The expression table is *profile-level* (1,775 rows). We keep the model's
    default profile (`IsDefaultEntryForModel == "Yes"`); if a ModelID has no
    default row we fall back to its MC-default, then to the first occurrence.
    Returns a float32 DataFrame indexed by ModelID (genes as columns).
    """
    df = pd.read_csv(path, low_memory=False)
    first_col = df.columns[0]
    if first_col not in EXPR_META_COLS and first_col.strip() == "":
        df = df.drop(columns=[first_col])

    def _truthy(series: pd.Series) -> pd.Series:
        return series.astype(str).str.strip().str.lower().isin({"yes", "true", "1"})

    is_model_default = _truthy(df["IsDefaultEntryForModel"])
    is_mc_default = _truthy(df["IsDefaultEntryForMC"])
    # Priority: model-default > mc-default > any. Stable sort keeps first.
    prio = np.where(is_model_default, 0, np.where(is_mc_default, 1, 2))
    order = np.lexsort((prio, df["ModelID"].values))  # sort by ModelID, then prio
    df = df.iloc[order]
    df = df.drop_duplicates(subset="ModelID", keep="first")

    gene_cols = [c for c in df.columns if c not in EXPR_META_COLS]
    out = df.set_index("ModelID")[gene_cols].copy()
    out.columns = clean_gene_columns(out.columns)
    out = dedupe_columns(out)
    out.index.name = "ModelID"
    return out.astype(np.float32)


def load_gene_effect(path: str | Path) -> pd.DataFrame:
    """Load CRISPR Chronos gene-effect, ModelID-indexed, clean gene symbols."""
    df = pd.read_csv(path, index_col=0, low_memory=False)
    df.index.name = "ModelID"
    df.columns = clean_gene_columns(df.columns)
    df = dedupe_columns(df)
    return df.astype(np.float32)


def load_model_meta(path: str | Path) -> pd.DataFrame:
    """Load Model.csv metadata, ModelID-indexed."""
    df = pd.read_csv(path, dtype=str, low_memory=False)
    df = df.set_index("ModelID")
    return df


def analyzable_models(
    expr: pd.DataFrame, effect: pd.DataFrame, meta: pd.DataFrame
) -> list[str]:
    """ModelIDs present in all three tables (order: meta order)."""
    common = set(expr.index) & set(effect.index) & set(meta.index)
    return [m for m in meta.index if m in common]
