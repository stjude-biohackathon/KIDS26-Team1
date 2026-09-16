"""Precomputed protein-complex-classifier predictions (see .pi/skills/protein-complex-classifier
and the session's investigation of user_data/protein_classifier/).

This does NOT run the classifier -- it only loads and looks up predictions that
were precomputed elsewhere and dropped into user_data/ as a CSV. Read-only;
never writes to user_data/.

Important context carried over from testing the live classifier pipeline this
session: its ProteomeLM-S feature block is context-dependent on which other
proteins were embedded in the same batch, so absolute probabilities/labels
from one batch are NOT comparable to another batch. This loader deliberately
uses ONLY ONE file (the caller-specified "extended list") per lookup, rather
than merging predictions from files that may have been generated as different
batches -- mixing batches into one table would silently combine incomparable
numbers.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

REQUIRED_COLS = ("Gene", "UniProt_ID", "Probability", "Prediction", "Confidence")


def load_predictions(path: Path) -> pd.DataFrame:
    """Load a precomputed predictions CSV, indexed by upper-cased gene symbol.

    Returns an empty DataFrame (not an error) if the file is missing or
    malformed, so the caller can degrade to "not generated/cached" for every
    gene rather than crash the tab.
    """
    if not path.exists():
        return pd.DataFrame(columns=REQUIRED_COLS).set_index("Gene")
    try:
        df = pd.read_csv(path)
    except Exception:
        return pd.DataFrame(columns=REQUIRED_COLS).set_index("Gene")
    if not set(REQUIRED_COLS).issubset(df.columns):
        return pd.DataFrame(columns=REQUIRED_COLS).set_index("Gene")
    df = df.copy()
    df["Gene"] = df["Gene"].astype(str).str.upper().str.strip()
    df = df.drop_duplicates(subset="Gene", keep="first").set_index("Gene")
    return df


def lookup(genes: list[str], preds: pd.DataFrame) -> pd.DataFrame:
    """One row per requested gene: prediction fields if cached, else flagged
    not-cached (never dropped, never fabricated)."""
    rows = []
    for g in genes:
        gu = g.upper().strip()
        if gu in preds.index:
            r = preds.loc[gu]
            rows.append({
                "gene": gu, "cached": True,
                "UniProt_ID": r.get("UniProt_ID", ""),
                "Probability": r.get("Probability", float("nan")),
                "Prediction": r.get("Prediction", ""),
                "Confidence": r.get("Confidence", ""),
            })
        else:
            rows.append({
                "gene": gu, "cached": False,
                "UniProt_ID": "", "Probability": float("nan"),
                "Prediction": "Not generated/cached", "Confidence": "",
            })
    return pd.DataFrame(rows).set_index("gene")
