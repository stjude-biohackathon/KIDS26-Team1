"""GTEx v11 normal-tissue median-TPM loader (generalized: any tissue, not just
adrenal/brain). Ported from the collaborator's `gtex_loader.py`, adapted to be
disease-agnostic and consistent with this app's caching conventions.

Cache location: data/gtex/gtex_v11_median_tpm.gct.gz (bundle-friendly, like
data/processed/ for DepMap). If the collaborator's own copy already exists
under user_data/gtex_other_expression/data/, we reuse it (read-only copy) to
avoid a redundant ~10 MB download -- user_data/ itself is never written.
"""
from __future__ import annotations

import gzip
import shutil
import urllib.request
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

GTEX_MEDIAN_TPM_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/"
    "GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz"
)

_SANDBOX = Path(__file__).resolve().parent.parent.parent
_CACHE_DIR = _SANDBOX / "data" / "gtex"
_CACHE_FILE = _CACHE_DIR / "gtex_v11_median_tpm.gct.gz"
_COLLAB_COPY = _SANDBOX / "user_data" / "gtex_other_expression" / "data" / "gtex_v11_median_tpm.gct.gz"


def _ensure_local_file() -> Path:
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    if _CACHE_FILE.exists() and _CACHE_FILE.stat().st_size > 1_000_000:
        return _CACHE_FILE
    if _COLLAB_COPY.exists() and _COLLAB_COPY.stat().st_size > 1_000_000:
        shutil.copyfile(_COLLAB_COPY, _CACHE_FILE)  # read-only source, write only our cache copy
        return _CACHE_FILE
    urllib.request.urlretrieve(GTEX_MEDIAN_TPM_URL, _CACHE_FILE)
    return _CACHE_FILE


@lru_cache(maxsize=1)
def _load_full_gct() -> pd.DataFrame:
    """Full GTEx median-TPM-by-tissue table, indexed by gene symbol (Description)."""
    path = _ensure_local_file()
    with gzip.open(path, "rt") as f:
        lines = f.readlines()
    header = lines[2].rstrip("\n").split("\t")
    rows = [ln.rstrip("\n").split("\t") for ln in lines[3:]]
    df = pd.DataFrame(rows, columns=header)
    df = df.drop(columns=["Name"]).rename(columns={"Description": "Gene"}).set_index("Gene")
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


@lru_cache(maxsize=1)
def gtex_tissues() -> list[str]:
    """All GTEx tissue column names available (for the UI multiselect)."""
    return list(_load_full_gct().columns)


def load_gtex(gene_symbols: tuple[str, ...], tissues: tuple[str, ...]) -> pd.DataFrame:
    """Median TPM for the requested genes across the requested GTEx tissues.

    Rows = gene symbol (only genes found in GTEx), columns = requested tissues.
    """
    full = _load_full_gct()
    cols = [t for t in tissues if t in full.columns]
    sub = full[cols].copy()
    wanted = {g.upper() for g in gene_symbols}
    sub = sub[sub.index.str.upper().isin(wanted)]
    sub.index = sub.index.str.upper()
    sub = sub[~sub.index.duplicated(keep="first")]
    return sub


def normal_gate(gtex: pd.DataFrame, cutoff: float, mode: str = "median") -> pd.Series:
    """Normal-tissue silence gate across the requested GTEx tissues.

    mode='median' (default): median TPM across the selected tissues < cutoff.
      Matches the original collaborator pipeline's per-organ-aggregate approach
      (e.g. Adrenal + brain-region median), tolerant of one atypical subregion.
    mode='all': EVERY selected tissue individually < cutoff (stricter).
    """
    if gtex.empty:
        return pd.Series(dtype=bool)
    if mode == "all":
        return (gtex < cutoff).all(axis=1)
    return gtex.median(axis=1) < cutoff
