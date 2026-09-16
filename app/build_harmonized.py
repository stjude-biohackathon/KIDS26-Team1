"""One-time build: raw DepMap CSVs (read-only) -> float32 parquets under data/processed/.

Run:  uv run python app/build_harmonized.py
Never writes under user_data/.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import pandas as pd

from core.harmonize import (
    load_expression,
    load_gene_effect,
    load_model_meta,
    analyzable_models,
)

SANDBOX = Path(__file__).resolve().parent.parent
RAW = SANDBOX / "user_data"
OUT = SANDBOX / "data" / "processed"

EXPR_CSV = RAW / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
EFFECT_CSV = RAW / "CRISPRGeneEffect.csv"
MODEL_CSV = RAW / "Model.csv"


def _digest(path: Path, nbytes: int = 50 * 1024 * 1024) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        h.update(f.read(nbytes))
    return h.hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print("Loading expression (collapse to ModelID)...")
    expr = load_expression(EXPR_CSV)
    print(f"  expr: {expr.shape[0]} models x {expr.shape[1]} genes")

    print("Loading CRISPR gene-effect...")
    effect = load_gene_effect(EFFECT_CSV)
    print(f"  effect: {effect.shape[0]} models x {effect.shape[1]} genes")

    print("Loading Model.csv metadata...")
    meta = load_model_meta(MODEL_CSV)
    print(f"  meta: {meta.shape[0]} models x {meta.shape[1]} fields")

    common = analyzable_models(expr, effect, meta)
    print(f"  analyzable (in all three): {len(common)} models")

    expr.to_parquet(OUT / "expr_by_model.parquet")
    effect.to_parquet(OUT / "effect_by_model.parquet")
    meta.to_parquet(OUT / "model_meta.parquet")

    manifest = {
        "built_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "source": "DepMap 26Q1 Public (user_data/, read-only)",
        "expr": {"models": int(expr.shape[0]), "genes": int(expr.shape[1]),
                 "src_md5_50mb": _digest(EXPR_CSV)},
        "effect": {"models": int(effect.shape[0]), "genes": int(effect.shape[1]),
                   "src_md5_50mb": _digest(EFFECT_CSV)},
        "meta": {"models": int(meta.shape[0]), "fields": int(meta.shape[1]),
                 "src_md5_50mb": _digest(MODEL_CSV)},
        "analyzable_models": len(common),
        "genes_in_both_matrices": int(len(set(expr.columns) & set(effect.columns))),
        "build_seconds": round(time.time() - t0, 1),
    }
    (OUT / "data_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2))
    print(f"Done in {manifest['build_seconds']}s -> {OUT}")


if __name__ == "__main__":
    main()
