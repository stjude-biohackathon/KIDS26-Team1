"""Verify SCRAP-AI is installed and data is ready.

Checks:
  1. dependencies import
  2. seed files present
  3. either prebuilt parquets exist, OR raw DepMap CSVs exist (then offers rebuild)

Run:  uv run python check_setup.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "user_data"
PROC = ROOT / "data" / "processed"

EXPECTED_CSVS = {
    "CRISPRGeneEffect.csv": "Chronos gene-effect (1208 lines × ~18.5k genes)",
    "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv": "log2(TPM+1) expression",
    "Model.csv": "cell-line metadata",
    "CRISPRInferredCommonEssentials.csv": "common-essential list",
}
SEED_FILES = ["ino80.proteins.txt", "srcap.proteins.txt"]
PARQUETS = ["expr_by_model.parquet", "effect_by_model.parquet", "model_meta.parquet"]

ok = True


def check(cond, good, bad):
    global ok
    print(("  ✓ " if cond else "  ✗ ") + (good if cond else bad))
    if not cond:
        ok = False


print("SCRAP-AI setup check\n" + "-" * 40)

print("Dependencies:")
for mod in ("numpy", "pandas", "scipy", "statsmodels", "pyarrow", "streamlit", "plotly"):
    try:
        __import__(mod)
        check(True, f"{mod} importable", "")
    except Exception as e:  # noqa: BLE001
        check(False, "", f"{mod} missing ({e}) — run `uv sync`")

print("Seed files:")
for f in SEED_FILES:
    check((RAW / f).is_file(), f"user_data/{f}", f"user_data/{f} MISSING")

print("Data (need EITHER prebuilt parquets OR raw CSVs):")
have_parquets = all((PROC / p).is_file() for p in PARQUETS)
have_csvs = all((RAW / c).is_file() for c in EXPECTED_CSVS)

if have_parquets:
    check(True, "prebuilt parquets present (data/processed/) — ready to run", "")
elif have_csvs:
    print("  • raw DepMap CSVs present but parquets not built.")
    print("    Build them with:  uv run python app/build_harmonized.py")
    ok = False
else:
    print("  ✗ No prebuilt parquets AND no raw CSVs.")
    print("    Download these 4 files from https://depmap.org/portal/data_page/?tab=allData")
    print("    (DepMap 26Q1 Public) into user_data/ :")
    for c, desc in EXPECTED_CSVS.items():
        print(f"       - {c}   ({desc})")
    print("    then:  uv run python app/build_harmonized.py")
    ok = False

# quick shape sanity if parquets exist
if have_parquets:
    try:
        import pandas as pd
        e = pd.read_parquet(PROC / "effect_by_model.parquet")
        x = pd.read_parquet(PROC / "expr_by_model.parquet")
        m = pd.read_parquet(PROC / "model_meta.parquet")
        both = len(set(e.index) & set(x.index) & set(m.index))
        print(f"  • effect {e.shape}, expr {x.shape}, meta {m.shape}; "
              f"{both} analyzable models")
    except Exception as ex:  # noqa: BLE001
        check(False, "", f"parquets unreadable ({ex})")

print("-" * 40)
if ok:
    print("READY ✓   launch with:  uv run streamlit run app/streamlit_app.py")
    sys.exit(0)
else:
    print("NOT READY — resolve the ✗ items above.")
    sys.exit(1)
