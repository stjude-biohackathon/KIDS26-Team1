#!/usr/bin/env bash
# Setup script for the protein-complex-classifier skill.
# Installs the EXACT validated dependency versions via `uv sync` + `uv.lock`
# (shipped alongside pyproject.toml at the skill root) for a reproducible
# environment — not an open-ended `uv add` that could resolve to newer,
# untested package versions.
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$SKILL_ROOT"

echo "=================================================="
echo "  Protein-Complex Classifier — Environment Setup"
echo "=================================================="

if [ ! -f pyproject.toml ] || [ ! -f uv.lock ]; then
    echo "ERROR: pyproject.toml and/or uv.lock not found at $SKILL_ROOT."
    echo "Both files ship with this skill for reproducible installs —"
    echo "copy the entire skill directory (not just scripts/) before running setup."
    exit 1
fi

echo ""
echo "Installing pinned dependencies via 'uv sync' (uses uv.lock)..."
uv sync

echo ""
echo "Installing esm (ESM-C) with --no-deps..."
# torchtext/torchvision (declared esm deps) lack Python 3.13 wheels and are
# NOT needed for ESM-C embedding generation. Pinned to the validated version.
uv pip install "esm==3.4.0" --no-deps

echo ""
echo "Verifying installation..."
uv run python -c "
import numpy, pandas, sklearn, torch, transformers, joblib
from esm.models.esmc import ESMC
import esm
from transformers import T5EncoderModel, BertModel
from imblearn.over_sampling import SMOTE
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from bioservices import UniProt
print(f'esm={esm.__version__} torch={torch.__version__} transformers={transformers.__version__}')
print('✓ All dependencies verified (pinned versions)')
" 2>&1 | grep -v "^ESMC:" || true

echo ""
echo "=================================================="
echo "  Setup complete."
echo "=================================================="
echo ""
echo "Next: run the full pipeline"
echo "  uv run python scripts/run_pipeline.py \\"
echo "      --positive <positive_genes> --negative <negative_genes> \\"
echo "      --output-dir output/"
