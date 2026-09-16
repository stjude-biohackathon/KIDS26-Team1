#!/usr/bin/env bash
# SCRAP-AI launcher (macOS/Linux). Requires uv: https://docs.astral.sh/uv/
set -e
cd "$(dirname "$0")"
command -v uv >/dev/null 2>&1 || { echo "uv not found. Install: curl -LsSf https://astral.sh/uv/install.sh | sh"; exit 1; }
echo "Syncing dependencies (first run downloads them)…"
uv sync
echo "Checking setup…"
uv run python check_setup.py || { echo "Setup incomplete — see messages above."; exit 1; }
echo "Launching SCRAP-AI at http://localhost:8501 …"
uv run streamlit run app/streamlit_app.py
