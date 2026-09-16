@echo off
REM SCRAP-AI launcher (Windows). Requires uv: https://docs.astral.sh/uv/
cd /d "%~dp0"
where uv >nul 2>nul || (echo uv not found. Install from https://docs.astral.sh/uv/ && exit /b 1)
echo Syncing dependencies (first run downloads them)...
uv sync || exit /b 1
echo Checking setup...
uv run python check_setup.py || (echo Setup incomplete - see messages above. && exit /b 1)
echo Launching SCRAP-AI at http://localhost:8501 ...
uv run streamlit run app/streamlit_app.py
