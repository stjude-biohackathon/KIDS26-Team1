<h1 align="center"><b><u>Team 1 | Bio Hackathon 2026</u></b></h1>
  <!-- 2. The Title (Centered, Bolded, Underlined) -->
  <h1 align="center"><b>Agentic AI Workflows to Discover Pediatric-Specific INO80/SRCAP Chromatin Remodeler Subunits and Adaptors for Cancer Dependencies</b></h1>
  <!-- 3. The Tagline -->

$${\Large {\color{#DC143C}\underline{\textbf{\textsf{AI}}}}\ \textsf{assisted}\ {\color{#DC143C}\underline{\textbf{\textsf{S}}}}\textsf{elective}\ {\color{#DC143C}\underline{\textbf{\textsf{C}}}}\textsf{ancer-dependency}\ {\color{#DC143C}\underline{\textbf{\textsf{R}}}}\textsf{anking}\ \textsf{of}\ {\color{#DC143C}\underline{\textbf{\textsf{A}}}}\textsf{ssociated}\ {\color{#DC143C}\underline{\textbf{\textsf{P}}}}\textsf{roteins}}$$

<p align="center">
  <img src="assets/SCRAP-AI_logo_recolored.svg" alt="Project Logo" width="400">
  <br>

</p>

Multi-agent Streamlit app that ranks genes **co-expressed with the INO80/SRCAP
chromatin-remodeling complex** (or any other proteins of user interest) by their **cohort-selective essentiality**, folds
in **STRING physical/functional interaction** evidence, and builds grounded
**target dossiers** (literature + protein/drug annotations) for candidate genes.

## Team members

- Pandurang Kolekar (Co-Lead)
- Vinesh Vinayachandran (Co-Lead)
- Satish Sati
- Abhinandita Dash (Computational Biology, St. Jude Children's Research Hospital)
- Sheetal Bhatara
- Yuyu Zhang
- Sreerag M
- Vishakan Deivanayagam
- Dhanus Chinnakalipatti Thilagaa Ramesh
- Nithish Marimuthu

## About

Each human cell packs roughly two meters of DNA into a nucleus only a few micrometers wide, wrapping the genome around histone proteins into nucleosomes — the repeating beads-on-a-string units of chromatin. That packaging is far from static: ATP-dependent chromatin remodelers like INO80 and SRCAP continuously reposition nucleosomes and swap in the histone variant H2A.Z at gene promoters, deciding which genes a cell can even reach. These remodelers are universal — every cell depends on them, which is exactly why their core catalytic machinery is essentially undruggable. But both complexes are modular, built from swappable accessory subunits, and a specific subset of these appears to be one that pediatric tumors are quietly addicted to while healthy cells barely need them — a pediatric-specific weak spot hidden inside a universal machine. This matters because most childhood cancers carry few DNA mutations at all; their vulnerabilities are written into chromatin state, not sequence. Pinpointing these cancer-specific accessory submodules opens the door to disabling a tumor's remodeling machinery while leaving the essential core untouched in healthy tissue — precisely the kind of once-undruggable target that modern degraders like PROTACs and molecular glues can now reach, pointing toward safer, more selective treatments for children with catastrophic cancers.

## 1. Prerequisites

- **[uv](https://docs.astral.sh/uv/)** (installs Python + all dependencies in an isolated venv):
  ```bash
  # macOS / Linux
  curl -LsSf https://astral.sh/uv/install.sh | sh
  # Windows (PowerShell)
  powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
  ```
- ~1 GB free disk (dependencies + bundled data). No admin rights needed.

## 2. One command

```bash
# macOS / Linux
./run.sh
# Windows
run.cmd
```

That syncs the pinned dependencies, verifies the setup, and launches the app at
**http://localhost:8501**.

### Or step by step
```bash
uv sync                                   # create venv + install exact pinned deps
uv run python check_setup.py              # verify install + data
uv run streamlit run app/streamlit_app.py # launch
```

## 3. Data (REQUIRED for the lite package)

This **lite package ships code only** — no prebuilt data. You must download the
DepMap files once and build the parquets before first launch.

Download these 4 files from the
[DepMap portal](https://depmap.org/portal/data_page/?tab=allData) (26Q1 Public)
into `user_data/`:

| File | Content |
|---|---|
| `CRISPRGeneEffect.csv` | Chronos gene-effect |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | log2(TPM+1) expression |
| `Model.csv` | cell-line metadata |
| `CRISPRInferredCommonEssentials.csv` | common-essential list  |

Then build the harmonized parquets (~10s):

```bash
uv run python app/build_harmonized.py
```

`check_setup.py` tells you exactly what's missing and what to run. Cite DepMap
when using these data.

> Prefer a zero-setup try? Use the **ready-to-run** bundle (`scrap-ai.tar.gz`),
> which ships the prebuilt parquets — no download needed.

## 4. Verify / smoke-test

```bash
uv run python check_setup.py     # deps + data readiness
uv run pytest app/tests -q       # 136 tests (network-gated STRING/RAG tests skip offline)
```

## 5. Optional: LLM narrative & external evidence

The core app is fully **offline and keyless** — all ranking, analysis, and the
**deterministic** gene narrative work with no key. Two optional features reach the
network:

- **Insights → literature/drug/STRING dossiers** use keyless public APIs
  (Europe PMC, UniProt, ChEMBL, STRING) — no key required.
- **Insights → BioMCP** use keyless BioMCP APIs
    - BioMCP (https://biomcp.org/) can be installed using `uv tool install biomcp-cli`
    - If VPN is enabled, BioMCP may not work.
- **Insights → "Generate with LLM"** uses a cloud model if you provide a key,
  else local **Ollama** (`llama3.1:8b`) as per availability.
    - Refer to Ollama download/installation instructions here: https://docs.ollama.com/quickstart
    - Ollama server should be running at http://localhost:11434

### Using your own LLM key

Bring your own key (BYOK) any **one** of three ways — nothing is bundled, and your
key is only sent to the provider you choose:

1. **`.env` file** (best for a shared/persistent setup):
   ```bash
   cp app/.env.example .env      # in the package root
   # edit .env and uncomment ONE line, e.g.:
   #   ANTHROPIC_API_KEY=sk-ant-...
   ```
   The app auto-loads `.env` from the package root (or `app/`) at startup.
   `.env` is git-ignored so it won't leak if you re-share the folder.

2. **In-app field** (nothing to edit): launch → **Insights** tab →
   **"LLM backend / API key"** → pick provider, paste key, **Use key**
   (kept in memory for the session only), then **Test LLM connection**.

3. **Shell environment variable**:
   ```bash
   export ANTHROPIC_API_KEY="sk-ant-..."   # macOS/Linux
   set ANTHROPIC_API_KEY=sk-ant-...        # Windows cmd
   ```

**Provider priority** (first key present wins), then local Ollama:
```
OPENROUTER_API_KEY → OPENAI_API_KEY → ANTHROPIC_API_KEY → Ollama (llama3.1:8b)
```
Override the model with `OPENROUTER_MODEL` / `OPENAI_MODEL` / `ANTHROPIC_MODEL`
(the in-app dropdown also lists valid models for a cloud key). No key and no
Ollama? The deterministic "Generate narrative" button still works.

> Provider notes: a bad model id shows the real error + a live model list; models
> that reject `temperature` (e.g. Claude Sonnet 5) are retried automatically.

## 6. What's included

```
scrap-ai/
├── run.sh / run.cmd          # one-command launcher
├── check_setup.py            # install + data verifier
├── pyproject.toml            # pinned dependencies
├── .python-version           # Python 3.14
├── .streamlit/config.toml    # brand theme
├── app/                      # application code + tests + README
├── assets/                   # SCRAP-AI logo + icon
├── data/processed/           # (lite: empty; built by app/build_harmonized.py)
├── user_data/                # seed gene lists (+ where DepMap CSVs go)
├── GLOSSARY.md               # Glossary of terms and scores used in the app
└── ARCHITECTURE.md           # stack + diagrams
```

## Troubleshooting

- **`uv: command not found`** → install uv (step 1), reopen the terminal.
- **`check_setup.py` says data missing** → the parquets should be bundled; if you
  removed them, download the CSVs (step 3) and run `app/build_harmonized.py`.
- **Port in use** → `uv run streamlit run app/streamlit_app.py --server.port 8600`.
- **Windows** requires Git Bash only if you use `run.sh`; `run.cmd` works in plain cmd.

See `app/README.md` for feature details and `ARCHITECTURE.md` for the design.
