# SCRAP-AI — Selective Cancer-dependency Ranking of Associated Proteins

Multi-agent Streamlit app that ranks genes **co-expressed with the INO80/SRCAP
chromatin-remodeling complex** by their **cohort-selective essentiality**, folds
in **STRING physical/functional interaction** evidence, and builds grounded
**target dossiers** (literature + protein/drug annotations) for candidate genes.
DepMap 26Q1 Public. `user_data/` is treated as **read-only**.

## Setup & run

```bash
# 1. One-time: build harmonized parquets from the read-only CSVs (~10s)
uv run python app/build_harmonized.py

# 2. Launch the app
uv run streamlit run app/streamlit_app.py

# tests
uv run pytest app/tests -q            # 136 tests (network-gated ones skip offline)
```

## How it works

Four role-based agents threaded through a single `AnalysisContext`, orchestrated
by Streamlit (heavy stages cached with `@st.cache_data`; the orchestrator uses
**no LLM** — all science is deterministic pandas/numpy):

| Agent | Role |
|---|---|
| **DataAgent** | load parquets, slice cohort (age+disease), variance-filter candidates (seeds exempt) |
| **CoexpressionAgent** | candidate × seed Pearson/Spearman on log2(TPM+1) expression + BH-FDR — selects co-expressed candidates |
| **EssentialityAgent** | per-cohort-cell-line metrics on co-expressed candidates: % co-expressed, % essential, **% joint**, % essential\|co-expr; + pan-cancer baseline and selectivity score |
| **SynthesisAgent** | common-essential flag, **selectivity-weighted score**, deterministic narrative |

### Ranking

A candidate is **co-expressed** because its expression correlates with the
INO80/SRCAP seed proteins across cell lines. Then, within the selected cancer
cohort, per cell line L:
- *co-expressed in L* = candidate expressed (TPM > threshold) **and** ≥ X% of seed genes expressed;
- *essential in L* = Chronos gene-effect ≤ threshold.

**Selectivity-weighted score** (primary rank):

```
selectivity_score = pct_joint × (1 + selectivity_Δ/100) × common_essential_penalty
```

where `selectivity_Δ` = cohort essential% − pan-cancer essential% (rewards
cohort-specific dependencies; pan-essential housekeeping genes are down-weighted).
The pan-cancer baseline is over all analyzable lines; a sidebar toggle can
**exclude the selected cohort** from that baseline for a purer in- vs
out-of-cohort contrast.

**Optional STRING-weighted re-rank** (top-N, normalized weighted-sum, all ∈ [0,1]):

```
string_score = w_sel·sel_norm + w_phys·physical_max + w_func·functional_max
```

`sel_norm` = min-max of `selectivity_score` across the top-N; the three weights
auto-normalize to sum 1 (defaults 0.5 / 0.35 / 0.15). `physical_max` /
`functional_max` are the candidate's strongest STRING scores to any seed.

## UI

- **Sidebar:** theme toggle (light/dark); 🧬 Cohort (age, Fetus toggle, disease with
  live sizes); 🌱 Seeds (INO80 / SRCAP / Both / Custom — paste or upload a `.txt` gene
  list, editable); ⚙️ Method & thresholds
  (Pearson/Spearman, cohort/pan-cancer, min |r|, FDR q, TPM, Chronos, seed-on
  fraction, common-essential penalty); 🔗 STRING interactions (opt-in fold-into-ranking
  with selectivity/physical/functional weights + cutoff).
- **Tabs:**
  - **Overview** — deterministic synthesis of the run
  - **Co-expression** — candidate selection: volcano + click-to-drill per-cell-line scatter
  - **Essentiality & ranking** — % co-expressed vs % essential scatter, ranked bar,
    table, hide-common-essential toggle
  - **Heatmap** — selective candidate × seed correlation, score-ordered
  - **Interactions** — STRING candidate × seed heatmap + interaction network graph
    (physical/functional toggle; node size = selectivity, edge width = STRING) + re-ranked table
  - **Tissue Specificity** — GTEx normal-tissue silence + TCGA pan-cancer + a matched
    tumor cohort (TCGA/TARGET), for gene(s) you type. Tissue and tumor-cohort picks are
    best-effort auto-suggested from the current cancer type (always overridable), with
    configurable thresholds and a Plotly heatmap. On-demand only (Xena network queries);
    does not affect ranking. Generalized from a collaborator's Neuroblastoma-specific
    pipeline (see `user_data/gtex_other_expression/`); zero new dependencies (the Xena
    TOIL-hub query is hand-reimplemented over `urllib`, not the `xenaPython` package).
    A fetal/CELLxGENE-Census confirmation layer exists in that standalone pipeline but is
    **not integrated** here (heavy install, ~1–2 GB one-time download) — run it separately
    if needed.
  - **Drug Targets** — batch scoring for gene(s) you type (default: top-ranked candidates):
    queries MyGene.info, ChEMBL, DGIdb, and Open Targets for known drugs, then scores each
    gene-drug pair by clinical maturity + cancer relevance (3-tier cancer-agnostic, or
    5-tier cancer-specific pre-filled from the current disease). Optional PRISM/DepMap
    sensitivity enrichment (off by default, ~100 MB download). Interactive Plotly heatmap
    + score-tier bar chart, CSV export. Ported from a collaborator's pipeline (see
    `user_data/drug_targets/README.md` for full scoring rationale); the original's flaky
    per-drug ChEMBL name-lookup was replaced with one Open Targets call per gene (4–5×
    faster, fully reliable — see Notes).
  - **Insights** — type gene(s) → deterministic narrative (always) + optional LLM-generated
    **structured report** (6 fixed sections: Executive Summary, Cohort-Selective Dependency
    Evidence, Protein Complex/STRING Association, Drug-Design & Therapeutic Relevance,
    Literature Evidence, Experimental Validation & Testing Guidelines), grounded in
    literature / protein / drug / STRING evidence, plus a **deterministic Comparative
    Summary Table** (code-computed, never LLM-authored) appended after the LLM's text —
    tissue-specificity is deliberately excluded from the model's own prompt (its gate/score
    nuance is too easy to flatten into an overconfident sentence) but still shown via that
    table. Each narrative has its own **\U0001f4c4 Export as PDF** button.
  - **Model Predictions** — displays precomputed predictions from an external protein-language-model
    classifier (`.pi/skills/protein-complex-classifier`), loaded read-only from a CSV export under
    `user_data/top_neuroblastoma_predictions/`. This tab never runs the classifier live. A gene not
    present in the loaded export shows **Not generated/cached** rather than a guessed value.
    ⚠️ that classifier has a feature block sensitive to which other proteins were embedded in the
    same batch — probabilities are only comparable within one export, not across a different file
    or an ad-hoc small-batch re-run (empirically inflates scores; see Notes).
- All plots interactive (Plotly); every table is CSV-downloadable.

## Layout

```
app/
├── streamlit_app.py        # orchestrator + UI
├── build_harmonized.py     # CSV -> parquet (one-time)
├── ui.py                   # theme (light/dark), header, cards, chips, Plotly templates
├── agents/
│   ├── base.py             # AnalysisContext, SeedSelection, Params, Agent protocol
│   ├── data_agent.py       # load parquets, cohort slice, variance filter
│   ├── correlation_agents.py  # CoexpressionAgent (+ legacy CoessentialityAgent)
│   ├── essentiality_agent.py  # per-line % metrics + pan-cancer baseline
│   ├── synthesis_agent.py  # common-essential flag + selectivity score + narrative
│   ├── narrative.py        # deterministic + LLM narrative (6-section report), provider
│   │                        # backend, RAG, comparative-table builder
│   ├── knowledge.py        # Europe PMC / UniProt / ChEMBL dossiers
│   ├── string_db.py        # STRING functional/physical scores + network image
│   ├── string_rank.py      # STRING-weighted top-N re-rank + matrices/edges
│   ├── xena_client.py      # minimal UCSC Xena client (plain urllib, no xenaPython)
│   ├── gtex_agent.py       # GTEx v11 normal-tissue median-TPM loader
│   ├── tissue_specificity.py  # GTEx + TCGA/TARGET triage: gates + selectivity score
│   └── drug_targets.py     # ChEMBL/DGIdb/Open Targets/PRISM batch drug-target scoring
├── pdf_export.py           # markdown-subset -> PDF renderer for both narrative types
├── core/                   # harmonize, correlation, cohort, drilldown
├── viz/plots.py            # Plotly figures (volcano, heatmap, scatter, network)
└── tests/                  # 136 tests (uv run pytest app/tests -q)
assets/                     # SCRAP-AI logo + icon (svg/png)
.streamlit/config.toml      # brand theme
data/processed/             # harmonized parquets (gitignored)
```

## LLM backend (optional — Insights narrative only)

The gene-relevance narrative uses a cloud model when a key is present (priority
`OPENROUTER_API_KEY` → `OPENAI_API_KEY` → `ANTHROPIC_API_KEY`), else local
**Ollama** (`llama3.1:8b`). Override the cloud model with `OPENROUTER_MODEL` /
`OPENAI_MODEL` / `ANTHROPIC_MODEL`. The deterministic "Generate narrative" needs
**no model at all**.

- Provide a key three ways: a **`.env`** file in the sandbox root or `app/`
  (auto-loaded; copy `app/.env.example`), a **shell export**, or the **in-UI key
  field** on the Insights tab → "LLM backend / API key". Use **Test LLM
  connection** to confirm. Existing shell env vars win over `.env`.
- Robust to provider quirks: bad model ids surface the real error + a live model
  list; models that reject `temperature` (e.g. Claude Sonnet 5) are retried
  without it.
- Each cloud provider's endpoint can be redirected with a `{PROVIDER}_BASE_URL`
  env var (e.g. `ANTHROPIC_BASE_URL`) — this is how to use **Claude via Microsoft
  Azure AI Foundry** instead of api.anthropic.com directly: set
  `ANTHROPIC_API_KEY` to your Foundry API key, `ANTHROPIC_BASE_URL` to
  `https://{your-resource}.services.ai.azure.com/anthropic/v1`, and
  `ANTHROPIC_MODEL` to your Foundry **deployment name** (not necessarily the
  public Claude model id). Azure Foundry speaks the same Anthropic Messages API
  wire format and accepts the same `x-api-key` header, so no other change is
  needed. (Microsoft Entra ID token auth is not supported — API keys only.)

## Evidence sources (keyless, live, cached)

For the gene(s) entered on the Insights tab:

- **Europe PMC** — recent papers for `gene × cancer type` (title/year/PMID/link)
- **UniProt** — function, domains, active/binding sites, subcellular location, PDB structures
- **ChEMBL** — known drugs/mechanisms (mapped via UniProt accession) + druggability badge
- **STRING** — functional & physical confidence to the seed proteins; also folds
  into ranking and drives the Interactions tab (optional `STRING_API_KEY` used
  only as a rate-limit fallback)
- **BioMCP** (optional, off by default) — Reactome/KEGG pathway membership +
  Human Protein Atlas tissue expression/subcellular localization, via the
  external `biomcp` CLI (biomcp.org, MIT-licensed; install separately with
  `uv tool install biomcp-cli` — not bundled with this app, so this checkbox
  shows an honest "not installed" message rather than a fabricated result if
  you enable it without installing biomcp). Shelled out to as a subprocess
  (`app/agents/biomcp_client.py`); no new Python dependency.

All dossiers are shown as linked cards **and** fed to the LLM (grounded RAG): the
prompt forbids inventing drugs, PDB ids, or PMIDs, and uses STRING to distinguish
a **direct complex member** (high physical) from a **pathway-level** association.

## Notes

### Corporate VPN / TLS-inspecting proxy (`CERTIFICATE_VERIFY_FAILED`)

If a network tab fails with:

```
[SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed:
self-signed certificate in certificate chain
```

...while your **browser and `curl` reach the same URL fine**, you are almost
certainly behind a TLS-inspecting proxy (Cloudflare Zero Trust/WARP, Zscaler,
Netskope, many corporate VPNs). It re-signs HTTPS with a private root CA that
your OS trusts but Python's bundled `certifi` does not. It can hit **one host
and not others**, so partial failure is normal (e.g. Xena fails, ChEMBL works).

The app ships with **`truststore`** and enables it at startup (`app/net_trust.py`),
which makes Python use the OS trust store and resolves this automatically.
Certificate verification stays **fully enabled** — only the source of trusted
roots changes; the app never disables verification.

If you still hit it (e.g. running the modules outside the app):

```bash
# point Python at your proxy's root CA (export it from Keychain / cert store)
export SSL_CERT_FILE=/path/to/corporate-root-ca.pem
export REQUESTS_CA_BUNDLE=/path/to/corporate-root-ca.pem
```

...or simply disconnect the VPN/Zero Trust client. **Do not** disable TLS
verification to work around this.

### Other notes

- INO80 and SRCAP share 4 subunits (RUVBL1, RUVBL2, ACTL6A, ACTB); under "Both" they appear once.
- "Custom" seed source empties the seed-gene box so you can paste your own genes or upload
  a `.txt` file (newline/comma/whitespace separated) instead of using the built-in
  INO80/SRCAP protein lists.
- Small cohorts (age × rarer disease) trigger a min-cohort-size warning; correlations across few cell lines are unstable.
- The Interactions tab makes ~N STRING API calls on first load (cached after); it is opt-in.
- The Tissue Specificity tab needs outbound internet for Xena/GTEx (first call per session
  can take ~15–20s to fetch TCGA/TARGET sample metadata, cached after); tissue/tumor-cohort
  auto-suggestions are best-effort keyword matches, not authoritative — always verify.
- The Drug Targets tab needs outbound internet (MyGene.info/ChEMBL/DGIdb/Open Targets) and
  is the one place this app uses `requests` (all other evidence sources are dependency-free
  `urllib`). It's sequential per gene, so batches are bounded and cached per gene+cancer.
  Well-studied genes are slower (e.g. ALK: ~5.5 min, 3,297 compounds) than typical genes
  (~30s–2 min); PRISM enrichment (opt-in) adds ~100 MB + ~5 min. The original collaborator
  pipeline's per-drug ChEMBL name-lookup step (used to fetch cancer indications for
  DGIdb-only drugs) was hitting an intermittently-failing ChEMBL endpoint and was replaced
  with one Open Targets `drugAndClinicalCandidates` call per gene — verified 4–5× faster
  and fully reliable, with identical scoring results.
- Cite DepMap when using these data. See `ARCHITECTURE.md` for the full stack and diagrams.
- PDF export (Insights tab, both narrative types) uses `fpdf2` — pure Python, no system
  libraries (unlike e.g. weasyprint's Cairo/Pango dependency). `app/pdf_export.py` renders
  only the markdown subset this app's own narratives actually produce (headings, bold/
  italic, bullet lists, GFM pipe tables, and a bounded `==highlight==` span), not general
  markdown.
- Both narrative types support a small, app-owned emphasis vocabulary beyond **bold**/
  *italic*: `==highlighted text==` (not standard markdown) is parsed/rendered entirely by
  this app's own code — `ui.markdown_with_highlights()` on-screen, `pdf_export.py` in the
  PDF — never raw HTML/CSS from the LLM. The deterministic narrative highlights each
  gene's headline `selectivity_delta` clause; the LLM is instructed to use it sparingly,
  at most once or twice per gene, only around facts already present in the JSON.
- The Insights tab's "Include plots from other tabs" checkbox (off by default) reuses
  each tab's own cached last-run result to show the Co-expression volcano, Essentiality
  scatter, ranked-selectivity bar (requested gene(s) circled/recolored crimson via a new
  optional `highlight_genes` param on `viz/plots.py`'s functions), STRING physical heatmap,
  Tissue Specificity heatmap, Drug Targets bar/heatmap, and Model Predictions classifier
  bar — never triggers new computation from Insights. On-screen this needs nothing extra;
  embedding the same plots in the PDF export additionally uses the `kaleido` package to
  rasterize each Plotly figure to PNG. **`kaleido` v1.x no longer bundles a browser** — it
  needs a real Chrome/Chromium already installed on the machine running the app. If that
  prerequisite is missing, `_fig_to_png_bytes()` returns `None` and the PDF gets one clear
  note ("kaleido needs a real Chrome/Chromium...") instead of a crash; the on-screen plots
  are unaffected either way.
