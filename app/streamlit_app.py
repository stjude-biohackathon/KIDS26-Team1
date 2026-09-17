"""INO80 DepMap Explorer — Streamlit orchestrator.

Run:  uv run streamlit run app/streamlit_app.py
Reads only harmonized parquets under data/processed/ (built from the read-only
user_data/ CSVs by app/build_harmonized.py). Never writes user_data/.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

APP = Path(__file__).resolve().parent
sys.path.insert(0, str(APP))
RAW = APP.parent / "user_data"

# Route TLS through the OS trust store BEFORE importing anything that may open
# an HTTPS connection (Xena/GTEx, STRING, ChEMBL). Needed behind TLS-inspecting
# proxies whose private root CA lives in the OS store but not in certifi.
# Verification is never disabled -- see app/net_trust.py.
from net_trust import enable_os_trust_store  # noqa: E402

_TRUST_OK, _TRUST_MSG = enable_os_trust_store()

from agents.base import AnalysisContext, SeedSelection, Params  # noqa: E402
from agents.data_agent import DataAgent, _load_all  # noqa: E402
from agents.correlation_agents import CoexpressionAgent  # noqa: E402
from agents.essentiality_agent import EssentialityAgent  # noqa: E402
from agents.synthesis_agent import SynthesisAgent  # noqa: E402
from core.cohort import CohortSpec, disease_options  # noqa: E402
from core.drilldown import cohort_gene_vectors  # noqa: E402
from agents import knowledge  # noqa: E402
from agents import tissue_specificity  # noqa: E402
from agents import gtex_agent  # noqa: E402
from agents import drug_targets  # noqa: E402
from agents import protein_predictions  # noqa: E402
from viz import plots  # noqa: E402
import ui  # noqa: E402
import pdf_export  # noqa: E402

ui.register_plotly_templates()
_ICON = APP.parent / "assets" / "scrap_ai_icon.png"
st.set_page_config(page_title="SCRAP-AI · DepMap Explorer", layout="wide",
                   page_icon=str(_ICON) if _ICON.exists() else "\U0001f9ec",
                   initial_sidebar_state="expanded")


# ---------- seed loading ----------
@st.cache_data
def load_seed_file(name: str) -> list[str]:
    return [g.strip() for g in (RAW / name).read_text().splitlines() if g.strip()]


def resolve_seeds(source: str) -> list[str]:
    ino80 = load_seed_file("ino80.proteins.txt")
    srcap = load_seed_file("srcap.proteins.txt")
    if source == "INO80":
        return ino80
    if source == "SRCAP":
        return srcap
    if source == "Custom":
        return []  # user pastes/uploads their own list; nothing pre-filled
    return list(dict.fromkeys(ino80 + srcap))  # Both, de-duplicated union


def _sync_default_genes(key: str, default_text: str) -> None:
    """Keep a gene-list text_input's auto-populated default in sync with the
    current top-ranked genes across reruns.

    Streamlit only honors a keyed widget's `value=` the first time it's
    created; afterwards the widget is fully controlled by session_state, so a
    changed ranking (from different thresholds/params) would otherwise never
    show up. This re-seeds session_state[key] ONLY when the underlying
    computed default has actually changed since we last set it -- so it
    updates when the ranked-gene set changes, but doesn't clobber a manual
    edit the user made while the ranking stayed the same. Call BEFORE creating
    the widget (Streamlit forbids mutating session_state after instantiation).
    """
    sig_key = f"_{key}_autoval"
    if st.session_state.get(sig_key) != default_text:
        st.session_state[key] = default_text
        st.session_state[sig_key] = default_text


@st.cache_data(show_spinner=False)
def meta_frame() -> pd.DataFrame:
    return _load_all()[2]


# ---------- cached pipeline ----------
@st.cache_data(show_spinner="Running agents…")
def run_pipeline(age: str, disease: str, include_fetus: bool,
                 source: str, genes: tuple[str, ...], r_min: float, q_max: float,
                 min_n: int, var_pctile: float, positive_only: bool,
                 method: str, population: str,
                 expr_tpm: float, ess_thresh: float, seed_on_frac: float, ce_penalty: float,
                 pan_exclude_cohort: bool):
    spec = CohortSpec(age_group=age, disease=disease, include_fetus=include_fetus)
    seeds = SeedSelection(source=source, genes=list(genes))
    params = Params(r_min=r_min, q_max=q_max, min_cohort_n=min_n,
                    var_pctile=var_pctile, positive_only=positive_only,
                    method=method, population=population,
                    expr_tpm=expr_tpm, ess_thresh=ess_thresh, seed_on_frac=seed_on_frac,
                    ce_penalty=ce_penalty, pan_exclude_cohort=pan_exclude_cohort)
    ctx = AnalysisContext(spec=spec, seeds=seeds, params=params)
    for agent in (DataAgent(), CoexpressionAgent(), EssentialityAgent(), SynthesisAgent()):
        ctx = agent.run(ctx)
    return ctx


# ---------- sidebar ----------
st.sidebar.markdown("### \U0001f9ec  SCRAP-AI")
_theme = st.sidebar.radio("Theme", ["Light", "Dark"],
                         index=0 if st.session_state.get("ui_theme", "Light") == "Light" else 1,
                         horizontal=True, key="ui_theme")
ui.apply_theme()
ui.force_sidebar_expanded_on_load()
ui.force_dark_mode_text_colors()
ui.sidebar_credits("Team1 \u00b7 St Jude Bio Hackathon 2026",
                  github_url="https://github.com/stjude-biohackathon/KIDS26-Team1")

ui.sidebar_section("\U0001f9ec  Cohort")
age = st.sidebar.radio("Age group", ["Pediatric", "Adult"], horizontal=True)
include_fetus = False
if age == "Pediatric":
    include_fetus = st.sidebar.checkbox("Include Fetus", value=False)

age_cats = ("Adult",) if age == "Adult" else (("Pediatric", "Fetus") if include_fetus else ("Pediatric",))
counts = disease_options(meta_frame(), age_cats)
dis_labels = [f"{d}  (n={n})" for d, n in counts.items()]
dis_map = {f"{d}  (n={n})": d for d, n in counts.items()}
disease_label = st.sidebar.selectbox("Cancer type (OncotreePrimaryDisease)", dis_labels)
disease = dis_map[disease_label]

ui.sidebar_section("\U0001f331  Seeds")
source = st.sidebar.radio("Seed source", ["INO80", "SRCAP", "Both", "Custom"], index=2,
                          horizontal=True,
                          help="'Custom' clears the seed list below so you can paste your "
                               "own genes or upload a text file, instead of using the "
                               "built-in INO80/SRCAP proteins.")

if source == "Custom":
    uploaded_seed_file = st.sidebar.file_uploader(
        "Upload gene list (.txt)", type=["txt"],
        help="One gene symbol per line, or comma/whitespace separated. Populates the box "
             "below — review and edit before running.")
    if uploaded_seed_file is not None:
        _raw = uploaded_seed_file.getvalue().decode("utf-8", errors="ignore")
        _upload_genes = [g.strip() for g in _raw.replace(",", " ").split() if g.strip()]
        _seed_default_text = "\n".join(_upload_genes)
    else:
        _seed_default_text = ""  # blank canvas until the user pastes or uploads something
else:
    _seed_default_text = "\n".join(resolve_seeds(source))

_sync_default_genes("seed_genes_text", _seed_default_text)
edited = st.sidebar.text_area("Seed genes (editable)", key="seed_genes_text", height=140)
genes = [g.strip() for g in edited.splitlines() if g.strip()]
st.sidebar.caption(f"{len(genes)} seed genes · source: {source}")

ui.sidebar_section("\u2699\ufe0f  Method & thresholds")
with st.sidebar.expander("Method & population", expanded=False):
    method = st.radio("Correlation method", ["pearson", "spearman"], horizontal=True)
    population_label = st.radio(
        "Correlation population",
        ["Cohort (age + disease)", "Pan-cancer (all analyzable lines)"],
    )
    population = "pan_cancer" if population_label.startswith("Pan") else "cohort"

with st.sidebar.expander("Co-expression thresholds", expanded=False):
    r_min = st.slider("min |r|", 0.0, 0.9, 0.30, 0.05,
                      help="Aggregate mean|r| vs the seed set; cohort values rarely exceed ~0.3.")
    q_max = st.select_slider("max FDR q", options=[0.001, 0.01, 0.05, 0.1, 0.25], value=0.05)
    min_n = st.number_input("min cohort size", 3, 100, 15)
    var_pctile = st.slider("variance filter percentile", 0.0, 75.0, 25.0, 5.0)
    positive_only = st.checkbox("positive correlations only", value=True)

with st.sidebar.expander("Expression & essentiality", expanded=True):
    expr_tpm = st.slider("expressed if TPM >", 0.0, 10.0, 3.0, 0.5,
                         help="Per cell line; matrix is log2(TPM+1).")
    ess_thresh = st.slider("essential if Chronos ≤", -2.0, 0.0, -0.5, 0.05,
                           help="More negative = stronger dependency.")
    seed_on_frac = st.slider("seed set 'on' if ≥ this frac of seeds expressed",
                             0.0, 1.0, 0.5, 0.05)
    ce_penalty = st.slider("common-essential penalty (ranking)", 0.0, 1.0, 0.25, 0.05,
                           help="Multiplier on DepMap common-essential genes in the "
                                "selectivity score. 0 = exclude them from the top, "
                                "1 = no penalty.")
    pan_exclude_cohort = st.checkbox("Pan-cancer baseline excludes the cohort", value=False,
                                     help="selectivity_delta = cohort essential% − pan-cancer "
                                          "essential%. On: pan-cancer = all analyzable lines "
                                          "EXCLUDING the selected cohort (purer in- vs "
                                          "out-of-cohort contrast). Off: includes the cohort.")

show_ce = st.sidebar.checkbox("Show common-essential genes", value=False,
                             help="Uncheck to hide DepMap common-essential (pan-essential) "
                                  "genes everywhere — tables, plots, and heatmap.")

ui.sidebar_section("\U0001f517  STRING interactions")
fold_string = st.sidebar.checkbox("Fold STRING into ranking", value=True,
                                  help="Re-rank the top genes using STRING physical + "
                                       "functional association to the seed complex "
                                       "(one API call per gene; opt-in).")
string_topn = 25
sw_sel, sw_phys, sw_func = 0.5, 0.35, 0.15
string_rank_cut = 0.4
if fold_string:
    with st.sidebar.expander("STRING weights", expanded=True):
        string_topn = st.slider("top genes to re-rank", 10, 50, 25, 5)
        st.caption("Weights auto-normalize to sum 1 — string_score ∈ [0,1].")
        sw_sel = st.slider("selectivity weight", 0.0, 1.0, 0.5, 0.05)
        sw_phys = st.slider("physical weight", 0.0, 1.0, 0.35, 0.05)
        sw_func = st.slider("functional weight", 0.0, 1.0, 0.15, 0.05)
        string_rank_cut = st.slider("edge cutoff", 0.15, 0.9, 0.4, 0.05)

run = st.sidebar.button("Run analysis", type="primary", use_container_width=True)


# ---------- main ----------
ui.header("AI-assisted Selective Cancer-dependency Ranking of Associated Proteins",
          "INO80 / SRCAP dependency explorer",
          "By Team 1 (St Jude Bio Hackathon 2026)")

if "has_run" not in st.session_state:
    st.session_state["has_run"] = True  # auto-run once on first page load, using defaults
if run:
    st.session_state["has_run"] = True
if not st.session_state.get("has_run"):
    st.info("\U0001f44b  Set the cohort and seed source in the sidebar, then click "
            "**Run analysis** to explore selective dependencies of the INO80/SRCAP complex.")
    st.stop()

with st.spinner("Running the 4-agent pipeline…"):
    ctx = run_pipeline(age, disease, include_fetus, source, tuple(genes),
                       r_min, q_max, int(min_n), var_pctile, positive_only, method, population,
                       expr_tpm, ess_thresh, seed_on_frac, ce_penalty, pan_exclude_cohort)

# Single source of truth for "which ranked genes is the user actually being shown".
# The Essentiality tab hides common-essential genes when show_ce is unchecked; the
# auto-populated gene lists (Tissue Specificity / Drug Targets / Insights) must use
# the SAME view, otherwise hidden genes leak into those boxes and the user sees
# genes that appear nowhere in the ranked table.
if (ctx.ranked is not None and not ctx.ranked.empty
        and "common_essential" in ctx.ranked.columns and not show_ce):
    ranked_visible = ctx.ranked[~ctx.ranked["common_essential"]]
else:
    ranked_visible = ctx.ranked

# cohort metric cards
ui.metric_cards([
    (str(ctx.ess_denom), "cohort lines (both data)"),
    ("0" if ctx.coexpr is None else str(len(ctx.coexpr)), "co-expr candidates"),
    ("0" if ranked_visible is None else str(len(ranked_visible)),
     "ranked genes" if show_ce else "ranked genes (shown)"),
    (str(ctx.corr_n), "lines correlated"),
])

# status chips
_chips = [
    (f"{age} · {disease}", "info"),
    (f"seeds: {source} ({len(genes)})", "info"),
    (f"{method} · {population.replace('_', '-')}", ""),
    (f"TPM&gt;{expr_tpm:g} · Chronos≤{ess_thresh:g}", ""),
]
if ctx.ess_denom and ctx.ess_denom < int(min_n):
    _chips.append((f"⚠ small cohort (n={ctx.ess_denom})", "warn"))
if ctx.seeds_missing:
    _chips.append((f"{len(ctx.seeds_missing)} seed(s) missing", "warn"))
ui.chips(_chips)
for w in ctx.warnings:
    st.warning(w)

tab_ov, tab_cx, tab_ess, tab_hm, tab_int, tab_gtex, tab_drug, tab_clf, tab_ins = st.tabs(
    ["Overview", "Co-expression", "Essentiality & ranking", "Heatmap",
     "Interactions", "Tissue Specificity", "Drug Targets", "Model Predictions", "Insights"])

with tab_ov:
    st.subheader(f"{age} · {disease} · seeds={source} — {ctx.ess_denom} cohort lines with both data")
    st.markdown(ctx.narrative or "_Run to see synthesis._")


@st.cache_data(show_spinner="Fetching STRING & re-ranking…")
def string_rerank_cached(genes_index: tuple, sel_scores: tuple, seeds: tuple,
                         top_n: int, w_sel: float, w_phys: float, w_func: float, cutoff: float):
    import pandas as _pd
    from agents.string_rank import string_rerank
    df = _pd.DataFrame({"selectivity_score": list(sel_scores)}, index=list(genes_index))
    return string_rerank(df, seeds, top_n, w_sel, w_phys, w_func, cutoff)


@st.cache_data(show_spinner=False)
def drilldown_vectors(age, disease, include_fetus, gene_a, gene_b, modality):
    spec = CohortSpec(age_group=age, disease=disease, include_fetus=include_fetus)
    return cohort_gene_vectors(spec, gene_a, gene_b, modality)


@st.cache_data(show_spinner=False)
def _fetch_dossiers(genes: tuple, disease: str, n_papers: int) -> dict:
    return {g: knowledge.dossier(g, disease, n_papers) for g in genes}


@st.cache_data(show_spinner=False)
def _fetch_string(genes: tuple, seeds: tuple, cutoff: float) -> dict:
    from agents import string_db
    return {g: string_db.string_scores(g, seeds, cutoff) for g in genes}


@st.cache_data(show_spinner=False)
def _fetch_biomcp(genes: tuple) -> dict:
    from agents import biomcp_client
    return {g: biomcp_client.gene_pathways_hpa(g) for g in genes}


def _render_biomcp(genes, prefetched=None):
    from agents import biomcp_client
    data = prefetched if prefetched is not None else _fetch_biomcp(tuple(genes))
    st.divider()
    st.subheader("BioMCP: pathways + tissue atlas")
    if not biomcp_client.is_available() and all(v is None for v in data.values()):
        st.info("The optional `biomcp` CLI isn't installed (or isn't on PATH) — "
                "install with `uv tool install biomcp-cli` (biomcp.org) to enable "
                "this section. Nothing was fabricated in its place.")
        return
    for g in genes:
        d = data.get(g)
        with st.expander(f"{g} — BioMCP", expanded=False):
            if d is None:
                st.caption("Not available for this gene (unresolved symbol, or the "
                           "underlying APIs were unreachable).")
                continue
            pathways = d.get("pathways") or []
            if pathways:
                st.markdown(f"**Pathways ({len(pathways)}):**")
                st.dataframe(
                    pd.DataFrame(pathways)[["source", "id", "name"]],
                    use_container_width=True, height=min(38 * (len(pathways) + 1), 300),
                )
            else:
                st.caption("No Reactome/KEGG pathway membership found.")
            hpa = d.get("hpa")
            if hpa:
                st.markdown(f"**HPA protein summary:** {hpa.get('protein_summary', '—')}")
                st.markdown(f"**HPA RNA summary:** {hpa.get('rna_summary', '—')} "
                            f"(reliability: {hpa.get('reliability', '—')})")
                loc = ", ".join(hpa.get("subcellular_main_location", []) or [])
                if loc:
                    st.markdown(f"**Subcellular location:** {loc}")
                tissues = hpa.get("tissues") or []
                high = [t["tissue"] for t in tissues if t.get("level") == "High"]
                if high:
                    st.caption(f"High expression in {len(high)} tissue(s): "
                               f"{', '.join(high[:10])}{'…' if len(high) > 10 else ''}")
            else:
                st.caption("No Human Protein Atlas data found.")


def _render_string(genes, seeds, cutoff, prefetched=None):
    from agents import string_db
    data = prefetched or _fetch_string(tuple(genes), tuple(seeds), cutoff)
    st.divider()
    st.subheader("STRING association with the seed complex")
    st.caption(f"Functional & physical confidence (0–1) vs the {len(seeds)} seed proteins; "
               f"cutoff {cutoff:g}. Physical ≈ direct complex membership; functional ≈ pathway.")
    for g in genes:
        sd = data.get(g)
        if not sd:
            continue
        fn, ph = sd["functional"], sd["physical"]
        with st.expander(f"{g} — STRING", expanded=False):
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Functional max", f"{fn['max']:.3f}", help=f"best seed: {fn['best_seed'] or '—'}")
            m2.metric("Functional mean", f"{fn['mean']:.3f}")
            m3.metric("Physical max", f"{ph['max']:.3f}", help=f"best seed: {ph['best_seed'] or '—'}")
            m4.metric("Physical mean", f"{ph['mean']:.3f}")
            interp = ("strong physical link → likely direct complex member"
                      if ph["max"] >= 0.7 else
                      ("functional link only → pathway-level association"
                       if fn["max"] >= cutoff else "no STRING evidence above cutoff"))
            st.caption(f"**{fn['n_seeds_ge_cutoff']}** seeds functional ≥{cutoff:g} · "
                       f"**{ph['n_seeds_ge_cutoff']}** physical ≥{cutoff:g} — _{interp}_")
            if fn["linked_seeds"]:
                top = list(fn["linked_seeds"].items())[:8]
                st.markdown("Top functional links: " + ", ".join(f"{k} ({v:.2f})" for k, v in top))
            img = string_db.network_image(g, tuple(seeds[:8]))
            if img:
                st.image(img, caption=f"STRING network: {g} + seeds", width=440)
            st.markdown(f"[Open in STRING]({sd['string_url']})")


def _render_dossiers(genes, disease, n_papers, prefetched=None):
    doss = prefetched or _fetch_dossiers(tuple(genes), disease, n_papers)
    st.divider()
    st.subheader("Literature & drug-target annotations")
    for g in genes:
        d = doss.get(g)
        if not d:
            continue
        up = d.get("uniprot") or {}
        ch = d.get("chembl") or {}
        with st.expander(f"{g} — {up.get('protein_name','(no UniProt match)')}", expanded=True):
            st.markdown(ui.druggability_badge(ch.get("n_mechanisms", 0),
                                              (up or {}).get("n_pdb", 0)),
                        unsafe_allow_html=True)
            if up:
                links = []
                if up.get("uniprot_url"):
                    links.append(f"[UniProt {up.get('accession')}]({up['uniprot_url']})")
                if ch.get("chembl_url"):
                    links.append(f"[ChEMBL]({ch['chembl_url']})")
                if up.get("pdb_ids"):
                    pid = up["pdb_ids"][0]
                    links.append(f"[PDB {pid}](https://www.rcsb.org/structure/{pid})")
                if links:
                    st.markdown(" · ".join(links))
                if up.get("function"):
                    st.markdown(f"**Function.** {up['function'][:600]}")
                cols = st.columns(2)
                with cols[0]:
                    if up.get("domains"):
                        st.markdown("**Domains:** " + "; ".join(up["domains"]))
                    if up.get("active_sites"):
                        st.markdown("**Active sites:** " + "; ".join(up["active_sites"]))
                    if up.get("binding_sites"):
                        st.markdown("**Binding sites:** " + "; ".join(up["binding_sites"]))
                with cols[1]:
                    if up.get("subcellular_location"):
                        st.markdown("**Location:** " + ", ".join(up["subcellular_location"]))
                    st.markdown(f"**PDB structures:** {up.get('n_pdb',0)}")
                    if up.get("keywords"):
                        st.caption("Keywords: " + ", ".join(up["keywords"]))
            if ch.get("n_mechanisms", 0) > 0:
                st.markdown(f"**Known drugs / mechanisms ({ch['n_mechanisms']}):**")
                st.dataframe(pd.DataFrame(ch["mechanisms"]), use_container_width=True,
                             hide_index=True)
            else:
                st.caption("No approved-drug mechanism in ChEMBL.")
            lit = d.get("literature", [])
            if lit:
                st.markdown(f"**Literature ({g} × {disease}):**")
                for p in lit:
                    cite = f"{p['title']} — *{p['journal']}* {p['year']}"
                    if p.get("url"):
                        cite += f" · [{'PMID '+p['pmid'] if p['pmid'] else 'link'}]({p['url']})"
                    st.markdown(f"- {cite}")
            else:
                st.caption(f"No Europe PMC hits for {g} × {disease}.")


def _coexpr_tab(agg, label="Co-expression", modality="expr"):
    if agg is None or agg.empty:
        st.info("No co-expressed candidates (cohort too small or no seeds found).")
        return
    ev = st.plotly_chart(
        plots.volcano(agg, r_min, q_max, f"{label} volcano — {disease}/{age}", positive_only),
        use_container_width=True, key=f"volcano_{modality}", on_select="rerun",
    )
    picked = None
    try:
        pts = ev.selection["points"]  # type: ignore[attr-defined]
        if pts:
            picked = pts[0]["customdata"][0]
    except (AttributeError, KeyError, IndexError, TypeError):
        picked = None

    seeds_found = ctx.seeds_found_expr
    cc1, cc2 = st.columns([1, 1])
    default_gene = picked or (agg.index[0] if len(agg) else None)
    gene_sel = cc1.selectbox(f"Candidate gene ({label})", list(agg.index),
                             index=(list(agg.index).index(default_gene) if default_gene in agg.index else 0),
                             key=f"cand_{modality}")
    best_seed = agg.loc[gene_sel, "best_seed"] if gene_sel in agg.index else (seeds_found[0] if seeds_found else None)
    seed_sel = cc2.selectbox("Seed gene", seeds_found,
                             index=(seeds_found.index(best_seed) if best_seed in seeds_found else 0) if seeds_found else 0,
                             key=f"seed_{modality}")
    if picked:
        st.caption(f"Clicked **{picked}** → scatter below.")
    if gene_sel and seed_sel:
        vec = drilldown_vectors(age, disease, include_fetus, seed_sel, gene_sel, modality)
        if vec is None or vec.empty:
            st.info("No paired data for this gene/seed in the cohort.")
        else:
            unit = "log2(TPM+1)"
            st.plotly_chart(
                plots.drilldown_scatter(vec["x"], vec["y"], vec["cell_line"],
                                        f"{seed_sel} ({unit})", f"{gene_sel} ({unit})", label),
                use_container_width=True, key=f"dd_{modality}")
    st.dataframe(agg.round(4), use_container_width=True, height=360)
    st.download_button(f"Download {label} CSV", agg.to_csv().encode(),
                       file_name=f"{label.lower()}_{disease}_{age}_{source}.csv".replace(" ", "_"),
                       mime="text/csv")


with tab_cx:
    st.caption("Genes whose expression correlates with the seed complex across the "
               "correlation population. These are the candidates carried into ranking.")
    _coexpr_tab(ctx.coexpr)

with tab_ess:
    ranked = ctx.ranked
    if ranked is None or ranked.empty:
        st.info("No ranked genes. Ensure co-expressed candidates exist and the cohort "
                "has cell lines with both expression and CRISPR data; try relaxing the "
                "TPM or |r| thresholds.")
    else:
        view = ranked_visible
        st.markdown(f"**{len(view)} genes** (of {len(ranked)} co-expressed) ranked by a "
                    "**selectivity-weighted score** (joint % × cohort-vs-pan-cancer "
                    f"essentiality × common-essential penalty {ce_penalty:g}×), across "
                    f"{ctx.ess_denom} cohort cell lines."
                    + ("" if show_ce else "  — _common-essential genes hidden._"))
        cola, colb = st.columns(2)
        with cola:
            st.plotly_chart(plots.essentiality_scatter(
                view, "% co-expressed vs % essential"), use_container_width=True)
        with colb:
            top_k = st.slider("Top genes (bar)", 5, 40, 20, 5)
            st.plotly_chart(plots.ranked_bar(view, top_k, "Top genes by selectivity score"),
                            use_container_width=True)
        show_cols = ["selectivity_score", "pct_joint", "selectivity_delta",
                     "pct_essential", "pct_essential_pan", "pct_coexpressed",
                     "pct_essential_given_coexpr", "best_seed", "common_essential"]
        show_cols = [c for c in show_cols if c in view.columns]
        st.dataframe(view[show_cols].round(3), use_container_width=True, height=420)
        st.download_button("Download ranked CSV", view.to_csv().encode(),
                           file_name=f"ranked_{disease}_{age}_{source}.csv".replace(" ", "_"),
                           mime="text/csv")

with tab_hm:
    st.caption("Selective genes only (DepMap common-essential excluded), ordered by "
               "selectivity-weighted score. Rows = candidates, columns = seed genes, colour = r.")
    top_k = st.slider("Top selective genes", 10, 60, 30, 5, key="hm_k")
    if ctx.ranked is None or ctx.ranked.empty:
        st.info("No ranked genes to show.")
    else:
        selective_order = list(ctx.ranked[~ctx.ranked["common_essential"]].index)
        if not selective_order:
            st.info("No selective (non-common-essential) genes at current thresholds.")
        else:
            st.plotly_chart(plots.heatmap(ctx.coexpr_matrix, top_k, ctx.coexpr,
                            "Selective candidate × seed r (ordered by selectivity score)",
                            order=selective_order), use_container_width=True)

with tab_int:
    st.subheader("STRING interactions with the seed complex")
    if not fold_string:
        st.info("Enable **Fold STRING into ranking** in the sidebar to fetch STRING "
                "physical & functional associations and re-rank the top genes.")
    elif ranked_visible is None or ranked_visible.empty:
        st.info("No ranked genes to score.")
    else:
        _top = ranked_visible.head(string_topn)
        rr, phys_m, func_m, edges = string_rerank_cached(
            tuple(_top.index), tuple(_top["selectivity_score"]), tuple(genes),
            string_topn, sw_sel, sw_phys, sw_func, string_rank_cut)
        wts = rr.attrs.get("weights", {})
        st.caption(f"Top {len(rr)} genes re-ranked by a normalized weighted-sum — "
                   f"string_score = {wts.get('selectivity','?')}·sel_norm + "
                   f"{wts.get('physical','?')}·physical + {wts.get('functional','?')}·functional "
                   f"(all ∈ [0,1]; weights normalized). Edge cutoff {string_rank_cut:g}: STRING "
                   "links below it count as 0 in the score and are hidden from the network.")
        kind = st.radio("Network / heatmap layer", ["physical", "functional"],
                        horizontal=True, key="int_kind")
        colA, colB = st.columns([1, 1])
        with colA:
            mat = phys_m if kind == "physical" else func_m
            st.plotly_chart(plots.string_heatmap(mat, f"STRING {kind} — candidate × seed"),
                            use_container_width=True)
        with colB:
            node_size = {g: float(rr.loc[g, "selectivity_score"]) for g in rr.index}
            st.plotly_chart(plots.string_network(edges, node_size, kind,
                            "Candidate–seed interactions"), use_container_width=True)
        cols = ["string_score", "sel_norm", "string_phys", "string_func",
                "selectivity_score", "string_best_phys_seed", "pct_joint", "common_essential"]
        cols = [c for c in cols if c in rr.columns]
        st.dataframe(rr[cols].round(3), use_container_width=True, height=420)
        st.download_button("Download STRING-ranked CSV", rr.to_csv().encode(),
                           file_name=f"string_ranked_{disease}_{age}_{source}.csv".replace(" ", "_"),
                           mime="text/csv")


@st.cache_data(show_spinner=False)
def _lineage_for_disease(disease_name: str) -> str:
    m = meta_frame()
    sub = m[m["OncotreePrimaryDisease"] == disease_name]
    if sub.empty or "OncotreeLineage" not in sub.columns:
        return ""
    vc = sub["OncotreeLineage"].value_counts()
    return vc.index[0] if len(vc) else ""


@st.cache_data(show_spinner="Running tissue-specificity triage…")
def tissue_triage_cached(genes: tuple, tissues: tuple, targets: tuple, cutoff: float,
                        pancancer_frac: float, target_median: float, target_frac: float,
                        gate_mode: str):
    return tissue_specificity.triage(list(genes), list(tissues), list(targets),
                                     normal_cutoff=cutoff, pancancer_fraction=pancancer_frac,
                                     target_median_tpm=target_median,
                                     target_fraction_expressed=target_frac,
                                     normal_gate_mode=gate_mode)


with tab_gtex:
    st.subheader("Tissue specificity — GTEx normal tissue + TCGA/TARGET pan-cancer")
    st.caption("For gene(s) you enter: silent in the normal tissue(s) you pick, low in most "
               "other TCGA cancer types, and highly expressed in the matched tumor cohort. "
               "On-demand (Xena network queries); not part of the ranking pipeline. "
               "Ported and generalized from a collaborator's Neuroblastoma-specific pipeline "
               "(see user_data/gtex_other_expression/README.md); fetal/CELLxGENE confirmation "
               "layer not integrated — run that script standalone if needed.")

    _gtex_default = ", ".join(list(ranked_visible.index[:10])) if (ranked_visible is not None and not ranked_visible.empty) else ""
    _sync_default_genes("gtex_genes", _gtex_default)
    gtex_genes_txt = st.text_input("Gene(s), comma or space separated", key="gtex_genes",
                                   help="Auto-populated with the current top-10 ranked genes; "
                                        "updates when parameters change the ranking, unless "
                                        "you've edited it since the last change.")
    gtex_genes = [g.strip().upper() for g in gtex_genes_txt.replace(",", " ").split() if g.strip()]

    lineage_guess = _lineage_for_disease(disease)
    try:
        tissue_suggestions = tissue_specificity.suggest_gtex_tissues(lineage_guess)
        all_tissues = gtex_agent.gtex_tissues()
        all_targets = tissue_specificity.target_label_options()
        target_suggestions = (tissue_specificity.suggest_target_labels(disease)
                              or tissue_specificity.suggest_target_labels(lineage_guess))
        options_ok = True
    except Exception as e:  # noqa: BLE001
        _msg = str(e)
        if "CERTIFICATE_VERIFY_FAILED" in _msg or "SSL" in type(e).__name__:
            st.error(
                "**TLS certificate verification failed** contacting the Xena hub "
                f"({_msg}).\n\nThis usually means a TLS-inspecting proxy or VPN "
                "(Cloudflare Zero Trust/WARP, Zscaler, Netskope, corporate VPN) is "
                "re-signing HTTPS traffic with a private root CA that your OS trusts "
                "but Python's bundled `certifi` does not — your browser and `curl` "
                "will work fine, which is the giveaway.\n\n"
                "**Fixes:** install `truststore` (`uv add truststore`) so Python uses "
                "the OS trust store; or export the proxy's root CA and set "
                "`SSL_CERT_FILE=/path/to/ca.pem`; or disconnect the VPN/Zero Trust client.")
        else:
            st.error(f"Could not load GTEx/Xena tissue & tumor-type options ({e}). "
                     "This tab needs outbound internet on first use (GTEx file + a live Xena "
                     "lookup). Retry after checking connectivity.")
        tissue_suggestions, all_tissues, all_targets, target_suggestions = [], [], [], []
        options_ok = False
    gcol1, gcol2 = st.columns([1, 1])
    with gcol1:
        chosen_tissues = st.multiselect(
            f"Normal GTEx tissue(s) — suggested for '{lineage_guess or disease}'",
            options=all_tissues, default=tissue_suggestions,
            help="Best-effort suggestion from the cancer lineage; not authoritative — add/remove freely.")
    with gcol2:
        chosen_targets = st.multiselect(
            f"Matched tumor cohort (TCGA/TARGET) — suggested for '{disease}'",
            options=all_targets, default=target_suggestions,
            help="Best-effort keyword match against live Xena TCGA/TARGET labels; verify/adjust.")

    with st.expander("Thresholds", expanded=False):
        tc1, tc2 = st.columns(2)
        normal_cutoff = tc1.slider("normal/pan-cancer TPM cutoff", 1.0, 20.0, 5.0, 0.5)
        pancancer_frac = tc2.slider("pan-cancer low fraction", 0.0, 1.0, 0.20, 0.05,
                                    help="Minimum fraction of OTHER TCGA/TARGET cancer types "
                                         "in which the gene must be below the TPM cutoff to "
                                         "pass the pan-cancer gate. 1.0 = silent in every other "
                                         "cancer type (most selective); 0.0 = gate effectively "
                                         "off. Also multiplies into selectivity_score.")
        target_median = tc1.slider("target median TPM (positive gate)", 1.0, 50.0, 10.0, 1.0)
        target_frac = tc2.slider("target fraction expressed (positive gate)", 0.1, 1.0, 0.50, 0.05)
        gate_mode = st.radio("Normal-tissue aggregation", ["median", "all"], horizontal=True,
                             help="'median': tolerant of one atypical subregion (matches the "
                                  "original pipeline). 'all': every selected tissue must be below "
                                  "cutoff (stricter).")

    run_gtex = st.button("Run tissue-specificity triage", type="primary", disabled=not options_ok)
    gtex_table = None
    if run_gtex:
        if not gtex_genes:
            st.warning("Enter at least one gene.")
        elif not chosen_tissues:
            st.warning("Pick at least one normal GTEx tissue.")
        elif not chosen_targets:
            st.warning("Pick at least one matched tumor cohort label.")
        else:
            try:
                gtex_table, gtex_heat = tissue_triage_cached(
                    tuple(gtex_genes), tuple(chosen_tissues), tuple(chosen_targets),
                    normal_cutoff, pancancer_frac, target_median, target_frac, gate_mode)
                st.session_state["gtex_last"] = (gtex_table, gtex_heat, chosen_tissues, chosen_targets)
            except (OSError, TimeoutError) as e:
                st.error(f"Xena/GTEx network error: {e}. This needs outbound internet "
                         "access; retry, or reduce the gene/sample count.")
            except Exception as e:  # noqa: BLE001
                st.error(f"Tissue-specificity query failed: {e}. This is likely a data/code "
                         "issue (e.g. an unresolved gene symbol), not a connectivity problem "
                         "— please report it. Try removing the least-common gene(s) from the "
                         "list and retry.")

    if gtex_table is None and "gtex_last" in st.session_state:
        gtex_table, gtex_heat, chosen_tissues, chosen_targets = st.session_state["gtex_last"]

    if gtex_table is not None and not gtex_table.empty:
        n_pass = int(gtex_table["pass_all"].sum())
        tier_counts = gtex_table["specificity_tier"].value_counts()
        tier_summary = " · ".join(
            f"{tier_counts.get(t, 0)} {t}" for t in ("high", "medium", "low", "fail"))
        st.markdown(f"**{n_pass}/{len(gtex_table)} genes pass all gates** "
                    f"(silent in {', '.join(chosen_tissues[:3])}{'…' if len(chosen_tissues) > 3 else ''}; "
                    f"low in ≥{pancancer_frac:.0%} of other TCGA types; "
                    f"high in {', '.join(chosen_targets)}).")
        st.caption(f"Graded breakdown (by number of gates passed, 0–3): {tier_summary}. "
                   "'high' = pass_all; 'fail'/'low'/'medium' distinguish genes that pass_all "
                   "alone would otherwise lump together as a single 'False'.")
        show_cols = ["rank", "normal_max", "pass_normal", "pancancer_low_fraction", "pass_pancancer",
                    "target_median_tpm", "target_fraction_expressed", "pass_target", "pass_all",
                    "n_gates_passed", "specificity_tier", "selectivity_score"]
        st.dataframe(gtex_table[show_cols].round(3), use_container_width=True, height=360)
        st.plotly_chart(plots.tissue_heatmap(gtex_heat, "Expression: normal tissue → TCGA → target cohort"),
                        use_container_width=True)
        st.download_button("Download tissue-specificity CSV", gtex_table.to_csv().encode(),
                           file_name=f"tissue_specificity_{disease}.csv".replace(" ", "_"),
                           mime="text/csv")


@st.cache_data(show_spinner=False)
def _resolve_cancer_preview(cancer_name: str):
    try:
        resolver = drug_targets.CancerResolver()
        return resolver.resolve(cancer_name)
    except Exception:
        return None


def _run_drug_targets(genes, cancer_type, min_score, max_drugs_per_gene, use_prism):
    import tempfile
    prog = st.progress(0.0, text="Starting…")
    status = st.empty()

    def _cb(i, total, gene):
        prog.progress(i / max(total, 1), text=f"[{i}/{total}] Querying {gene}…")
        status.caption(f"Resolving {gene} against ChEMBL / DGIdb / Open Targets…")

    with tempfile.TemporaryDirectory() as tmp:
        df = drug_targets.run_pipeline(
            genes, output_dir=tmp, use_prism=use_prism, min_score=min_score,
            max_drugs_per_gene=max_drugs_per_gene, cancer_type=cancer_type or None,
            skip_plots=True, progress_cb=_cb)
    prog.empty(); status.empty()
    return df


@st.cache_data(show_spinner=False)
def _load_classifier_predictions_cached(path_str: str, mtime: float) -> pd.DataFrame:
    """Cached load of the precomputed protein-classifier CSV. `mtime` is part of
    the cache key so an updated file (e.g. a fresh export dropped in place) is
    picked up without restarting the app."""
    return protein_predictions.load_predictions(Path(path_str))


with tab_drug:
    st.subheader("Drug targets — known drugs, FDA status & cancer-specific scoring")
    st.caption("For a batch of gene(s): queries MyGene.info, ChEMBL, DGIdb, and Open Targets "
               "for known drugs, scores each gene-drug pair by clinical maturity + cancer "
               "relevance, and (optionally) PRISM/DepMap sensitivity. On-demand; not part of "
               "the ranking pipeline. Ported from a collaborator's pipeline "
               "(see `user_data/drug_targets/README.md` for full scoring rationale).")

    _dt_default_n = min(10, len(ranked_visible)) if (ranked_visible is not None and not ranked_visible.empty) else 0
    _dt_default = ", ".join(list(ranked_visible.index[:_dt_default_n])) if _dt_default_n else ""
    _sync_default_genes("dt_genes", _dt_default)
    dt_genes_txt = st.text_input("Gene(s), comma or space separated", key="dt_genes",
                                 help="Auto-populated with the current top-10 ranked genes; "
                                      "updates when parameters change the ranking, unless "
                                      "you've edited it since the last change.")
    dt_genes = [g.strip().upper() for g in dt_genes_txt.replace(",", " ").split() if g.strip()]

    dc1, dc2 = st.columns([1, 1])
    with dc1:
        dt_top_n = st.slider("max genes to score", 1, 50, 15,
                             help="Bounded batch size — each gene queries 3-4 external "
                                  "APIs sequentially. ~1-6 min/gene depending on how well-"
                                  "studied it is (e.g. ALK ~5 min; most genes are faster).")
    with dc2:
        cancer_mode = st.checkbox("Cancer-specific 5-tier scoring", value=True,
                                  help="On: score drugs by relevance to your selected cancer "
                                       "type (5 tiers). Off: cancer-agnostic 3-tier scoring "
                                       "(just 'is this druggable at all').")
    dt_cancer = ""
    if cancer_mode:
        dt_cancer = st.text_input("Cancer type (Open Targets-resolved)", value=disease,
                                  help="Pre-filled from the selected disease; edit if you "
                                       "want a different/broader term (e.g. 'Lung cancer').")

    with st.expander("Thresholds & options", expanded=False):
        dt_min_score = st.slider("min score to include in main results", 1,
                                 5 if cancer_mode else 3, 2)
        dt_max_drugs = st.slider("max drugs per gene shown in heatmap", 5, 50, 15, 5)
        use_prism = st.checkbox("Enable PRISM/DepMap sensitivity enrichment", value=False,
                                help="Adds cell-line drug-sensitivity evidence. Downloads "
                                     "~100 MB on first use and adds ~5 min. Off by default.")

    if dt_genes[:dt_top_n]:
        n_shown = len(dt_genes[:dt_top_n])
        est_lo, est_hi = n_shown * 1, n_shown * 6
        st.caption(f"Will score **{n_shown}** gene(s). Rough first-run estimate: "
                   f"**{est_lo}-{est_hi} min** (cached per gene+cancer after the first run).")

    run_dt = st.button("Run drug-target scoring", type="primary", disabled=not dt_genes)
    dt_table = None
    if run_dt:
        genes_batch = dt_genes[:dt_top_n]
        try:
            dt_table = _run_drug_targets(genes_batch, dt_cancer if cancer_mode else None,
                                         dt_min_score, dt_max_drugs, use_prism)
            st.session_state["drugtargets_last"] = (dt_table, cancer_mode, dt_cancer)
        except Exception as e:  # noqa: BLE001
            st.error(f"Drug-target query failed: {e}. This needs outbound internet "
                     "(MyGene.info / ChEMBL / DGIdb / Open Targets); retry, or reduce the "
                     "gene count.")

    if dt_table is None and "drugtargets_last" in st.session_state:
        dt_table, cancer_mode, dt_cancer = st.session_state["drugtargets_last"]

    if dt_table is not None and not dt_table.empty:
        n_genes = dt_table["Gene"].nunique()
        st.markdown(f"**{len(dt_table)} gene-drug pairs** across **{n_genes} genes**.")
        if cancer_mode:
            for s, label in [(5, "Approved for your cancer"), (4, "Clinical trial for your cancer"),
                             (3, "Approved for other cancers")]:
                n = int((dt_table["Score"] == s).sum())
                if n:
                    st.caption(f"Score {s} ({label}): **{n}**")
        st.plotly_chart(plots.drug_target_heatmap(
            dt_table, dt_max_drugs, vmax=5 if cancer_mode else 3,
            title="Genes × drugs (color = score)"), use_container_width=True)
        if cancer_mode:
            st.plotly_chart(plots.drug_target_score_bar(
                dt_table, f"Score tiers per gene — {dt_cancer}"), use_container_width=True)
        show_cols = [c for c in ["Gene", "Drug_Name", "Max_Phase", "FDA_Approved", "Score",
                                 "Score_Label", "Cancer_Types", "Best_Potency_nM"] if c in dt_table.columns]
        st.dataframe(dt_table[show_cols].sort_values(["Gene", "Score"], ascending=[True, False]),
                    use_container_width=True, height=380)
        st.download_button("Download drug-target results CSV", dt_table.to_csv(index=False).encode(),
                           file_name=f"drug_targets_{disease}.csv".replace(" ", "_"), mime="text/csv")


with tab_clf:
    st.subheader("Precomputed protein-complex-classifier predictions")
    _clf_path = RAW / "top_neuroblastoma_predictions" / "neuroblastoma_extended_list_predictions.csv"
    st.caption(
        "Precomputed P(interactor with the INO80/SRCAP complex) from a protein-language-model "
        "ensemble classifier (see `.pi/skills/protein-complex-classifier`), NOT computed live by "
        "this app — this tab only loads and displays a fixed CSV export "
        f"(`{_clf_path.relative_to(RAW.parent)}`). One of the model's feature blocks is sensitive "
        "to which other proteins were embedded in the same batch, so these numbers are only "
        "comparable to each other (all scored together in one batch) — not to numbers from a "
        "different export or a smaller ad-hoc run. A gene absent from this file shows **Not "
        "generated/cached**, never a fabricated score.")

    _clf_mtime = _clf_path.stat().st_mtime if _clf_path.exists() else 0.0
    _clf_preds = _load_classifier_predictions_cached(str(_clf_path), _clf_mtime)
    if not _clf_path.exists():
        st.warning(f"No predictions file found at `{_clf_path.relative_to(RAW.parent)}`. "
                   "Every requested gene below will show as not generated/cached.")

    _clf_default = ", ".join(list(ranked_visible.index[:10])) if (ranked_visible is not None and not ranked_visible.empty) else ""
    _sync_default_genes("clf_genes", _clf_default)
    clf_genes_txt = st.text_input("Gene(s), comma or space separated", key="clf_genes",
                                  help="Auto-populated with the current top-10 ranked genes; "
                                       "updates when parameters change the ranking, unless "
                                       "you've edited it since the last change.")
    clf_genes = [g.strip().upper() for g in clf_genes_txt.replace(",", " ").split() if g.strip()]

    if not clf_genes:
        st.info("Enter at least one gene.")
    else:
        clf_lookup = protein_predictions.lookup(clf_genes, _clf_preds)
        n_cached = int(clf_lookup["cached"].sum())
        st.markdown(f"**{n_cached}/{len(clf_lookup)} genes** have a cached prediction "
                    f"({len(_clf_preds)} genes total in the loaded file).")
        st.plotly_chart(plots.classifier_bar(clf_lookup, "Precomputed classifier P(interactor)"),
                        use_container_width=True)
        show_clf_cols = ["Prediction", "Probability", "Confidence", "UniProt_ID", "cached"]
        st.dataframe(clf_lookup[show_clf_cols], use_container_width=True, height=360)
        st.download_button("Download predictions CSV", clf_lookup.to_csv().encode(),
                           file_name=f"classifier_predictions_{disease}.csv".replace(" ", "_"),
                           mime="text/csv")


with tab_ins:
    st.subheader("Gene relevance narrative")
    st.caption("Type one or more ranked genes; the agent generates a scientific narrative "
               "on their relevance (co-expression with the seed complex + cohort-selective "
               "essentiality), grounded in the computed metrics.")
    from agents.narrative import (gene_narrative_deterministic, gene_narrative_llm,
                                  list_ollama_models, list_cloud_models,
                                  active_backend, set_api_key, set_base_url, llm_selftest)

    def _fig_to_png_bytes(fig, width: int = 900, height: int = 520) -> bytes | None:
        """Rasterize a Plotly figure to PNG bytes for PDF embedding, via kaleido.
        kaleido v1.x needs a real Chrome/Chromium ALREADY installed on this machine
        (it no longer bundles one) -- returns None on ANY failure (kaleido missing,
        no Chrome, a renderer crash) so PDF export degrades to "see the on-screen
        plots instead" rather than crashing; see app/README.md for the prerequisite.
        """
        try:
            return fig.to_image(format="png", width=width, height=height, scale=2)
        except Exception:
            return None

    def _gather_insight_plots(req: list[str]) -> list[tuple[str, object]]:
        """Assemble whichever plots from the other 8 tabs have data covering `req` --
        reuses each tab's own cached/session-state result; never triggers a new
        network fetch or heavy recomputation from here."""
        out: list[tuple[str, object]] = []
        if ctx.coexpr is not None and not ctx.coexpr.empty:
            out.append(("Co-expression volcano (requested gene(s) circled)",
                       plots.volcano(ctx.coexpr, r_min, q_max,
                                    f"Co-expression volcano — {disease}/{age}",
                                    positive_only, highlight_genes=req)))
        if ranked_visible is not None and not ranked_visible.empty:
            out.append(("Essentiality: % co-expressed vs % essential (requested gene(s) circled)",
                       plots.essentiality_scatter(ranked_visible, "% co-expressed vs % essential",
                                                  highlight_genes=req)))
            out.append(("Top genes by selectivity score (requested gene(s) in crimson)",
                       plots.ranked_bar(ranked_visible, 20, "Top genes by selectivity score",
                                       highlight_genes=req)))
        if fold_string and ranked_visible is not None and not ranked_visible.empty:
            _top = ranked_visible.head(string_topn)
            _rr, _phys_m, _func_m, _edges = string_rerank_cached(
                tuple(_top.index), tuple(_top["selectivity_score"]), tuple(genes),
                string_topn, sw_sel, sw_phys, sw_func, string_rank_cut)
            out.append(("STRING physical association — candidate × seed (from Interactions tab)",
                       plots.string_heatmap(_phys_m, "STRING physical — candidate × seed")))
        if "gtex_last" in st.session_state:
            _tt, _gtex_heat, _tissues, _targets = st.session_state["gtex_last"]
            if _tt is not None and not _tt.empty and any(g in _tt.index for g in req):
                out.append(("Tissue-specificity heatmap (from Tissue Specificity tab's last run)",
                           plots.tissue_heatmap(
                               _gtex_heat, "Expression: normal tissue → TCGA → target cohort")))
        if "drugtargets_last" in st.session_state:
            _dtb, _dt_cancer_mode, _dt_cancer = st.session_state["drugtargets_last"]
            if _dtb is not None and not _dtb.empty and any(g in set(_dtb["Gene"]) for g in req):
                if _dt_cancer_mode:
                    out.append(("Drug-target score tiers per gene (from Drug Targets tab's last run)",
                               plots.drug_target_score_bar(_dtb, f"Score tiers per gene — {_dt_cancer}")))
                else:
                    out.append(("Genes × drugs heatmap (from Drug Targets tab's last run)",
                               plots.drug_target_heatmap(_dtb, 15, vmax=3,
                                                         title="Genes × drugs (color = score)")))
        if _clf_path.exists():
            _lk = protein_predictions.lookup(req, _clf_preds)
            out.append(("Precomputed protein-classifier P(interactor) (Model Predictions tab)",
                       plots.classifier_bar(_lk, "Precomputed classifier P(interactor)")))
        return out

    def _render_insight_plots(req: list[str], key_prefix: str) -> None:
        with st.expander("📊 Plots from other tabs", expanded=False):
            figs = _gather_insight_plots(req)
            if not figs:
                st.caption("No plots available yet — run the other tabs first "
                          "(Co-expression/Essentiality always show once ranking exists).")
            for i, (caption, fig) in enumerate(figs):
                st.caption(caption)
                st.plotly_chart(fig, use_container_width=True, key=f"{key_prefix}_insight_plot_{i}")

    def _insight_plot_pdf_images(req: list[str]) -> tuple[list[tuple[str, bytes]], bool]:
        """Returns (images, any_failed). `any_failed` is True if at least one available
        plot could not be rasterized (almost always: kaleido installed but no Chrome/
        Chromium on this machine) -- the caller adds one clear note to the PDF body
        rather than a per-image failure message."""
        images: list[tuple[str, bytes]] = []
        any_failed = False
        for caption, fig in _gather_insight_plots(req):
            png = _fig_to_png_bytes(fig)
            if png is not None:
                images.append((caption, png))
            else:
                any_failed = True
        return images, any_failed

    ranked = ranked_visible
    if ranked is None or ranked.empty:
        st.info("No ranked genes yet — run an analysis that produces ranked candidates.")
    else:
        # --- LLM backend config (optional cloud key) ---
        with st.expander("LLM backend / API key", expanded=False):
            prov = st.selectbox("Cloud provider (optional)",
                                ["anthropic", "openai", "openrouter"], index=0)
            key_in = st.text_input(f"{prov} API key", type="password",
                                   help="Kept in this process only (not written to disk). "
                                        "Or set it in a .env file / shell env.")
            base_in = st.text_input(
                f"{prov} base URL (optional)", value="",
                help="Leave blank for the provider's default endpoint. Set this to use a "
                     "Claude-compatible endpoint other than api.anthropic.com — e.g. "
                     "Microsoft Azure AI Foundry: "
                     "https://{your-resource}.services.ai.azure.com/anthropic/v1 "
                     "(then the model field below is your Foundry *deployment name*, and "
                     "the API key is your Foundry key — same x-api-key header, no other "
                     "change needed).")
            kc1, kc2 = st.columns([1, 2])
            if kc1.button("Use key"):
                if key_in:
                    set_api_key(prov, key_in)
                    if base_in.strip():
                        set_base_url(prov, base_in)
                    st.success(f"{prov} key set for this session.")
                else:
                    st.warning("Enter a key first.")
            if kc2.button("Test LLM connection"):
                with st.spinner("Testing…"):
                    ok, msg = llm_selftest()
                st.success(f"✓ {msg}") if ok else st.error(f"✗ {msg}")
            b = active_backend()
            st.caption(f"Active backend: **{b['provider']}** · `{b['model']}` · `{b['base']}`. "
                       "Priority OpenRouter → OpenAI → Anthropic → local Ollama.")

        backend = active_backend()
        default_top = ", ".join(list(ranked.index[:3]))
        _sync_default_genes("insight_genes", default_top)
        gene_text = st.text_input("Gene(s), comma or space separated", key="insight_genes",
                                  help="Auto-populated with the current top-3 ranked genes; "
                                       "updates when parameters change the ranking, unless "
                                       "you've edited it since the last change.")
        req_genes = [g.strip().upper() for g in gene_text.replace(",", " ").split() if g.strip()]
        oc1, oc2 = st.columns([1, 1])
        fetch_kb = oc1.checkbox("Fetch literature + protein/drug annotations", value=True,
                                help="Europe PMC papers + UniProt protein annotation + "
                                     "ChEMBL known drugs for each gene; also fed to the LLM.")
        n_papers = oc2.slider("papers per gene", 3, 10, 6)
        os1, os2 = st.columns([1, 1])
        fetch_string = os1.checkbox("Fetch STRING association with seeds", value=True,
                                    help="Functional + physical STRING scores vs the seed "
                                         "complex; shown and fed to the LLM. Uses the SAME "
                                         "edge cutoff as the sidebar's Interactions settings "
                                         "(the sidebar's selectivity/physical/functional "
                                         "weights don't apply here: that combined string_score "
                                         "is normalized within the Interactions tab's top-N "
                                         "batch, so it isn't a portable per-gene number — "
                                         "this shows the raw functional/physical evidence "
                                         "instead).")
        string_cut = string_rank_cut
        os2.metric("STRING edge cutoff", f"{string_cut:g}",
                   help="From the sidebar's STRING interactions section — change it there "
                        "to affect this too." + ("" if fold_string else 
                        " (sidebar default — 'Fold STRING into ranking' is currently off)"))
        use_tissue = st.checkbox(
            "Include last Tissue Specificity triage in the comparative table (if available)",
            value=True,
            help="Reuses the most recent result from the Tissue Specificity tab for any "
                 "typed genes included in that run — does not trigger a new Xena query here. "
                 "Feeds only the deterministic comparative table's 'Tissue gates' column, not "
                 "the LLM's own prose — the pass/fail gates mix a lenient median-based check "
                 "with a stricter max-based ranking score, a nuance too easy for a model to "
                 "flatten into an oversimplified sentence, so it narrates only the plain "
                 "code-computed pass/fail via the table instead.")
        use_drug = st.checkbox(
            "Use last Drug Targets scoring as evidence (if available)", value=True,
            help="Reuses the most recent result from the Drug Targets tab for any typed "
                 "genes included in that run — does not trigger a new query here.")
        use_biomcp = st.checkbox(
            "Fetch BioMCP pathways + tissue atlas (if installed)", value=False,
            help="Reactome/KEGG pathway membership + Human Protein Atlas (HPA) tissue "
                 "expression and subcellular localization, via the optional `biomcp` CLI "
                 "(biomcp.org, MIT-licensed) — install separately with "
                 "`uv tool install biomcp-cli`; not bundled with this app. Shown in-tab "
                 "and fed to the LLM. Off by default since it needs an external binary; "
                 "shows a clear 'not installed' message rather than silently omitting "
                 "the section if you enable it without installing biomcp.")
        use_plots = st.checkbox(
            "Include plots from other tabs (Co-expression/Essentiality/Interactions/"
            "Tissue Specificity/Drug Targets/Model Predictions) in this report", value=False,
            help="Reuses each tab's own cached last-run result (no new computation "
                 "triggered here); Co-expression/Essentiality/Model Predictions are "
                 "always available once ranking exists, the other three only if you've "
                 "run that tab. Shown on-screen always. Embedding them in the PDF export "
                 "additionally needs the `kaleido` package AND a real Chrome/Chromium "
                 "already installed on this machine (kaleido v1.x no longer bundles one) "
                 "— if that prerequisite is missing, the PDF notes it and you still get "
                 "the on-screen plots.")
        cola, colb = st.columns([1, 1])
        gen_det = cola.button("Generate narrative", type="primary", use_container_width=True)
        model = None
        gen_llm = False
        if backend["provider"] == "ollama":
            models = list_ollama_models()
            if models:
                default_ix = models.index("llama3.1:8b") if "llama3.1:8b" in models else 0
                model = colb.selectbox("Ollama model (optional LLM)", models, index=default_ix,
                                       key="insight_model")
                gen_llm = colb.button("Generate with LLM", use_container_width=True)
            else:
                colb.caption("No local Ollama models detected for the optional LLM narrative.")
        else:
            cloud_models = list_cloud_models()
            if cloud_models:
                dflt = backend["model"] if backend["model"] in cloud_models else cloud_models[0]
                model = colb.selectbox(f"{backend['provider']} model", cloud_models,
                                       index=cloud_models.index(dflt), key="insight_cloud_model")
            else:
                colb.caption(f"LLM: **{backend['provider']}** · `{backend['model']}` "
                             "(could not list models; using default)")
            gen_llm = colb.button("Generate with LLM", use_container_width=True)

        if req_genes:
            missing = [g for g in req_genes if g not in ranked.index]
            if missing:
                st.caption(f"Not in ranked set: {', '.join(missing)}")

        if gen_det and req_genes:
            st.session_state["det_narrative_result"] = {
                "genes": tuple(req_genes),
                "text": gene_narrative_deterministic(ctx, req_genes),
            }

        det_result = st.session_state.get("det_narrative_result")
        if det_result and det_result["genes"] == tuple(req_genes):
            ui.markdown_with_highlights(det_result["text"])
            det_images, det_imgs_failed = ([], False)
            if use_plots:
                det_images, det_imgs_failed = _insight_plot_pdf_images(req_genes)
            det_body = det_result["text"]
            if use_plots and det_imgs_failed:
                det_body += ("\n\n*(Some plots could not be embedded in this PDF \u2014 "
                            "kaleido needs a real Chrome/Chromium installed on this "
                            "machine; see the on-screen plots below/in the app.)*")
            det_pdf = pdf_export.markdown_to_pdf_bytes(
                "Gene relevance narrative \u2014 deterministic",
                f"{ctx.spec.age_group} / {ctx.spec.disease} \u00b7 seeds={ctx.seeds.source} \u00b7 "
                f"genes: {', '.join(det_result['genes'])}",
                det_body, images=det_images or None,
            )
            st.download_button(
                "\U0001f4c4 Export as PDF", data=det_pdf, key="det_pdf_dl",
                file_name=f"scrap-ai_deterministic_{'-'.join(det_result['genes'])}.pdf",
                mime="application/pdf",
            )
            if fetch_string:
                _render_string(req_genes, genes, string_cut)
            if fetch_kb:
                _render_dossiers(req_genes, disease, n_papers)
            if use_biomcp:
                _render_biomcp(req_genes)
            if use_plots:
                _render_insight_plots(req_genes, "det")
            if use_tissue and "gtex_last" in st.session_state:
                _tt, _, _, _ = st.session_state["gtex_last"]
                hit = [g for g in req_genes if g in _tt.index]
                if hit:
                    st.caption(f"Tissue-specificity evidence available for: {', '.join(hit)} "
                               "(from the Tissue Specificity tab's last run).")
            if use_drug and "drugtargets_last" in st.session_state:
                _dtb, _, _ = st.session_state["drugtargets_last"]
                hit = [g for g in req_genes if g in set(_dtb["Gene"])]
                if hit:
                    st.caption(f"Drug-target evidence available for: {', '.join(hit)} "
                               "(from the Drug Targets tab's last run).")
        if gen_llm and req_genes:
            dossiers = None
            string_data = None
            tissue_data = None
            drug_data = None
            biomcp_data = None
            if fetch_kb:
                with st.spinner("Fetching literature + annotations…"):
                    dossiers = _fetch_dossiers(tuple(req_genes), disease, n_papers)
            if fetch_string:
                with st.spinner("Fetching STRING scores…"):
                    string_data = _fetch_string(tuple(req_genes), tuple(genes), string_cut)
            if use_biomcp:
                with st.spinner("Fetching BioMCP pathways + tissue atlas…"):
                    biomcp_data = _fetch_biomcp(tuple(req_genes))
            if use_tissue and "gtex_last" in st.session_state:
                _tt, _, _tissues, _targets = st.session_state["gtex_last"]
                tissue_data = {}
                for g in req_genes:
                    if g in _tt.index:
                        row = _tt.loc[g]
                        tissue_data[g] = {
                            "gtex_tissues_checked": _tissues, "target_cohort": _targets,
                            "normal_max_tpm": round(float(row["normal_max"]), 2),
                            "pass_normal_silent": bool(row["pass_normal"]),
                            "pancancer_low_fraction": round(float(row["pancancer_low_fraction"]), 3),
                            "pass_pancancer": bool(row["pass_pancancer"]),
                            "target_median_tpm": round(float(row["target_median_tpm"]), 2),
                            "pass_target": bool(row["pass_target"]),
                            "pass_all_gates": bool(row["pass_all"]),
                            "n_gates_passed": int(row["n_gates_passed"]),
                            "specificity_tier": str(row["specificity_tier"]),
                        }
            if use_drug and "drugtargets_last" in st.session_state:
                _dtb, _dt_cancer_mode, _dt_cancer_name = st.session_state["drugtargets_last"]
                drug_data = {}
                for g in req_genes:
                    grows = _dtb[_dtb["Gene"] == g].sort_values("Score", ascending=False)
                    if grows.empty:
                        continue
                    top = grows.head(5)
                    cols = [c for c in ["Drug_Name", "Max_Phase", "FDA_Approved", "Score",
                                        "Score_Label", "Cancer_Types"] if c in top.columns]
                    drug_data[g] = {
                        "cancer_specific_mode": bool(_dt_cancer_mode),
                        "query_cancer": _dt_cancer_name if _dt_cancer_mode else None,
                        "n_drugs_found": int(len(grows)),
                        "top_drugs": top[cols].to_dict("records"),
                    }
            label = model or backend["model"]
            with st.spinner(f"Generating with {backend['provider']} / {label}…"):
                try:
                    txt = gene_narrative_llm(ctx, req_genes, model=model, timeout=600,
                                             dossiers=dossiers, string_data=string_data,
                                             tissue_data=tissue_data, drug_data=drug_data,
                                             biomcp_data=biomcp_data)
                    if not txt.strip():
                        st.session_state.pop("llm_narrative_result", None)
                        st.error("LLM returned an empty response — it likely hit its output-"
                                 "token budget (e.g. an extended-thinking model spending the "
                                 "whole budget on internal reasoning) rather than failing "
                                 "outright. Try fewer genes at once, or use **Generate "
                                 "narrative** for the deterministic version.")
                    else:
                        st.session_state["llm_narrative_result"] = {
                            "genes": tuple(req_genes), "text": txt,
                            "string_data": string_data, "dossiers": dossiers,
                            "biomcp_data": biomcp_data,
                        }
                except (TimeoutError, OSError) as e:
                    st.session_state.pop("llm_narrative_result", None)
                    st.error(f"LLM timed out/failed ({e}). Use **Generate narrative** for the "
                             "deterministic version, or retry.")
                except ValueError as e:
                    st.session_state.pop("llm_narrative_result", None)
                    st.error(str(e))
                except Exception as e:  # noqa: BLE001
                    st.session_state.pop("llm_narrative_result", None)
                    st.error(f"LLM call failed: {e}")

        llm_result = st.session_state.get("llm_narrative_result")
        if llm_result and llm_result["genes"] == tuple(req_genes):
            with ui.card_container("llm-narrative-card", "LLM narrative", icon="\U0001f916"):
                st.warning("LLM output — grounded in the computed metrics + retrieved "
                           "sources; verify every claim and citation before use.")
                ui.markdown_with_highlights(llm_result["text"])
            llm_images, llm_imgs_failed = ([], False)
            if use_plots:
                llm_images, llm_imgs_failed = _insight_plot_pdf_images(req_genes)
            llm_body = llm_result["text"]
            if use_plots and llm_imgs_failed:
                llm_body += ("\n\n*(Some plots could not be embedded in this PDF \u2014 "
                            "kaleido needs a real Chrome/Chromium installed on this "
                            "machine; see the on-screen plots below/in the app.)*")
            llm_pdf = pdf_export.markdown_to_pdf_bytes(
                "Gene relevance narrative \u2014 LLM narrative",
                f"{ctx.spec.age_group} / {ctx.spec.disease} \u00b7 seeds={ctx.seeds.source} \u00b7 "
                f"genes: {', '.join(llm_result['genes'])} \u00b7 {backend['provider']}/{backend['model']}",
                llm_body, images=llm_images or None,
            )
            st.download_button(
                "\U0001f4c4 Export as PDF", data=llm_pdf, key="llm_pdf_dl",
                file_name=f"scrap-ai_llm-narrative_{'-'.join(llm_result['genes'])}.pdf",
                mime="application/pdf",
            )
            if fetch_string and llm_result.get("string_data"):
                _render_string(req_genes, genes, string_cut, prefetched=llm_result["string_data"])
            if fetch_kb and llm_result.get("dossiers"):
                _render_dossiers(req_genes, disease, n_papers, prefetched=llm_result["dossiers"])
            if use_biomcp and llm_result.get("biomcp_data") is not None:
                _render_biomcp(req_genes, prefetched=llm_result["biomcp_data"])
            if use_plots:
                _render_insight_plots(req_genes, "llm")
