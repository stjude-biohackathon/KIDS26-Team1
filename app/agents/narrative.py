"""Narrative synthesis: a deterministic report from the computed numbers, plus an
optional local-Ollama rephrasing that is strictly grounded in those numbers.

The deterministic text is the source of truth. The LLM step only reorganizes the
facts we hand it; the prompt forbids introducing new genes or numbers.
"""
from __future__ import annotations

import json
import os
import urllib.request
import urllib.error

import pandas as pd

from agents.base import AnalysisContext

OLLAMA_BASE = os.environ.get("OLLAMA_BASE_URL", "http://localhost:11434").rstrip("/")
DEFAULT_OLLAMA_MODEL = "llama3.1:8b"


def load_dotenv(path: str | None = None) -> list[str]:
    """Minimal .env loader (no dependency). Sets keys not already in the env.

    Looks at `path`, else the sandbox root and app/ dirs. Returns the names of
    the vars it set. Existing environment values always win.
    """
    from pathlib import Path as _P
    candidates = [path] if path else []
    here = _P(__file__).resolve()
    candidates += [here.parent.parent.parent / ".env", here.parent.parent / ".env"]
    set_keys: list[str] = []
    for c in candidates:
        if not c:
            continue
        p = _P(c)
        if not p.is_file():
            continue
        for line in p.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            k = k.strip()
            v = v.strip().strip('"').strip("'")
            if k and k not in os.environ:
                os.environ[k] = v
                set_keys.append(k)
    return set_keys


# auto-load a .env on import so the app/tests pick up keys without shell exports
load_dotenv()

# Cloud providers are used automatically when a key is present, in this priority.
# Each entry: (env var, default model, base url, style). style: "openai" | "anthropic".
_CLOUD_PROVIDERS = [
    ("OPENROUTER_API_KEY", "openai/gpt-4o-mini",
     "https://openrouter.ai/api/v1", "openai"),
    ("OPENAI_API_KEY", "gpt-4o-mini",
     "https://api.openai.com/v1", "openai"),
    ("ANTHROPIC_API_KEY", "claude-haiku-4-5-20251001",
     "https://api.anthropic.com/v1", "anthropic"),
]


def active_backend() -> dict:
    """Which LLM backend is used by default.

    Returns {provider, model, style, base, env}. A present cloud key wins (in
    priority order); otherwise local Ollama with DEFAULT_OLLAMA_MODEL.

    Each cloud provider's base URL can be overridden with a `{PROVIDER}_BASE_URL`
    env var (e.g. `ANTHROPIC_BASE_URL`) — this is how you point the Anthropic
    entry at a Claude-compatible endpoint that isn't api.anthropic.com, such as
    Microsoft Azure AI Foundry (`https://{resource}.services.ai.azure.com/anthropic/v1`).
    Azure Foundry uses the same Messages API wire format and accepts the same
    `x-api-key` header, so no other code changes are needed — just set
    ANTHROPIC_API_KEY to your Azure key, ANTHROPIC_BASE_URL to the Foundry URL
    above, and ANTHROPIC_MODEL to your Foundry *deployment name* (not
    necessarily the public Claude model id).
    """
    for env, model, base, style in _CLOUD_PROVIDERS:
        if os.environ.get(env):
            provider = env.replace("_API_KEY", "").lower()
            model = os.environ.get(f"{provider.upper()}_MODEL", model)
            base = os.environ.get(f"{provider.upper()}_BASE_URL", base)
            return {"provider": provider, "model": model, "style": style,
                    "base": base, "env": env}
    return {"provider": "ollama", "model": DEFAULT_OLLAMA_MODEL, "style": "ollama",
            "base": OLLAMA_BASE, "env": None}


def set_api_key(provider: str, key: str) -> None:
    """Set a provider key into the process env at runtime (from the UI field).

    provider: 'anthropic' | 'openai' | 'openrouter'.
    """
    envmap = {"anthropic": "ANTHROPIC_API_KEY", "openai": "OPENAI_API_KEY",
              "openrouter": "OPENROUTER_API_KEY"}
    env = envmap.get(provider)
    if not env:
        raise ValueError(f"unknown provider {provider!r}")
    if key:
        os.environ[env] = key.strip()


def set_base_url(provider: str, url: str) -> None:
    """Set a provider's base-URL override into the process env at runtime.

    provider: 'anthropic' | 'openai' | 'openrouter'. Empty `url` clears the
    override (reverting to that provider's default endpoint). This is how the
    Anthropic entry gets pointed at a Claude-compatible endpoint that isn't
    api.anthropic.com — e.g. Microsoft Azure AI Foundry
    (https://{resource}.services.ai.azure.com/anthropic/v1).
    """
    envmap = {"anthropic": "ANTHROPIC_BASE_URL", "openai": "OPENAI_BASE_URL",
              "openrouter": "OPENROUTER_BASE_URL"}
    env = envmap.get(provider)
    if not env:
        raise ValueError(f"unknown provider {provider!r}")
    url = (url or "").strip()
    if url:
        os.environ[env] = url
    else:
        os.environ.pop(env, None)


def llm_selftest(timeout: float = 60.0) -> tuple[bool, str]:
    """Tiny round-trip to the active backend. Returns (ok, message)."""
    b = active_backend()
    try:
        txt = llm_complete("You are a terse assistant.", "Reply with exactly: OK",
                           timeout=timeout, max_tokens=8)
        ok = "ok" in txt.lower()
        return ok, f"{b['provider']} / {b['model']}: {txt.strip()[:60] or '(empty)'}"
    except Exception as e:  # noqa: BLE001
        return False, f"{b['provider']} / {b['model']} failed: {e}"


def _post_json(url: str, payload: dict, headers: dict, timeout: float) -> dict:
    req = urllib.request.Request(url, data=json.dumps(payload).encode(), headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode())
    except urllib.error.HTTPError as e:
        # Surface the provider's own error body (explains bad model vs bad path).
        try:
            detail = e.read().decode()[:400]
        except Exception:
            detail = ""
        raise RuntimeError(f"HTTP {e.code} from {url} — {detail}") from None


# Params some newer models reject; dropped and retried when the API complains.
_OPTIONAL_PARAMS = ("temperature", "top_p", "top_k")


def _post_retry(url: str, payload: dict, headers: dict, timeout: float,
                _tries: int = 4) -> dict:
    """POST, and if the API rejects an optional sampling param as deprecated/
    unsupported, drop that param and retry (some Claude/OpenAI models forbid it)."""
    p = dict(payload)
    for _ in range(_tries):
        try:
            return _post_json(url, p, headers, timeout)
        except RuntimeError as e:
            msg = str(e).lower()
            if "http 400" not in msg:
                raise
            dropped = None
            for param in _OPTIONAL_PARAMS:
                if param in p and (param in msg):
                    del p[param]
                    dropped = param
                    break
            if dropped is None:
                raise
    return _post_json(url, p, headers, timeout)


def _strip_think(text: str) -> str:
    text = (text or "").strip()
    if "</think>" in text:
        text = text.split("</think>", 1)[1].strip()
    return text


def _join_endpoint(base: str, suffix: str) -> str:
    """Join a configured base URL with a path suffix (e.g. 'messages', 'models'),
    tolerating a base that already ends with that exact suffix.

    Users copying an Azure AI Foundry / APIM "Target URI" often paste the full
    endpoint (already ending in `/messages`) rather than the bare base the
    OpenAI/Anthropic SDKs expect — this avoids silently building a broken
    `.../messages/messages` URL in that case.
    """
    b = base.rstrip("/")
    if b.endswith("/" + suffix):
        return b
    return f"{b}/{suffix}"


def llm_complete(system: str, prompt: str, model: str | None = None,
                 timeout: float = 600.0, temperature: float = 0.2,
                 max_tokens: int = 500) -> str:
    """Unified completion. Cloud (OpenRouter/OpenAI/Anthropic) if a key is set,
    else local Ollama. `model` overrides the backend default."""
    b = active_backend()
    mdl = model or b["model"]
    style = b["style"]

    if style == "openai":
        key = os.environ[b["env"]]
        headers = {"Content-Type": "application/json", "Authorization": f"Bearer {key}"}
        if b["provider"] == "openrouter":
            headers["HTTP-Referer"] = "http://localhost"
            headers["X-Title"] = "INO80 DepMap Explorer"
        payload = {
            "model": mdl, "temperature": temperature, "max_tokens": max_tokens,
            "messages": [{"role": "system", "content": system},
                         {"role": "user", "content": prompt}],
        }
        out = _post_retry(_join_endpoint(b["base"], "chat/completions"), payload, headers, timeout)
        return _strip_think(out["choices"][0]["message"]["content"])

    if style == "anthropic":
        key = os.environ[b["env"]]
        headers = {"Content-Type": "application/json", "x-api-key": key,
                   "anthropic-version": "2023-06-01"}
        payload = {
            "model": mdl, "max_tokens": max_tokens, "temperature": temperature,
            "system": system,
            "messages": [{"role": "user", "content": prompt}],
        }
        out = _post_retry(_join_endpoint(b["base"], "messages"), payload, headers, timeout)
        parts = out.get("content", [])
        return _strip_think("".join(p.get("text", "") for p in parts))

    # ollama
    out = _post_json(f"{OLLAMA_BASE}/api/generate", {
        "model": mdl, "prompt": f"{system}\n\n{prompt}", "stream": False,
        "keep_alive": "10m", "options": {"temperature": temperature, "num_predict": max_tokens},
    }, {"Content-Type": "application/json"}, timeout)
    return _strip_think(out.get("response", ""))


# ---------- deterministic ----------
def deterministic_narrative(ctx: AnalysisContext) -> str:
    p = ctx.params
    pop = "all analyzable lines (pan-cancer)" if p.population == "pan_cancer" else \
        f"{ctx.spec.disease} / {ctx.spec.age_group} lines"
    lines: list[str] = []
    lines.append(
        f"Co-expression with the {ctx.seeds.source} seed proteins: **{p.method}** "
        f"correlation across **{ctx.corr_n}** {pop}. Essentiality metrics are then "
        f"measured across **{ctx.ess_denom}** {ctx.spec.disease}/{ctx.spec.age_group} "
        "cell lines that have both expression and CRISPR data."
    )
    lines.append(
        f"Definitions per cell line: *co-expressed* = candidate expressed (TPM > {p.expr_tpm:g}) "
        f"AND ≥{p.seed_on_frac:.0%} of seed genes expressed; *essential* = Chronos ≤ {p.ess_thresh:g}."
    )

    ranked = ctx.ranked
    if ranked is None or ranked.empty:
        lines.append(
            "\n**No co-expressed candidates** pass the correlation thresholds "
            f"(|r| ≥ {p.r_min}, q ≤ {p.q_max}) with usable essentiality data. "
            "Consider lowering min |r|, raising q, or relaxing the TPM threshold."
        )
        return _with_caveats(ctx, lines)

    n = len(ranked)
    lines.append(
        f"\n**{n} genes** are co-expressed with the seed set. They are ranked by a "
        "**selectivity-weighted score** = joint % × (1 + [cohort−pan-cancer essential%]/100) "
        f"× common-essential penalty ({p.ce_penalty:g}×). This favors genes that are more "
        "essential in this cohort than across all cancers, and down-weights pan-essential genes."
    )
    top = ranked.head(8)
    rows = []
    for g, r in top.iterrows():
        tag = " *(pan-essential)*" if r.get("common_essential", False) else ""
        rows.append(
            f"- **{g}**{tag}: score {r['selectivity_score']:.0f} — joint {r['pct_joint']:.0f}%, "
            f"selectivity Δ {r['selectivity_delta']:+.0f} pts "
            f"(essential {r['pct_essential']:.0f}% here vs {r['pct_essential_pan']:.0f}% pan-cancer)"
        )
    lines.append("\nTop candidates (selectivity-weighted):\n" + "\n".join(rows))

    if "common_essential" in ranked.columns:
        sel = ranked[~ranked["common_essential"]]
        ce = ranked[ranked["common_essential"]]
        if len(ce):
            lines.append(
                f"\n{len(ce)} ranked genes are DepMap **common-essential** "
                f"({', '.join(ce.index[:6])}{'…' if len(ce) > 6 else ''}) — essential across "
                "most lines, so a high joint % partly reflects general essentiality rather "
                "than cohort-selective dependency."
            )
        if len(sel):
            top_sel = sel.head(6)
            lines.append(
                f"\n{len(sel)} are **not** pan-essential — the more selective candidates. "
                f"Top by joint %: {', '.join(f'{g} ({r.pct_joint:.0f}%)' for g, r in top_sel.iterrows())}."
            )
    return _with_caveats(ctx, lines)


def _with_caveats(ctx: AnalysisContext, lines: list[str]) -> str:
    p = ctx.params
    caveats = []
    if ctx.ess_denom and ctx.ess_denom < p.min_cohort_n:
        caveats.append(
            f"the cohort has only {ctx.ess_denom} cell lines with both data types "
            f"(< {p.min_cohort_n}); percentages are coarse"
        )
    if p.population == "pan_cancer":
        caveats.append(
            "candidate co-expression was selected pan-cancer, but the % metrics are still "
            "computed only within the selected cohort"
        )
    if ctx.seeds_missing:
        caveats.append(f"seeds absent from the data: {', '.join(ctx.seeds_missing)}")
    if caveats:
        lines.append("\n**Caveats:** " + "; ".join(caveats) + ".")
    return "\n".join(lines)


# ---------- optional local-LLM enhance ----------
def _facts_payload(ctx: AnalysisContext) -> dict:
    ranked = ctx.ranked
    top = [] if ranked is None or ranked.empty else [
        {
            "gene": g,
            "selectivity_score": round(float(r["selectivity_score"]), 1),
            "pct_joint": round(float(r["pct_joint"]), 1),
            "selectivity_delta": round(float(r["selectivity_delta"]), 1),
            "pct_essential_cohort": round(float(r["pct_essential"]), 1),
            "pct_essential_pan": round(float(r["pct_essential_pan"]), 1),
            "common_essential": bool(r.get("common_essential", False)),
        }
        for g, r in ranked.head(10).iterrows()
    ]
    return {
        "age_group": ctx.spec.age_group,
        "disease": ctx.spec.disease,
        "seed_source": ctx.seeds.source,
        "method": ctx.params.method,
        "population": ctx.params.population,
        "expr_tpm_threshold": ctx.params.expr_tpm,
        "essentiality_threshold": ctx.params.ess_thresh,
        "common_essential_penalty": ctx.params.ce_penalty,
        "n_cohort_lines_with_both_data": ctx.ess_denom,
        "n_ranked": 0 if ranked is None else len(ranked),
        "ranking_metric": "selectivity_score = pct_joint x (1 + (cohort-pan essential%)/100) x common_essential_penalty",
        "top_ranked": top,
    }


def ollama_enhance(ctx: AnalysisContext, model: str | None = None,
                   timeout: float = 600.0) -> str:
    """Rephrase the deterministic facts into prose via the active LLM backend
    (cloud if a key is set, else local Ollama). Grounded in the facts JSON."""
    facts = _facts_payload(ctx)
    system = (
        "You are a cautious cancer-genomics research assistant. You will be given "
        "a JSON object of ALREADY-COMPUTED results from a DepMap analysis: genes "
        "co-expressed with a seed protein complex, ranked by the % of the cohort's "
        "cell lines where each gene is BOTH co-expressed and essential. Write a "
        "concise (120-180 word) interpretation for a researcher. STRICT RULES: use "
        "ONLY the genes and numbers in the JSON; do NOT invent gene names, numbers, "
        "pathways, or citations; if you mention a gene it must appear in top_ranked; "
        "treat common_essential=true genes as poor selective targets. End with one "
        "sentence on an appropriate next validation step. Do not claim clinical relevance."
    )
    prompt = f"FACTS (JSON):\n{json.dumps(facts, indent=2)}\n\nWrite the interpretation now."
    return llm_complete(system, prompt, model=model, timeout=timeout, max_tokens=320)


def list_ollama_models() -> list[str]:
    try:
        with urllib.request.urlopen(f"{OLLAMA_BASE}/api/tags", timeout=4) as resp:
            data = json.loads(resp.read().decode())
        return [m["name"] for m in data.get("models", [])]
    except Exception:
        return []


def list_cloud_models(timeout: float = 15.0) -> list[str]:
    """List available model ids for the active cloud backend (empty on failure)."""
    b = active_backend()
    if b["provider"] == "ollama" or not b["env"]:
        return []
    key = os.environ.get(b["env"], "")
    try:
        if b["style"] == "anthropic":
            req = urllib.request.Request(
                _join_endpoint(b["base"], "models"),
                headers={"x-api-key": key, "anthropic-version": "2023-06-01"})
        else:  # openai / openrouter
            req = urllib.request.Request(
                _join_endpoint(b["base"], "models"), headers={"Authorization": f"Bearer {key}"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode())
        return [m["id"] for m in data.get("data", [])]
    except Exception:
        return []


# ---------- per-gene narrative (Insights tab) ----------
def gene_facts(ctx: AnalysisContext, gene: str) -> dict | None:
    """Structured facts for one gene from the ranked table, or None if absent."""
    ranked = ctx.ranked
    if ranked is None or ranked.empty or gene not in ranked.index:
        return None
    r = ranked.loc[gene]
    rank_pos = int(ranked.index.get_loc(gene)) + 1
    return {
        "gene": gene,
        "rank": rank_pos,
        "of_total": int(len(ranked)),
        "selectivity_score": round(float(r["selectivity_score"]), 1),
        "pct_joint": round(float(r["pct_joint"]), 1),
        "selectivity_delta": round(float(r["selectivity_delta"]), 1),
        "pct_essential_cohort": round(float(r["pct_essential"]), 1),
        "pct_essential_pan": round(float(r["pct_essential_pan"]), 1),
        "pct_coexpressed": round(float(r["pct_coexpressed"]), 1),
        "pct_essential_given_coexpr": round(float(r["pct_essential_given_coexpr"]), 1),
        "best_seed": str(r.get("best_seed", "")),
        "mean_abs_r": round(float(r.get("mean_abs_r", float("nan"))), 3),
        "common_essential": bool(r.get("common_essential", False)),
    }


def gene_narrative_deterministic(ctx: AnalysisContext, genes: list[str]) -> str:
    """Deterministic per-gene relevance summary grounded in the ranked table."""
    out: list[str] = []
    n_found = 0
    for g in genes:
        f = gene_facts(ctx, g)
        if f is None:
            total = 0 if ctx.ranked is None else len(ctx.ranked)
            out.append(f"- **{g}**: not among the {total} ranked co-expressed genes "
                       "for this cohort/threshold set.")
            continue
        n_found += 1
        sel = "pan-essential (poor selective target)" if f["common_essential"] else "selective"
        direction = ("more essential in this cohort than pan-cancer"
                     if f["selectivity_delta"] > 0 else
                     "not preferentially essential here vs pan-cancer")
        out.append(
            f"- **{g}** — rank **#{f['rank']}/{f['of_total']}** by selectivity score "
            f"({f['selectivity_score']:.0f}). Co-expressed with the seed set in "
            f"{f['pct_coexpressed']:.0f}% of cohort lines (strongest seed {f['best_seed']}, "
            f"mean|r|={f['mean_abs_r']:.2f}); essential in {f['pct_essential_cohort']:.0f}% of "
            f"cohort lines vs {f['pct_essential_pan']:.0f}% pan-cancer (delta {f['selectivity_delta']:+.0f} pts, "
            f"{direction}). Joint co-expressed+essential in {f['pct_joint']:.0f}% of lines. "
            f"Classification: **{sel}**."
        )
    header = (f"### Relevance of {n_found} gene(s) in {ctx.spec.age_group} / "
              f"{ctx.spec.disease} (seeds={ctx.seeds.source})\n")
    return header + "\n".join(out)


def comparative_facts_table(ctx: AnalysisContext, genes: list[str],
                            dossiers: dict | None = None, string_data: dict | None = None,
                            tissue_data: dict | None = None, drug_data: dict | None = None,
                            ) -> tuple[list[str], list[list[str]]]:
    """Deterministic comparison table across the requested genes' already-computed
    metrics. Built directly from the same facts fed to the LLM (gene_facts() plus
    whatever RAG blocks were fetched) -- never authored by the model, so the numbers
    in it are exactly as trustworthy as the rest of the app's ranking. Returns
    (headers, rows) as plain strings, ready to render as markdown (on-screen) or a
    native table (PDF export) without re-parsing.
    """
    headers = ["Gene", "Rank", "Selectivity score", "% essential (cohort)",
              "% essential (pan-cancer)", "Selectivity delta", "Common-essential?",
              "STRING physical", "STRING functional", "Known drugs", "Tissue gates"]
    rows: list[list[str]] = []
    for g in genes:
        f = gene_facts(ctx, g)
        if f is None:
            rows.append([g, "not ranked", "\u2013", "\u2013", "\u2013", "\u2013", "\u2013",
                        "\u2013", "\u2013", "\u2013", "\u2013"])
            continue
        sd = (string_data or {}).get(g)
        phys = f"{sd['physical']['max']:.2f}" if sd else "n/a"
        func = f"{sd['functional']['max']:.2f}" if sd else "n/a"
        # Prefer the Drug Targets tab's richer DGIdb/ChEMBL/Open Targets scoring
        # (drug_data, keyed by n_drugs_found) when it was fetched -- it is the same
        # data the LLM's prose is grounded in for this gene. Fall back to the
        # simpler per-gene ChEMBL-only dossier count only when drug_data wasn't
        # fetched at all; a bug previously ignored drug_data entirely here, so a
        # gene with real Drug Targets hits but 0 ChEMBL mechanisms showed "0".
        dt = (drug_data or {}).get(g) if drug_data else None
        if dt is not None:
            n_drugs = str(dt.get("n_drugs_found", 0))
        elif dossiers:
            dd = (dossiers or {}).get(g) or {}
            ch = dd.get("chembl") or {}
            n_drugs = str(ch.get("n_mechanisms", 0))
        else:
            n_drugs = "n/a"
        td = (tissue_data or {}).get(g) if tissue_data else None
        if td is None:
            gates = "n/a"
        elif "specificity_tier" in td:
            # graded fail/low/medium/high (see tissue_specificity.py) -- distinguishes a
            # gene passing 2/3 gates from one passing 0/3, which the older pass_all_gates
            # boolean alone collapsed into the same "fail".
            gates = str(td["specificity_tier"])
        else:
            # backward compatible with callers still passing only the older
            # {"pass_all_gates": bool} shape.
            gates = "pass" if td.get("pass_all_gates") else "fail"
        rows.append([
            g, f"{f['rank']}/{f['of_total']}", f"{f['selectivity_score']:.1f}",
            f"{f['pct_essential_cohort']:.1f}%", f"{f['pct_essential_pan']:.1f}%",
            f"{f['selectivity_delta']:+.1f}", "yes" if f["common_essential"] else "no",
            phys, func, n_drugs, gates,
        ])
    return headers, rows


def comparative_table_markdown(headers: list[str], rows: list[list[str]]) -> str:
    """Render (headers, rows) as a GitHub-flavored markdown pipe table, with a
    heading + caption that makes clear it is computed, not model-authored."""
    def esc(s: str) -> str:
        return str(s).replace("|", "\\|").replace("\n", " ")
    lines = [
        "## Comparative Summary Table",
        "*(deterministic \u2014 computed directly from the same facts above, not authored "
        "by the model)*",
        "",
        "| " + " | ".join(esc(h) for h in headers) + " |",
        "|" + "|".join(["---"] * len(headers)) + "|",
    ]
    for row in rows:
        lines.append("| " + " | ".join(esc(c) for c in row) + " |")
    return "\n".join(lines)


def gene_narrative_llm(ctx: AnalysisContext, genes: list[str],
                       model: str | None = None, timeout: float = 600.0,
                       dossiers: dict | None = None, string_data: dict | None = None,
                       tissue_data: dict | None = None, drug_data: dict | None = None) -> str:
    """Scientific report for specific genes via the active LLM backend, grounded in the
    computed facts AND (optionally) retrieved literature + protein/drug annotations (RAG).
    `dossiers` maps gene -> knowledge.dossier() output.

    `tissue_data` is deliberately NOT included in the LLM's own prompt (see the comment
    where it's applied to the deterministic comparative table below) -- it is still used,
    unmodified, to fill that table's "Tissue gates" column after the LLM call returns.
    """
    facts = [gene_facts(ctx, g) for g in genes]
    facts = [f for f in facts if f is not None]
    if not facts:
        raise ValueError("none of the requested genes are in the ranked table")

    # attach retrieved knowledge (trimmed) per gene when available
    if dossiers:
        for f in facts:
            d = dossiers.get(f["gene"])
            if not d:
                continue
            up = d.get("uniprot") or {}
            ch = d.get("chembl") or {}
            f["protein"] = {
                "name": up.get("protein_name", ""),
                "function": (up.get("function", "") or "")[:400],
                "domains": up.get("domains", []),
                "active_sites": up.get("active_sites", []),
                "binding_sites": up.get("binding_sites", []),
                "subcellular_location": up.get("subcellular_location", []),
                "n_pdb": up.get("n_pdb", 0),
            }
            f["known_drugs"] = {
                "n_mechanisms": ch.get("n_mechanisms", 0),
                "mechanisms": [m.get("moa") for m in ch.get("mechanisms", [])][:6],
            }
            f["druggability"] = d.get("druggability", "")
            f["literature"] = [
                {"title": p["title"], "year": p["year"], "pmid": p["pmid"]}
                for p in d.get("literature", [])[:6]
            ]

    if string_data:
        for f in facts:
            sd = string_data.get(f["gene"])
            if not sd:
                continue
            f["string_db"] = {
                "functional_max": sd["functional"]["max"],
                "functional_mean": sd["functional"]["mean"],
                "functional_best_seed": sd["functional"]["best_seed"],
                "physical_max": sd["physical"]["max"],
                "physical_mean": sd["physical"]["mean"],
                "physical_best_seed": sd["physical"]["best_seed"],
                "n_seeds_functional_ge_cutoff": sd["functional"]["n_seeds_ge_cutoff"],
                "n_seeds_physical_ge_cutoff": sd["physical"]["n_seeds_ge_cutoff"],
                "cutoff": sd["cutoff"],
            }

    # NOTE: tissue_data is deliberately NOT attached to the JSON facts sent to the LLM's
    # prompt (unlike dossiers/string_data/drug_data below) -- the pass_normal/pass_pancancer/
    # pass_target gates mix a lenient median-based gate with a stricter max-based ranking
    # score (see comparative_facts_table()'s docstring and tissue_specificity.py), a nuance
    # too easy for an LLM to flatten into an overconfident or misleading sentence. tissue_data
    # is still used, deterministically, in the appended comparative table's "Tissue gates"
    # column below -- that's a plain code-computed pass/fail, not a model interpretation.
    if drug_data:
        for f in facts:
            dd = drug_data.get(f["gene"])
            if not dd:
                continue
            f["drug_targets"] = dd

    ctx_facts = {
        "age_group": ctx.spec.age_group, "disease": ctx.spec.disease,
        "seed_source": ctx.seeds.source, "seed_genes": ctx.seeds.genes,
        "essentiality_threshold": ctx.params.ess_thresh,
        "expr_tpm_threshold": ctx.params.expr_tpm,
        "n_cohort_lines": ctx.ess_denom, "genes": facts,
    }
    has_rag = bool(dossiers or string_data or drug_data)
    system = (
        "You are a cautious cancer-genomics research assistant producing a STRUCTURED "
        "REPORT (not free-form prose) from ALREADY-COMPUTED DepMap metrics for one or more "
        "genes, each found co-expressed with the INO80/SRCAP chromatin-remodeling seed "
        "complex, with per-cell-line essentiality in the selected cancer cohort. "
        + ("Each gene may also include retrieved UniProt protein annotation, ChEMBL "
           "known-drug mechanisms, Europe PMC literature (title/year/PMID), STRING-db "
           "association scores to the seed proteins (functional and physical, 0-1), and a "
           "drug_targets block (known drugs/compounds for this gene from ChEMBL/DGIdb/Open "
           "Targets, each with a clinical-maturity + cancer-relevance score and label). "
           if has_rag else "")
        + "Structure your response as EXACTLY these six markdown sections, each starting "
        "with a '## ' heading, in this exact order, each section addressing ALL requested "
        "genes together by name (do not repeat the section list once per gene):\n"
        "## Executive Summary\n"
        "One short paragraph per gene: is it a plausible cohort-selective dependency worth "
        "pursuing? State the headline selectivity_delta for each gene.\n"
        "## Cohort-Selective Dependency Evidence\n"
        "For each gene: pct_essential_cohort vs pct_essential_pan, selectivity_delta, and "
        "whether common_essential=true (a poor selective target).\n"
        "## Protein Complex / STRING Association\n"
        "For each gene: does STRING support a physical link (high physical_max = likely "
        "direct complex member) or only functional/pathway-level association? If no STRING "
        "data was provided for a gene, say so explicitly rather than omitting it.\n"
        "## Drug-Design & Therapeutic Relevance\n"
        "For each gene: protein annotation (domains, active/binding sites, PDB structures) "
        "and, if a drug_targets block is present, its top-scoring known drug(s) with their "
        "Score_Label. If no annotation/drug data was provided, say so explicitly.\n"
        "## Literature Evidence\n"
        "For each gene with a literature list, cite papers ONLY by the PMIDs present in the "
        "JSON. If none were retrieved for a gene, say so explicitly.\n"
        "## Experimental Validation & Testing Guidelines\n"
        "For each gene, propose 2-4 concrete NEXT laboratory steps grounded only in what's "
        "actually in the JSON (e.g. orthogonal CRISPR/RNAi/shRNA knockdown validation in "
        "the named cohort, specific assay types implied by the annotation such as "
        "co-crystallization or fragment screening if PDB structures exist, viability/"
        "apoptosis readouts, comparison against other genes in this same request as "
        "controls). Do not invent cell line names, reagents, vendors, or assays not "
        "implied by the data.\n"
        "Do NOT add a seventh section or a table \u2014 a comparative summary table is "
        "appended separately after your response, verbatim, by the application. "
        "STRICT RULES: use ONLY the genes, numbers, annotations, drugs, STRING scores, "
        "tissue names, and papers present in the JSON; do NOT invent numbers, gene names, "
        "tissue/cohort names, drugs, PDB ids, pathways as fact, or PMIDs. If a field is "
        "empty say the evidence is absent. Only cite PMIDs that appear in the JSON. Treat "
        "common_essential=true as a poor selective target. Do not claim clinical utility."
    )
    prompt = f"DATA (JSON):\n{json.dumps(ctx_facts, indent=2)}\n\nWrite the report now."
    # Some Claude deployments (observed: an extended-thinking "5"-generation model behind
    # Azure AI Foundry) spend part of max_tokens on an internal reasoning pass before any
    # visible text -- with a small fixed budget this can silently consume the ENTIRE budget,
    # returning an empty string with a normal 200 OK (stop_reason=max_tokens), not an error.
    # Empirically verified: a fixed 900 was too small even for 3 genes with literature/STRING
    # context on such a deployment (0 chars back). The seven-section structured report is
    # longer per gene than the original single-paragraph ask, so the budget is scaled up
    # further and re-verified live. Capped so an unusually long gene list can't runaway.
    max_tok = min(6000, max(2200, 1500 * len(facts)))
    text = llm_complete(system, prompt, model=model, timeout=timeout, max_tokens=max_tok)
    if not text.strip():
        return text
    headers, rows = comparative_facts_table(ctx, genes, dossiers, string_data,
                                            tissue_data, drug_data)
    table_md = comparative_table_markdown(headers, rows)
    return f"{text}\n\n{table_md}"
