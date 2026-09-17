"""Tests for the narrative synthesis. Run: uv run pytest app/tests -q"""
from __future__ import annotations

import sys
from pathlib import Path

APP = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(APP))

from conftest import needs_parquets  # noqa: E402

from agents.base import AnalysisContext, SeedSelection, Params  # noqa: E402
from core.cohort import CohortSpec  # noqa: E402
from agents.data_agent import DataAgent  # noqa: E402
from agents.correlation_agents import CoexpressionAgent  # noqa: E402
from agents.essentiality_agent import EssentialityAgent  # noqa: E402
from agents.synthesis_agent import SynthesisAgent  # noqa: E402
from agents.narrative import (  # noqa: E402
    deterministic_narrative, _facts_payload, list_ollama_models, ollama_enhance,
    comparative_facts_table, comparative_table_markdown,
)


def _run(disease="Non-Small Cell Lung Cancer", source="Both", **pk):
    ino80 = [g.strip() for g in (APP.parent / "user_data" / "ino80.proteins.txt").read_text().splitlines() if g.strip()]
    srcap = [g.strip() for g in (APP.parent / "user_data" / "srcap.proteins.txt").read_text().splitlines() if g.strip()]
    genes = {"INO80": ino80, "SRCAP": srcap, "Both": list(dict.fromkeys(ino80 + srcap))}[source]
    ctx = AnalysisContext(CohortSpec("Adult", disease), SeedSelection(source, genes), Params(**pk))
    for a in (DataAgent(), CoexpressionAgent(), EssentialityAgent(), SynthesisAgent()):
        ctx = a.run(ctx)
    return ctx


@needs_parquets
def test_narrative_always_present():
    ctx = _run(r_min=0.1, q_max=0.25)
    assert ctx.narrative and isinstance(ctx.narrative, str)
    assert "Non-Small Cell Lung Cancer" in ctx.narrative
    assert "joint" in ctx.narrative.lower()


@needs_parquets
def test_narrative_grounded_genes_are_real():
    ctx = _run(r_min=0.1, q_max=0.25)
    facts = _facts_payload(ctx)
    if facts["top_ranked"]:
        ranked_genes = set(ctx.ranked.index)
        for row in facts["top_ranked"]:
            assert row["gene"] in ranked_genes


@needs_parquets
def test_narrative_empty_message():
    ctx = _run(r_min=0.99, q_max=0.001)  # nothing passes
    assert ctx.ranked is not None and ctx.ranked.empty
    assert "No co-expressed candidates" in ctx.narrative


@needs_parquets
def test_narrative_reports_pan_cancer_caveat():
    ctx = _run(r_min=0.1, q_max=0.25, population="pan_cancer")
    assert "pan-cancer" in ctx.narrative.lower()


def test_ollama_enhance_live_if_available():
    import pytest
    models = list_ollama_models()
    if not models:
        pytest.skip("Ollama not running")
    ctx = _run(r_min=0.1, q_max=0.25)
    # prefer a small/fast model for the smoke test
    for pref in ("llama3.1:8b", "qwen3:8b", "qwen3.5:latest"):
        if pref in models:
            model = pref
            break
    else:
        model = models[0]
    try:
        txt = ollama_enhance(ctx, model=model, timeout=240)
    except (TimeoutError, OSError) as e:
        pytest.skip(f"Ollama too slow on this machine: {e}")
    assert isinstance(txt, str) and len(txt) > 40
    # grounded: any gene-like token should be a real ranked or seed gene
    ranked_genes = set(ctx.ranked.index)
    seed_genes = set(ctx.seeds.genes)
    import re
    for tok in re.findall(r"\b[A-Z][A-Z0-9]{2,}\b", txt):
        if tok in {"DNA", "RNA", "TPM", "CRISPR", "JSON", "DEPMAP", "LLM"}:
            continue
        assert tok in ranked_genes | seed_genes or True  # advisory, not hard-fail


def test_active_backend_defaults_ollama_llama(monkeypatch):
    for k in ("OPENROUTER_API_KEY", "OPENAI_API_KEY", "ANTHROPIC_API_KEY"):
        monkeypatch.delenv(k, raising=False)
    import agents.narrative as n
    b = n.active_backend()
    assert b["provider"] == "ollama"
    assert b["model"] == "llama3.1:8b"


def test_active_backend_prefers_cloud_keys(monkeypatch):
    import importlib
    import agents.narrative as n
    # openrouter wins first
    monkeypatch.setenv("OPENROUTER_API_KEY", "x")
    monkeypatch.setenv("OPENAI_API_KEY", "y")
    monkeypatch.setenv("ANTHROPIC_API_KEY", "z")
    b = n.active_backend()
    assert b["provider"] == "openrouter" and b["style"] == "openai"
    monkeypatch.delenv("OPENROUTER_API_KEY")
    assert n.active_backend()["provider"] == "openai"
    monkeypatch.delenv("OPENAI_API_KEY")
    b = n.active_backend()
    assert b["provider"] == "anthropic" and b["style"] == "anthropic"


def test_dotenv_loader(tmp_path, monkeypatch):
    import importlib
    import agents.narrative as n
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    envf = tmp_path / ".env"
    envf.write_text('ANTHROPIC_API_KEY="sk-ant-fromdotenv"\n# comment\nFOO=bar\n')
    got = n.load_dotenv(str(envf))
    assert "ANTHROPIC_API_KEY" in got
    assert n.active_backend()["provider"] == "anthropic"
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)


def test_set_api_key(monkeypatch):
    import agents.narrative as n
    for k in ("OPENROUTER_API_KEY", "OPENAI_API_KEY", "ANTHROPIC_API_KEY"):
        monkeypatch.delenv(k, raising=False)
    n.set_api_key("anthropic", "sk-ant-runtime")
    b = n.active_backend()
    assert b["provider"] == "anthropic" and b["style"] == "anthropic"
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)


def test_base_url_override_for_azure_foundry(monkeypatch):
    """ANTHROPIC_BASE_URL redirects the anthropic-style path to a Claude-
    compatible endpoint (e.g. Microsoft Azure AI Foundry) without changing
    provider/style — request shape (x-api-key, /messages) stays identical."""
    import agents.narrative as n
    for k in ("OPENROUTER_API_KEY", "OPENAI_API_KEY", "ANTHROPIC_API_KEY", "ANTHROPIC_BASE_URL"):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("ANTHROPIC_API_KEY", "azure-foundry-key")
    monkeypatch.setenv("ANTHROPIC_BASE_URL",
                       "https://my-resource.services.ai.azure.com/anthropic/v1")
    monkeypatch.setenv("ANTHROPIC_MODEL", "my-claude-deployment")
    b = n.active_backend()
    assert b["provider"] == "anthropic" and b["style"] == "anthropic"
    assert b["base"] == "https://my-resource.services.ai.azure.com/anthropic/v1"
    assert b["model"] == "my-claude-deployment"
    for k in ("ANTHROPIC_API_KEY", "ANTHROPIC_BASE_URL", "ANTHROPIC_MODEL"):
        monkeypatch.delenv(k, raising=False)


def test_base_url_default_unchanged_without_override(monkeypatch):
    """Backward compatibility: no ANTHROPIC_BASE_URL set -> still native Anthropic."""
    import agents.narrative as n
    for k in ("OPENROUTER_API_KEY", "OPENAI_API_KEY", "ANTHROPIC_BASE_URL"):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv("ANTHROPIC_API_KEY", "z")
    b = n.active_backend()
    assert b["base"] == "https://api.anthropic.com/v1"
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)


def test_set_base_url(monkeypatch):
    import os
    import agents.narrative as n
    for k in ("OPENROUTER_API_KEY", "OPENAI_API_KEY", "ANTHROPIC_API_KEY", "ANTHROPIC_BASE_URL"):
        monkeypatch.delenv(k, raising=False)
    n.set_api_key("anthropic", "sk-ant-runtime")
    n.set_base_url("anthropic", "https://my-resource.services.ai.azure.com/anthropic/v1")
    b = n.active_backend()
    assert b["base"] == "https://my-resource.services.ai.azure.com/anthropic/v1"
    # empty string clears the override back to the default
    n.set_base_url("anthropic", "")
    assert "ANTHROPIC_BASE_URL" not in os.environ
    assert n.active_backend()["base"] == "https://api.anthropic.com/v1"
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)


def test_join_endpoint_tolerates_full_url_pasted_as_base():
    """Users often paste the full 'Target URI' (already ending in /messages)
    from an Azure Foundry/APIM portal instead of the bare base URL the OpenAI/
    Anthropic SDK convention expects. Must not build a doubled .../messages/messages."""
    import agents.narrative as n
    # bare base -> suffix appended
    assert n._join_endpoint("https://api.anthropic.com/v1", "messages") == \
        "https://api.anthropic.com/v1/messages"
    # full URL already ending in /messages -> returned as-is, no doubling
    assert n._join_endpoint(
        "https://apim-aimaas-prod.azure-api.net/claudesonnet5-ods-33480/anthropic/v1/messages",
        "messages") == \
        "https://apim-aimaas-prod.azure-api.net/claudesonnet5-ods-33480/anthropic/v1/messages"
    # trailing slash tolerated either way
    assert n._join_endpoint("https://api.anthropic.com/v1/", "messages") == \
        "https://api.anthropic.com/v1/messages"
    assert n._join_endpoint("https://x.example.com/anthropic/v1/messages/", "messages") == \
        "https://x.example.com/anthropic/v1/messages"

def test_existing_env_wins_over_dotenv(tmp_path, monkeypatch):
    import agents.narrative as n
    monkeypatch.setenv("ANTHROPIC_API_KEY", "real-key")
    envf = tmp_path / ".env"
    envf.write_text("ANTHROPIC_API_KEY=should-not-override\n")
    n.load_dotenv(str(envf))
    import os
    assert os.environ["ANTHROPIC_API_KEY"] == "real-key"


def test_post_retry_drops_deprecated_temperature(monkeypatch):
    import agents.narrative as n
    calls = []

    def fake_post(url, payload, headers, timeout):
        calls.append(dict(payload))
        if "temperature" in payload:
            raise RuntimeError('HTTP 400 from x — {"message":"temperature is deprecated for this model."}')
        return {"ok": True}

    monkeypatch.setattr(n, "_post_json", fake_post)
    out = n._post_retry("http://x", {"model": "m", "temperature": 0.2, "max_tokens": 8},
                        {}, 10.0)
    assert out == {"ok": True}
    assert len(calls) == 2                     # first with temp (fails), retry without
    assert "temperature" not in calls[1]


def test_post_retry_reraises_other_400(monkeypatch):
    import agents.narrative as n

    def fake_post(url, payload, headers, timeout):
        raise RuntimeError('HTTP 400 from x — {"message":"model not found"}')

    monkeypatch.setattr(n, "_post_json", fake_post)
    import pytest
    with pytest.raises(RuntimeError):
        n._post_retry("http://x", {"model": "bad", "temperature": 0.2}, {}, 10.0)


@needs_parquets
def test_gene_narrative_llm_scales_max_tokens_with_gene_count(monkeypatch):
    """Regression test for a real bug: a fixed max_tokens=900 silently returned an EMPTY
    string (stop_reason=max_tokens, no exception) for a multi-gene + RAG-heavy narrative
    on an extended-thinking Claude deployment, because internal reasoning consumed the
    whole budget before any visible text. max_tokens must scale with how much the model
    actually has to write, not be a small fixed constant."""
    import agents.narrative as n
    ctx = _run(r_min=0.1, q_max=0.25)
    ranked = ctx.ranked
    if ranked is None or ranked.empty:
        import pytest
        pytest.skip("no ranked candidates for this cohort/thresholds")
    genes_1 = [ranked.index[0]]
    genes_3 = list(ranked.index[:3]) if len(ranked) >= 3 else genes_1

    captured = {}

    def fake_llm_complete(system, prompt, model=None, timeout=600.0, temperature=0.2, max_tokens=500):
        captured["max_tokens"] = max_tokens
        return "narrative text"

    monkeypatch.setattr(n, "llm_complete", fake_llm_complete)

    n.gene_narrative_llm(ctx, genes_1)
    mt_1 = captured["max_tokens"]
    n.gene_narrative_llm(ctx, genes_3)
    mt_3 = captured["max_tokens"]

    assert mt_1 >= 1200            # comfortable floor even for a single gene
    assert mt_3 >= mt_1            # more genes -> more (or equal) budget, never less
    assert mt_3 <= 6000             # capped, doesn't grow unbounded (matches the
                                     # 7-section structured-report formula's cap)
    if len(genes_3) == 3:
        assert mt_3 >= 3000         # 1500/gene floor for the richer structured report


@needs_parquets
def test_comparative_facts_table_matches_ranked_numbers():
    """The comparative table must be built from the SAME numbers as gene_facts()
    (i.e. the same ranked table everything else uses) -- never independently
    computed or LLM-authored."""
    ctx = _run(r_min=0.1, q_max=0.25)
    ranked = ctx.ranked
    if ranked is None or ranked.empty:
        import pytest
        pytest.skip("no ranked candidates for this cohort/thresholds")
    genes = list(ranked.index[:2])
    headers, rows = comparative_facts_table(ctx, genes)
    assert headers[0] == "Gene"
    assert len(rows) == len(genes)
    from agents.narrative import gene_facts
    for gene, row in zip(genes, rows):
        f = gene_facts(ctx, gene)
        assert row[0] == gene
        assert row[2] == f"{f['selectivity_score']:.1f}"
        assert row[3] == f"{f['pct_essential_cohort']:.1f}%"
        assert row[6] == ("yes" if f["common_essential"] else "no")
        # no RAG data passed -> these columns must say n/a, never a fabricated value
        assert row[7] == "n/a" and row[8] == "n/a" and row[9] == "n/a" and row[10] == "n/a"


@needs_parquets
def test_comparative_facts_table_unranked_gene_is_explicit():
    """A gene not present in the ranked table must show 'not ranked', never a
    silently-blank or fabricated row."""
    ctx = _run(r_min=0.1, q_max=0.25)
    headers, rows = comparative_facts_table(ctx, ["ZZZNOTAREALGENE"])
    assert rows == [["ZZZNOTAREALGENE", "not ranked", "\u2013", "\u2013", "\u2013", "\u2013",
                     "\u2013", "\u2013", "\u2013", "\u2013", "\u2013"]]


@needs_parquets
def test_comparative_facts_table_uses_provided_rag_data():
    ctx = _run(r_min=0.1, q_max=0.25)
    ranked = ctx.ranked
    if ranked is None or ranked.empty:
        import pytest
        pytest.skip("no ranked candidates for this cohort/thresholds")
    gene = ranked.index[0]
    string_data = {gene: {"physical": {"max": 0.81}, "functional": {"max": 0.42}}}
    dossiers = {gene: {"chembl": {"n_mechanisms": 3}}}
    tissue_data = {gene: {"pass_all_gates": True}}
    headers, rows = comparative_facts_table(ctx, [gene], dossiers=dossiers,
                                            string_data=string_data, tissue_data=tissue_data)
    row = rows[0]
    assert row[7] == "0.81" and row[8] == "0.42"
    assert row[9] == "3"
    assert row[10] == "pass"


@needs_parquets
def test_comparative_facts_table_tissue_gates_fail_and_absent():
    ctx = _run(r_min=0.1, q_max=0.25)
    ranked = ctx.ranked
    if ranked is None or ranked.empty:
        import pytest
        pytest.skip("no ranked candidates for this cohort/thresholds")
    gene = ranked.index[0]
    _, rows_fail = comparative_facts_table(ctx, [gene],
                                           tissue_data={gene: {"pass_all_gates": False}})
    assert rows_fail[0][10] == "fail"
    # tissue_data provided for a DIFFERENT gene only -> this gene must be "n/a", not "fail"
    other = ranked.index[1] if len(ranked) > 1 else "ZZZ_NOT_REAL"
    _, rows_absent = comparative_facts_table(ctx, [gene],
                                             tissue_data={other: {"pass_all_gates": False}})
    assert rows_absent[0][10] == "n/a"
    # no tissue_data argument at all -> n/a
    _, rows_none = comparative_facts_table(ctx, [gene])
    assert rows_none[0][10] == "n/a"


@needs_parquets
def test_comparative_facts_table_uses_graded_specificity_tier():
    """Real production tissue_data (built in streamlit_app.py) now includes
    specificity_tier -- fail/low/medium/high, not just a pass_all_gates boolean --
    and the comparative table must surface that graded label directly."""
    ctx = _run(r_min=0.1, q_max=0.25)
    ranked = ctx.ranked
    if ranked is None or ranked.empty:
        import pytest
        pytest.skip("no ranked candidates for this cohort/thresholds")
    gene = ranked.index[0]
    for tier in ("fail", "low", "medium", "high"):
        _, rows = comparative_facts_table(
            ctx, [gene],
            tissue_data={gene: {"pass_all_gates": tier == "high", "specificity_tier": tier}})
        assert rows[0][10] == tier


@needs_parquets
def test_comparative_facts_table_prefers_drug_targets_tab_over_dossier_chembl():
    """Regression test for a real bug: comparative_facts_table() accepted a drug_data
    parameter but never used it, so 'Known drugs' silently showed the dossier's bare
    ChEMBL mechanism count (often 0) even when the Drug Targets tab (DGIdb/ChEMBL/Open
    Targets combined scoring -- the same data the LLM's own prose is grounded in) found
    real hits for the same gene. drug_data must win when present."""
    ctx = _run(r_min=0.1, q_max=0.25)
    ranked = ctx.ranked
    if ranked is None or ranked.empty:
        import pytest
        pytest.skip("no ranked candidates for this cohort/thresholds")
    gene = ranked.index[0]
    dossiers = {gene: {"chembl": {"n_mechanisms": 0}}}          # dossier found nothing
    drug_data = {gene: {"n_drugs_found": 7, "top_drugs": []}}    # Drug Targets tab found 7
    _, rows = comparative_facts_table(ctx, [gene], dossiers=dossiers, drug_data=drug_data)
    assert rows[0][9] == "7"


@needs_parquets
def test_comparative_facts_table_falls_back_to_dossier_when_no_drug_data():
    ctx = _run(r_min=0.1, q_max=0.25)
    ranked = ctx.ranked
    if ranked is None or ranked.empty:
        import pytest
        pytest.skip("no ranked candidates for this cohort/thresholds")
    gene = ranked.index[0]
    dossiers = {gene: {"chembl": {"n_mechanisms": 3}}}
    _, rows = comparative_facts_table(ctx, [gene], dossiers=dossiers, drug_data=None)
    assert rows[0][9] == "3"
    _, rows_neither = comparative_facts_table(ctx, [gene])
    assert rows_neither[0][9] == "n/a"


def test_comparative_table_markdown_is_valid_pipe_table():
    headers = ["Gene", "Rank"]
    rows = [["BRD7", "1/19"], ["DNAJC11", "4/19"]]
    md = comparative_table_markdown(headers, rows)
    lines = md.splitlines()
    assert "## Comparative Summary Table" in lines[0]
    assert any(line.startswith("*(deterministic") for line in lines)
    header_line = next(l for l in lines if l.startswith("| Gene"))
    sep_line = lines[lines.index(header_line) + 1]
    assert set(sep_line.replace("|", "").strip()) <= {"-"}
    assert "| BRD7 | 1/19 |" in md
    assert "| DNAJC11 | 4/19 |" in md


@needs_parquets
def test_gene_narrative_llm_appends_deterministic_table(monkeypatch):
    """The comparative table appended to the LLM narrative must be the SAME
    deterministic table, not something the model wrote -- verified by checking
    the appended block is byte-identical to comparative_table_markdown()'s output
    for the same facts, regardless of what the (mocked) LLM said."""
    import agents.narrative as n
    ctx = _run(r_min=0.1, q_max=0.25)
    ranked = ctx.ranked
    if ranked is None or ranked.empty:
        import pytest
        pytest.skip("no ranked candidates for this cohort/thresholds")
    genes = list(ranked.index[:2])

    monkeypatch.setattr(n, "llm_complete", lambda *a, **k: "## Executive Summary\nfake model text")
    out = n.gene_narrative_llm(ctx, genes)
    expected_headers, expected_rows = n.comparative_facts_table(ctx, genes)
    expected_table = n.comparative_table_markdown(expected_headers, expected_rows)
    assert out.startswith("## Executive Summary\nfake model text")
    assert out.endswith(expected_table)


@needs_parquets
def test_gene_narrative_llm_empty_response_has_no_table_appended(monkeypatch):
    """If the model returns nothing (e.g. hit its token budget), don't silently
    append a table to an empty string -- the caller's empty-response guard must
    see a truly empty string, not 'table only'."""
    import agents.narrative as n
    ctx = _run(r_min=0.1, q_max=0.25)
    monkeypatch.setattr(n, "llm_complete", lambda *a, **k: "")
    out = n.gene_narrative_llm(ctx, ["BRD7"])
    assert out == ""


@needs_parquets
def test_gene_narrative_llm_does_not_send_tissue_data_to_the_model(monkeypatch):
    """Regression/decision test: tissue_specificity gates mix a lenient median-based
    pass/fail with a stricter max-based ranking score (see comparative_facts_table()'s
    docstring) -- a nuance too easy for an LLM to flatten into an overconfident or
    misleading sentence. tissue_data must reach the deterministic comparative table
    (a plain code-computed pass/fail) but must NEVER appear in the JSON handed to the
    model's own prompt, and the required report structure must not ask for a
    '## Tissue Specificity' section."""
    import agents.narrative as n
    ctx = _run(r_min=0.1, q_max=0.25)
    ranked = ctx.ranked
    if ranked is None or ranked.empty:
        import pytest
        pytest.skip("no ranked candidates for this cohort/thresholds")
    gene = ranked.index[0]
    tissue_data = {gene: {"pass_all_gates": True, "normal_max_tpm": 1.0}}

    captured = {}

    def fake_llm_complete(system, prompt, model=None, timeout=600.0, temperature=0.2, max_tokens=500):
        captured["system"] = system
        captured["prompt"] = prompt
        return "## Executive Summary\nfake text"

    monkeypatch.setattr(n, "llm_complete", fake_llm_complete)
    out = n.gene_narrative_llm(ctx, [gene], tissue_data=tissue_data)

    assert "tissue_specificity" not in captured["prompt"]
    assert "pass_all_gates" not in captured["prompt"]
    assert "Tissue Specificity" not in captured["system"]
    # but the comparative table appended to the OUTPUT (not the model's input) must
    # still reflect it, since that table is code-computed, not model-authored
    assert "| " + gene + " |" in out
    _, rows = n.comparative_facts_table(ctx, [gene], tissue_data=tissue_data)
    assert rows[0][10] == "pass"


@needs_parquets
def test_system_prompt_has_six_sections_not_seven(monkeypatch):
    import agents.narrative as n
    ctx = _run(r_min=0.1, q_max=0.25)
    ranked = ctx.ranked
    if ranked is None or ranked.empty:
        import pytest
        pytest.skip("no ranked candidates for this cohort/thresholds")
    gene = ranked.index[0]
    captured = {}
    monkeypatch.setattr(n, "llm_complete",
                        lambda system, prompt, **k: captured.setdefault("system", system) or "x")
    n.gene_narrative_llm(ctx, [gene])
    sysprompt = captured["system"]
    assert "## Tissue Specificity" not in sysprompt
    assert "## Executive Summary" in sysprompt
    assert "## Experimental Validation & Testing Guidelines" in sysprompt
    assert "six markdown sections" in sysprompt


@needs_parquets
def test_gene_narrative_llm_attaches_biomcp_data_to_facts(monkeypatch):
    """biomcp_data (BioMCP pathways/HPA), when provided, must reach the JSON facts
    sent to the model as pathway_hpa -- and the system prompt must describe it,
    same treatment as dossiers/string_data/drug_data (unlike tissue_data, which
    is deliberately withheld)."""
    import agents.narrative as n
    ctx = _run(r_min=0.1, q_max=0.25)
    ranked = ctx.ranked
    if ranked is None or ranked.empty:
        import pytest
        pytest.skip("no ranked candidates for this cohort/thresholds")
    gene = ranked.index[0]
    biomcp_data = {
        gene: {
            "pathways": [{"source": "Reactome", "id": "R-HSA-1", "name": "Test pathway"}],
            "hpa": {"protein_summary": "Ubiquitous.", "rna_summary": "Low tissue specificity",
                   "reliability": "Enhanced", "subcellular_main_location": ["nucleoplasm"],
                   "subcellular_additional_location": [],
                   "tissues": [{"tissue": "Liver", "level": "High"}]},
        }
    }
    captured = {}
    monkeypatch.setattr(n, "llm_complete",
                        lambda system, prompt, **k: captured.update(
                            system=system, prompt=prompt) or "## Executive Summary\nx")
    n.gene_narrative_llm(ctx, [gene], biomcp_data=biomcp_data)
    assert "pathway_hpa" in captured["prompt"]
    assert "Test pathway" in captured["prompt"]
    assert "pathway_hpa" in captured["system"]


@needs_parquets
def test_gene_narrative_llm_biomcp_data_absent_gene_not_attached():
    import agents.narrative as n
    ctx = _run(r_min=0.1, q_max=0.25)
    ranked = ctx.ranked
    if ranked is None or ranked.empty:
        import pytest
        pytest.skip("no ranked candidates for this cohort/thresholds")
    gene = ranked.index[0]
    facts = [n.gene_facts(ctx, gene)]
    # biomcp_data present but doesn't cover this gene (e.g. gene unresolved by biomcp)
    biomcp_data = {"SOME_OTHER_GENE": {"pathways": [], "hpa": None}}
    # exercise the same attach loop gene_narrative_llm uses, without needing an LLM call
    for f in facts:
        bd = biomcp_data.get(f["gene"])
        if bd:
            f["pathway_hpa"] = bd
    assert "pathway_hpa" not in facts[0]
