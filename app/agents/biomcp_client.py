"""BioMCP integration: Reactome/KEGG pathway membership + Human Protein Atlas
(HPA) tissue-expression/localization context, per gene, for the Insights tab.

BioMCP (https://biomcp.org, MIT-licensed, github.com/genomoncology/biomcp) is a
standalone CLI/MCP-server binary, not a Python package -- it is installed
separately (`uv tool install biomcp-cli` or the biomcp.org install script) and
is entirely OPTIONAL. This module shells out to the installed `biomcp` binary
(`biomcp get gene <SYMBOL> pathways hpa -j`) rather than reimplementing a
client for its ~30 upstream sources; that command grammar and its JSON schema
were verified directly against a live install this session (real output for
BRD7 captured and used to shape the parsing below), not guessed from docs.

If the binary isn't installed, a gene isn't resolved, or the underlying APIs
are unreachable, every function here returns None / a dict with
`available: False` -- never a fabricated value. This mirrors the app's
established pattern for other optional/precomputed evidence (e.g. Model
Predictions' "Not generated/cached" fallback).
"""
from __future__ import annotations

import json
import shutil
import subprocess

# Fallback search locations for a couple of common install methods, in case
# the launching process's PATH doesn't include them (e.g. `uv tool install`
# defaults to ~/.local/bin, which isn't always on PATH for GUI-launched apps).
_FALLBACK_PATHS = [
    "~/.local/bin/biomcp",
    "~/.cargo/bin/biomcp",
]


def _find_biomcp() -> str | None:
    """Locate the biomcp binary, or None if it isn't installed anywhere we look."""
    found = shutil.which("biomcp")
    if found:
        return found
    import os
    for p in _FALLBACK_PATHS:
        expanded = os.path.expanduser(p)
        if os.path.isfile(expanded) and os.access(expanded, os.X_OK):
            return expanded
    return None


def is_available() -> bool:
    """Whether the biomcp binary can be found at all (says nothing about
    network reachability of its upstream sources)."""
    return _find_biomcp() is not None


def gene_pathways_hpa(gene: str, timeout: float = 30.0) -> dict | None:
    """Reactome/KEGG pathway membership + Human Protein Atlas tissue/subcellular
    context for one gene symbol, via `biomcp get gene <SYMBOL> pathways hpa -j`.

    Returns None if biomcp isn't installed, the gene doesn't resolve, the
    upstream APIs are unreachable, or the binary's output can't be parsed as
    the expected JSON shape -- callers should treat None as "no evidence
    available" and say so honestly, not silently omit the section.

    On success, returns:
        {
            "pathways": [{"source": "Reactome"|"KEGG", "id": str, "name": str}, ...],
            "hpa": {
                "protein_summary": str, "rna_summary": str, "reliability": str,
                "subcellular_main_location": [str, ...],
                "subcellular_additional_location": [str, ...],
                "tissues": [{"tissue": str, "level": str}, ...],
            } | None,   # None if HPA has no data for this gene
        }
    """
    binary = _find_biomcp()
    if binary is None:
        return None
    try:
        proc = subprocess.run(
            [binary, "get", "gene", gene, "pathways", "hpa", "-j"],
            capture_output=True, text=True, timeout=timeout,
        )
    except (subprocess.TimeoutExpired, OSError):
        return None
    if proc.returncode != 0 or not proc.stdout.strip():
        return None
    try:
        data = json.loads(proc.stdout)
    except json.JSONDecodeError:
        return None
    if not isinstance(data, dict):
        return None

    pathways = []
    for p in data.get("pathways") or []:
        if not isinstance(p, dict):
            continue
        pathways.append({
            "source": p.get("source", ""),
            "id": p.get("id", ""),
            "name": p.get("name", ""),
        })

    hpa_raw = data.get("hpa")
    hpa = None
    if isinstance(hpa_raw, dict):
        hpa = {
            "protein_summary": hpa_raw.get("protein_summary", ""),
            "rna_summary": hpa_raw.get("rna_summary", ""),
            "reliability": hpa_raw.get("reliability", ""),
            "subcellular_main_location": hpa_raw.get("subcellular_main_location", []) or [],
            "subcellular_additional_location":
                hpa_raw.get("subcellular_additional_location", []) or [],
            "tissues": [
                {"tissue": t.get("tissue", ""), "level": t.get("level", "")}
                for t in (hpa_raw.get("tissues") or []) if isinstance(t, dict)
            ],
        }

    return {"pathways": pathways, "hpa": hpa}
