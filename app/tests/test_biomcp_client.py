"""Tests for the BioMCP integration (Reactome/KEGG pathways + HPA tissue context).
Run: uv run pytest app/tests -q

Offline: parsing logic tested against a REAL captured fixture (BRD7's actual
`biomcp get gene BRD7 pathways hpa -j` output, captured live this session --
see fixtures/biomcp_gene_brd7_pathways_hpa.json), not a guessed schema.
Network+binary-gated: a live call, skipped cleanly if biomcp isn't installed
or the sandbox is offline.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

APP = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(APP))

from agents import biomcp_client as bc  # noqa: E402

FIXTURE = Path(__file__).parent / "fixtures" / "biomcp_gene_brd7_pathways_hpa.json"


def _fake_run_factory(stdout: str = "", returncode: int = 0):
    def _fake_run(cmd, capture_output, text, timeout):
        return subprocess.CompletedProcess(cmd, returncode, stdout=stdout, stderr="")
    return _fake_run


def test_parses_real_captured_fixture(monkeypatch):
    """The exact JSON shape BioMCP returned for a real `biomcp get gene BRD7
    pathways hpa -j` call this session -- not a guessed/synthetic schema."""
    raw = FIXTURE.read_text()
    monkeypatch.setattr(bc, "_find_biomcp", lambda: "/fake/biomcp")
    monkeypatch.setattr(bc.subprocess, "run", _fake_run_factory(stdout=raw))

    out = bc.gene_pathways_hpa("BRD7")
    assert out is not None
    assert len(out["pathways"]) == 13
    assert out["pathways"][0] == {
        "source": "KEGG", "id": "hsa05225",
        "name": "Hepatocellular carcinoma - Homo sapiens (human)",
    }
    assert any(p["source"] == "Reactome" for p in out["pathways"])
    assert out["hpa"]["reliability"] == "Enhanced"
    assert out["hpa"]["subcellular_main_location"] == ["nucleoplasm"]
    assert {"tissue": "Adrenal gland", "level": "High"} in out["hpa"]["tissues"]


def test_binary_not_installed_returns_none(monkeypatch):
    monkeypatch.setattr(bc, "_find_biomcp", lambda: None)
    assert bc.gene_pathways_hpa("BRD7") is None
    assert bc.is_available() is False


def test_unresolved_gene_returns_none_not_fabricated(monkeypatch):
    """Real observed behavior: exit code 1, empty stdout, error on stderr."""
    monkeypatch.setattr(bc, "_find_biomcp", lambda: "/fake/biomcp")
    monkeypatch.setattr(bc.subprocess, "run", _fake_run_factory(stdout="", returncode=1))
    assert bc.gene_pathways_hpa("NOTAREALGENE12345") is None


def test_timeout_returns_none(monkeypatch):
    def _raise_timeout(cmd, capture_output, text, timeout):
        raise subprocess.TimeoutExpired(cmd, timeout)
    monkeypatch.setattr(bc, "_find_biomcp", lambda: "/fake/biomcp")
    monkeypatch.setattr(bc.subprocess, "run", _raise_timeout)
    assert bc.gene_pathways_hpa("BRD7") is None


def test_malformed_json_returns_none(monkeypatch):
    monkeypatch.setattr(bc, "_find_biomcp", lambda: "/fake/biomcp")
    monkeypatch.setattr(bc.subprocess, "run", _fake_run_factory(stdout="not json{{{"))
    assert bc.gene_pathways_hpa("BRD7") is None


def test_no_hpa_data_is_none_not_missing_key(monkeypatch):
    """A gene with pathways but no HPA entry must yield hpa=None explicitly,
    not a KeyError or a silently-omitted field."""
    payload = json.dumps({"pathways": [{"source": "Reactome", "id": "R-1", "name": "x"}]})
    monkeypatch.setattr(bc, "_find_biomcp", lambda: "/fake/biomcp")
    monkeypatch.setattr(bc.subprocess, "run", _fake_run_factory(stdout=payload))
    out = bc.gene_pathways_hpa("SOMEGENE")
    assert out == {"pathways": [{"source": "Reactome", "id": "R-1", "name": "x"}], "hpa": None}


def test_empty_pathways_list_when_absent(monkeypatch):
    payload = json.dumps({"hpa": {"reliability": "Low"}})
    monkeypatch.setattr(bc, "_find_biomcp", lambda: "/fake/biomcp")
    monkeypatch.setattr(bc.subprocess, "run", _fake_run_factory(stdout=payload))
    out = bc.gene_pathways_hpa("SOMEGENE")
    assert out["pathways"] == []
    assert out["hpa"]["reliability"] == "Low"


def _biomcp_reachable() -> bool:
    if not bc.is_available():
        return False
    try:
        result = bc.gene_pathways_hpa("BRD7", timeout=15.0)
        return result is not None
    except Exception:  # noqa: BLE001
        return False


@pytest.mark.skipif(not _biomcp_reachable(), reason="biomcp not installed or unreachable")
def test_live_biomcp_call_matches_expected_shape():
    """Real end-to-end call against the actual installed biomcp binary, when
    available -- not just the mocked fixture above."""
    out = bc.gene_pathways_hpa("BRD7", timeout=30.0)
    assert out is not None
    assert isinstance(out["pathways"], list) and len(out["pathways"]) > 0
    assert all({"source", "id", "name"} <= set(p.keys()) for p in out["pathways"])
    if out["hpa"] is not None:
        assert "tissues" in out["hpa"]
