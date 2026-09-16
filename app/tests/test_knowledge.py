"""Tests for the knowledge module (Europe PMC / UniProt / ChEMBL).

Network-dependent tests skip gracefully when offline. Run: uv run pytest app/tests -q
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

APP = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(APP))

from agents import knowledge  # noqa: E402


def _online() -> bool:
    try:
        knowledge._get("https://rest.uniprot.org/uniprotkb/search?query=gene_exact:EGFR"
                       "+AND+organism_id:9606+AND+reviewed:true&format=json&fields=accession&size=1",
                       timeout=10)
        return True
    except Exception:
        return False


def test_uniprot_egfr():
    if not _online():
        pytest.skip("offline")
    up = knowledge.uniprot("EGFR")
    assert up and up["accession"] == "P00533"
    assert "kinase" in (up["protein_name"] + " ".join(up["keywords"])).lower() or up["domains"]
    assert up["n_pdb"] > 0


def test_chembl_egfr_has_drugs():
    if not _online():
        pytest.skip("offline")
    ch = knowledge.chembl_drugs("P00533")
    assert ch["target_chembl_id"] == "CHEMBL203"
    assert ch["n_mechanisms"] > 0
    assert len(ch["mechanisms"]) > 0 and ch["mechanisms"][0]["moa"]


def test_chembl_empty_for_nondruggable():
    if not _online():
        pytest.skip("offline")
    # LDB1 (Q86U70) is a transcription cofactor -> no approved-drug mechanism
    ch = knowledge.chembl_drugs("Q86U70")
    assert ch["n_mechanisms"] == 0


def test_literature_returns_papers():
    if not _online():
        pytest.skip("offline")
    lit = knowledge.literature("AURKA", "Neuroblastoma", 4)
    assert len(lit) >= 1
    assert all("title" in p and "pmid" in p for p in lit)


def test_dossier_shape_and_druggability():
    if not _online():
        pytest.skip("offline")
    d = knowledge.dossier("AURKA", "Neuroblastoma", 3)
    assert d["gene"] == "AURKA"
    assert d["uniprot"] and d["uniprot"]["accession"]
    assert d["chembl"]["n_mechanisms"] > 0
    assert "PDB" in d["druggability"] or "structure" in d["druggability"]


def test_chembl_empty_accession_safe():
    ch = knowledge.chembl_drugs("")
    assert ch["n_mechanisms"] == 0 and ch["mechanisms"] == []
