#!/usr/bin/env python3
"""
Drug-Target Scoring Pipeline (ported into SCRAP-AI)
====================================================
Original collaborator pipeline, adapted for in-app use: `run_pipeline()` is
called directly (no CLI), CACHE_DIR is redirected under this sandbox's data/
directory, and a `skip_plots` flag avoids writing matplotlib PNGs since the
UI renders its own interactive Plotly views from the returned DataFrame.
See user_data/drug_targets/README.md and run.md for full design rationale;
core scoring/API logic below is unchanged from that original.

Takes a list of gene symbols, queries ChEMBL, DGIdb, Open Targets, and
optionally PRISM/DepMap to find drugs targeting each gene's protein product,
annotates FDA approval status and cancer indications, assigns a 3-tier score,
and outputs a CSV results table plus two heatmaps (genes x drugs, genes x cancers).

Data sources (all free, no authentication):
    - MyGene.info    : gene symbol -> UniProt / Ensembl ID resolution
    - ChEMBL (EBI)   : drug-target bioactivity, clinical phase, drug indications
    - DGIdb v5       : aggregated drug-gene interactions (includes DrugBank)
    - Open Targets   : gene -> disease (cancer) associations
    - PRISM / DepMap : optional drug sensitivity across cancer cell lines

License: ChEMBL CC BY-SA 3.0 | DGIdb / Open Targets / MyGene.info / PRISM see their terms.
"""

import argparse
import json
import os
import sys
import time
import hashlib
import warnings
from pathlib import Path
from typing import Optional

import pandas as pd
import numpy as np
import requests

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
_SANDBOX = Path(__file__).resolve().parent.parent.parent
CACHE_DIR = _SANDBOX / "data" / "drug_targets_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

API_DELAY = 0.4  # seconds between API calls per source
REQUEST_TIMEOUT = 30
MAX_RETRIES = 3

# Cancer EFO/MONDO keyword filters for Open Targets disease filtering
CANCER_KEYWORDS = [
    "cancer", "carcinoma", "sarcoma", "leukemia", "leukaemia", "lymphoma",
    "melanoma", "glioma", "blastoma", "neoplasm", "tumor", "tumour",
    "myeloma", "mesothelioma", "adenocarcinoma", "oncolog",
]


# ---------------------------------------------------------------------------
# Utility: cached HTTP GET / POST with retry
# ---------------------------------------------------------------------------
def _cache_key(url, payload=None):
    raw = url + (json.dumps(payload, sort_keys=True) if payload else "")
    return hashlib.md5(raw.encode()).hexdigest()


def cached_get(url, params=None, cache_prefix="get"):
    """HTTP GET with disk cache and retry."""
    ck = CACHE_DIR / f"{cache_prefix}_{_cache_key(url, str(params))}.json"
    if ck.exists():
        return json.loads(ck.read_text())

    for attempt in range(MAX_RETRIES):
        try:
            r = requests.get(url, params=params, timeout=REQUEST_TIMEOUT,
                             headers={"Accept": "application/json"})
            r.raise_for_status()
            data = r.json()
            ck.write_text(json.dumps(data))
            time.sleep(API_DELAY)
            return data
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                time.sleep(2 ** attempt)
            else:
                print(f"  [WARN] GET {url} failed after {MAX_RETRIES} retries: {e}")
                return None


def cached_post(url, payload, cache_prefix="post"):
    """HTTP POST (GraphQL) with disk cache and retry."""
    ck = CACHE_DIR / f"{cache_prefix}_{_cache_key(url, payload)}.json"
    if ck.exists():
        return json.loads(ck.read_text())

    for attempt in range(MAX_RETRIES):
        try:
            r = requests.post(url, json=payload, timeout=REQUEST_TIMEOUT,
                              headers={"Content-Type": "application/json"})
            r.raise_for_status()
            data = r.json()
            ck.write_text(json.dumps(data))
            time.sleep(API_DELAY)
            return data
        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                time.sleep(2 ** attempt)
            else:
                print(f"  [WARN] POST {url} failed after {MAX_RETRIES} retries: {e}")
                return None


# ---------------------------------------------------------------------------
# Step 1: Gene Resolver (MyGene.info)
# ---------------------------------------------------------------------------
class GeneResolver:
    """Resolve gene symbols to UniProt accessions and Ensembl gene IDs."""

    BASE = "https://mygene.info/v3/query"

    def resolve(self, gene_symbol):
        """Return dict with uniprot, ensembl_id, name, or None if not found."""
        data = cached_get(
            self.BASE,
            params={
                "q": f"symbol:{gene_symbol}",
                "fields": "symbol,name,uniprot,ensembl.gene",
                "species": "human",
                "size": 1,
            },
            cache_prefix="mygene",
        )
        if not data or not data.get("hits"):
            # Try as alias / query string
            data = cached_get(
                self.BASE,
                params={
                    "q": gene_symbol,
                    "fields": "symbol,name,uniprot,ensembl.gene",
                    "species": "human",
                    "size": 1,
                },
                cache_prefix="mygene",
            )
        if not data or not data.get("hits"):
            return None

        hit = data["hits"][0]
        uniprot = hit.get("uniprot", {})
        if isinstance(uniprot, dict):
            # MyGene.info uses "Swiss-Prot" for the primary accession
            uniprot_acc = uniprot.get("Swiss-Prot")
        else:
            uniprot_acc = None

        ensembl = hit.get("ensembl")
        ensembl_id = None
        if isinstance(ensembl, list) and ensembl:
            ensembl_id = ensembl[0].get("gene")
        elif isinstance(ensembl, dict):
            ensembl_id = ensembl.get("gene")

        return {
            "symbol": hit.get("symbol", gene_symbol),
            "name": hit.get("name", ""),
            "uniprot": uniprot_acc,
            "ensembl_id": ensembl_id,
        }


# ---------------------------------------------------------------------------
# Step 2: ChEMBL Client
# ---------------------------------------------------------------------------
class ChEMBLClient:
    """Query ChEMBL for drug-target bioactivity and drug indications."""

    BASE = "https://www.ebi.ac.uk/chembl/api/data"

    def _get(self, endpoint, params=None):
        url = f"{self.BASE}/{endpoint}"
        p = {"format": "json"}
        if params:
            p.update(params)
        return cached_get(url, params=p, cache_prefix="chembl")

    def find_target_by_uniprot(self, uniprot_acc):
        """Find ChEMBL target ID by UniProt accession."""
        if not uniprot_acc:
            return None
        data = self._get("target", {"target_components__accession": uniprot_acc})
        if not data or not data.get("targets"):
            return None
        # Prefer SINGLE PROTEIN targets
        for t in data["targets"]:
            if t.get("target_type") == "SINGLE PROTEIN":
                return t
        return data["targets"][0]

    def find_target_by_name(self, gene_symbol):
        """Fallback: search ChEMBL targets by gene name / synonym."""
        data = self._get("target", {"pref_name__icontains": gene_symbol})
        if not data or not data.get("targets"):
            return None
        for t in data["targets"]:
            if t.get("target_type") == "SINGLE PROTEIN":
                return t
        return data["targets"][0] if data["targets"] else None

    def get_activities_for_target(self, target_chembl_id):
        """Fetch all bioactivity records for a target. Returns list of activity dicts."""
        activities = []
        offset = 0
        limit = 100
        while True:
            data = self._get("activity", {
                "target_chembl_id": target_chembl_id,
                "limit": limit,
                "offset": offset,
            })
            if not data:
                break
            acts = data.get("activities", [])
            activities.extend(acts)
            # ChEMBL's total_count is unreliable (often 0); use returned batch length
            if len(acts) < limit:
                break
            offset += limit
            # Safety cap at 5000 activities per target
            if offset >= 5000:
                break
        return activities

    def get_molecule(self, molecule_chembl_id):
        """Get molecule details including max_phase."""
        data = self._get(f"molecule/{molecule_chembl_id}")
        if not data:
            return None
        mols = data.get("molecules", [])
        return mols[0] if mols else None

    def get_molecules_batch(self, molecule_chembl_ids):
        """Batch fetch molecule details. Returns dict {chembl_id -> molecule_dict}."""
        results = {}
        ids = list(set(molecule_chembl_ids))
        batch_size = 50
        for i in range(0, len(ids), batch_size):
            batch = ids[i:i+batch_size]
            ids_str = ",".join(batch)
            data = self._get("molecule", {
                "molecule_chembl_id__in": ids_str,
                "limit": batch_size,
            })
            if not data:
                continue
            for mol in data.get("molecules", []):
                mid = mol.get("molecule_chembl_id")
                if mid:
                    results[mid] = mol
        return results

    def get_drug_indications(self, molecule_chembl_id):
        """Get FDA/clinical indications for a drug molecule."""
        data = self._get(f"drug_indication", {
            "molecule_chembl_id": molecule_chembl_id,
            "limit": 100,
        })
        if not data:
            return []
        return data.get("drug_indications", [])

    def query_gene(self, gene_info):
        """
        Full ChEMBL query for one gene.
        Returns list of drug dicts:
            {drug_name, chembl_id, max_phase, molecule_type, best_potency_nM,
             potency_type, indications, target_chembl_id}
        """
        target = self.find_target_by_uniprot(gene_info.get("uniprot"))
        if not target:
            target = self.find_target_by_name(gene_info["symbol"])
        if not target:
            return []

        target_id = target["target_chembl_id"]
        activities = self.get_activities_for_target(target_id)
        if not activities:
            return []

        # Group by molecule, keep best potency
        mol_dict = {}  # chembl_id -> {info}
        for act in activities:
            mid = act.get("molecule_chembl_id")
            if not mid:
                continue
            std_type = act.get("standard_type", "")
            std_val = act.get("standard_value")
            std_units = act.get("standard_units", "")

            # Only keep nM potency values
            potency_nM = None
            if std_val and std_units == "nM" and std_type in ("IC50", "Ki", "Kd", "Kd(app)", "EC50", "Potency"):
                try:
                    potency_nM = float(std_val)
                except (ValueError, TypeError):
                    pass

            if mid not in mol_dict:
                mol_dict[mid] = {
                    "drug_name": act.get("molecule_pref_name") or mid,
                    "chembl_id": mid,
                    "max_phase": -1,
                    "molecule_type": "",
                    "best_potency_nM": None,
                    "potency_type": "",
                    "indications": [],
                    "target_chembl_id": target_id,
                }
            mol = mol_dict[mid]
            if potency_nM is not None:
                if mol["best_potency_nM"] is None or potency_nM < mol["best_potency_nM"]:
                    mol["best_potency_nM"] = potency_nM
                    mol["potency_type"] = std_type

        # Fetch molecule details (max_phase, names) in batch
        mol_ids = list(mol_dict.keys())
        mol_details = self.get_molecules_batch(mol_ids)

        for mid in mol_ids:
            entry = mol_dict[mid]
            mol = mol_details.get(mid)
            if mol:
                entry["max_phase"] = int(float(mol.get("max_phase", 0) or 0))
                entry["molecule_type"] = mol.get("molecule_type", "")
                if not entry["drug_name"] or entry["drug_name"] == mid:
                    entry["drug_name"] = mol.get("pref_name") or mid

            # Get indications for approved/late-stage drugs only (to limit API calls)
            if entry["max_phase"] >= 3:
                inds = self.get_drug_indications(mid)
                entry["indications"] = [
                    {"term": ind.get("efo_term", "") or ind.get("mesh_heading", ""),
                     "phase": ind.get("max_phase_for_ind", ""),
                     "mesh": ind.get("mesh_heading", ""),
                     "efo_id": ind.get("efo_id", "")}
                    for ind in inds
                ]

        return list(mol_dict.values())

    def find_molecule_by_name(self, name):
        """Search ChEMBL molecule by preferred name. Returns chembl_id or None."""
        data = self._get("molecule", {"pref_name__iexact": name, "limit": 1})
        if not data:
            return None
        mols = data.get("molecules", [])
        if mols:
            return mols[0].get("molecule_chembl_id")
        return None

    def get_indications_for_drug(self, drug_name, chembl_id=None):
        """
        Get cancer indications for a drug, searching by ChEMBL ID or name.
        Returns list of indication dicts (same format as query_gene indications).
        """
        mid = chembl_id
        if not mid:
            mid = self.find_molecule_by_name(drug_name)
        if not mid:
            return []
        inds = self.get_drug_indications(mid)
        return [
            {"term": ind.get("efo_term", "") or ind.get("mesh_heading", ""),
             "phase": ind.get("max_phase_for_ind", ""),
             "mesh": ind.get("mesh_heading", ""),
             "efo_id": ind.get("efo_id", "")}
            for ind in inds
        ]


# ---------------------------------------------------------------------------
# Step 3: DGIdb Client
# ---------------------------------------------------------------------------
class DGIdbClient:
    """Query DGIdb v5 GraphQL API for drug-gene interactions."""

    URL = "https://dgidb.org/api/graphql"

    QUERY = """
    query GetInteractions($geneNames: [String!]!) {
      genes(names: $geneNames) {
        nodes {
          name
          interactions {
            drug {
              name
              conceptId
              approved
            }
            interactionScore
            interactionTypes {
              type
              directionality
            }
            publications {
              pmid
            }
            sources {
              sourceDbName
            }
          }
        }
      }
    }
    """

    def query_gene(self, gene_symbol):
        """Return list of drug dicts from DGIdb for one gene."""
        payload = {"query": self.QUERY, "variables": {"geneNames": [gene_symbol]}}
        data = cached_post(self.URL, payload, cache_prefix="dgidb")
        if not data or not data.get("data"):
            return []

        genes = data["data"].get("genes", {})
        nodes = genes.get("nodes", []) if genes else []
        if not nodes:
            return []

        interactions = nodes[0].get("interactions", []) or []
        results = []
        for inter in interactions:
            drug = inter.get("drug", {}) or {}
            itypes = inter.get("interactionTypes", []) or []
            type_str = "; ".join(t.get("type", "") for t in itypes if t.get("type"))
            sources = inter.get("sources", []) or []
            source_str = "; ".join(s.get("sourceDbName", "") for s in sources if s.get("sourceDbName"))

            results.append({
                "drug_name": drug.get("name", ""),
                "concept_id": drug.get("conceptId", ""),
                "approved": drug.get("approved", False),
                "interaction_score": inter.get("interactionScore"),
                "interaction_type": type_str,
                "sources": source_str,
            })
        return results


# ---------------------------------------------------------------------------
# Step 4: Open Targets Client
# ---------------------------------------------------------------------------
class OpenTargetsClient:
    """Query Open Targets GraphQL API for gene-disease (cancer) associations."""

    URL = "https://api.platform.opentargets.org/api/v4/graphql"

    # Known-drugs-for-target query: replaces the old per-drug ChEMBL
    # molecule-by-name reverse lookup (which hit a flaky/500-erroring ChEMBL
    # endpoint). ONE call per gene returns every known drug for that target,
    # each with its disease indications + trial phase, already keyed by gene
    # (no fragile name-matching needed).
    DRUGS_QUERY = """
    query TargetDrugs($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        drugAndClinicalCandidates {
          rows {
            maxClinicalStage
            drug {
              id
              name
              drugType
              indications {
                rows {
                  disease { id name }
                  maxClinicalStage
                }
              }
            }
          }
        }
      }
    }
    """

    _STAGE_TO_PHASE = {
        "APPROVAL": 4.0, "PHASE_4": 4.0, "PHASE_3": 3.0, "PHASE_2_3": 2.5,
        "PHASE_2": 2.0, "PHASE_1_2": 1.5, "PHASE_1": 1.0, "EARLY_PHASE_1": 0.5,
    }

    @classmethod
    def _stage_to_phase(cls, stage):
        if not stage:
            return None
        return cls._STAGE_TO_PHASE.get(str(stage).upper())

    def query_known_drugs(self, ensembl_id):
        """Return {drug_name_lower: {max_stage, indications: [{term,phase,mesh,efo_id}]}}
        for every known drug against this target. One cached call per gene.
        """
        if not ensembl_id:
            return {}
        payload = {"query": self.DRUGS_QUERY, "variables": {"ensemblId": ensembl_id}}
        data = cached_post(self.URL, payload, cache_prefix="ot_known_drugs")
        if not data or not data.get("data"):
            return {}
        target = data["data"].get("target")
        if not target:
            return {}
        rows = (target.get("drugAndClinicalCandidates") or {}).get("rows", []) or []
        out = {}
        for row in rows:
            d = row.get("drug") or {}
            name = (d.get("name") or "").strip()
            if not name:
                continue
            inds = []
            for ind in ((d.get("indications") or {}).get("rows", []) or []):
                disease = ind.get("disease") or {}
                inds.append({
                    "term": disease.get("name", ""),
                    "phase": self._stage_to_phase(ind.get("maxClinicalStage")),
                    "mesh": "",
                    "efo_id": disease.get("id", ""),
                })
            out[name.lower()] = {
                "max_stage": self._stage_to_phase(row.get("maxClinicalStage")),
                "indications": inds,
            }
        return out

    QUERY = """
    query AssociatedDiseases($ensemblId: String!, $size: Int!) {
      disease(ensemblId: $ensemblId) {
        id
        associatedDiseases(ensemblId: $ensemblId, BSize: $size) {
          rows {
            disease {
              id
              name
              therapeuticAreas {
                id
                name
              }
            }
            score
            datatypeScores {
              id
              score
            }
          }
        }
      }
    }
    """

    # Open Targets v4 API: associatedDiseases takes Bs (list) and page (Pagination)
    QUERY_V4 = """
    query AssociatedDiseases($ensemblId: String!, $size: Int!) {
      target(ensemblId: $ensemblId) {
        associatedDiseases(Bs: [], page: {size: $size, index: 0}) {
          rows {
            disease {
              id
              name
              therapeuticAreas {
                id
                name
              }
            }
            score
            datatypeScores {
              id
              score
            }
          }
        }
      }
    }
    """

    def _is_cancer(self, disease_name, therapeutic_areas):
        """Check if a disease is a cancer type."""
        name_lower = (disease_name or "").lower()
        for kw in CANCER_KEYWORDS:
            if kw in name_lower:
                return True
        for ta in (therapeutic_areas or []):
            ta_name = (ta.get("name", "") or "").lower()
            ta_id = (ta.get("id", "") or "")
            # EFO:0002618 = cancer or "neoplasm" therapeutic area
            if "cancer" in ta_name or "neoplasm" in ta_name or ta_id == "EFO_0002618" or ta_id == "OT_0000134":
                return True
        return False

    def query_gene(self, ensembl_id):
        """Return list of cancer associations: {cancer_name, score, disease_id}."""
        if not ensembl_id:
            return []
        payload = {
            "query": self.QUERY_V4,
            "variables": {"ensemblId": ensembl_id, "size": 250},
        }
        data = cached_post(self.URL, payload, cache_prefix="opentargets")
        if not data or not data.get("data"):
            return []

        target = data["data"].get("target")
        if not target:
            return []
        assoc = target.get("associatedDiseases", {})
        rows = assoc.get("rows", []) if assoc else []
        if not rows:
            return []

        cancers = []
        for row in rows:
            disease = row.get("disease", {}) or {}
            name = disease.get("name", "")
            areas = disease.get("therapeuticAreas", []) or []
            if self._is_cancer(name, areas):
                cancers.append({
                    "cancer_name": name,
                    "disease_id": disease.get("id", ""),
                    "score": row.get("score", 0),
                })
        return cancers


# ---------------------------------------------------------------------------
# Step 4b: Cancer Type Resolver (Open Targets disease search + ontology)
# ---------------------------------------------------------------------------
class CancerResolver:
    """Resolve a cancer type name to a set of matching terms (including subtypes)."""

    URL = "https://api.platform.opentargets.org/api/v4/graphql"

    SEARCH_QUERY = """
    query SearchDisease($queryString: String!) {
      search(queryString: $queryString, entityNames: ["disease"]) {
        hits {
          id
          name
          entity
        }
      }
    }
    """

    CHILDREN_QUERY = """
    query DiseaseChildren($efoId: String!) {
      disease(efoId: $efoId) {
        id
        name
        children {
          id
          name
        }
      }
    }
    """

    def resolve(self, cancer_name):
        """
        Resolve a cancer type name to a set of terms for matching.
        Returns dict: {disease_id, disease_name, cancer_terms: set of lowercased names}
        """
        cancer_name_lower = cancer_name.lower().strip()

        # Step 1: Search for the disease in Open Targets
        payload = {"query": self.SEARCH_QUERY, "variables": {"queryString": cancer_name}}
        data = cached_post(self.URL, payload, cache_prefix="ot_cancer_search")

        disease_id = None
        disease_name = cancer_name  # fallback to user input
        cancer_terms = {cancer_name_lower}

        if data and data.get("data", {}).get("search"):
            hits = data["data"]["search"].get("hits", [])
            if hits:
                # Pick the first disease hit
                disease_id = hits[0].get("id")
                disease_name = hits[0].get("name", cancer_name)
                cancer_terms.add(disease_name.lower())

        # Step 2: Fetch children (subtypes) if we found a disease ID
        if disease_id:
            payload2 = {"query": self.CHILDREN_QUERY, "variables": {"efoId": disease_id}}
            data2 = cached_post(self.URL, payload2, cache_prefix="ot_cancer_children")
            if data2 and data2.get("data", {}).get("disease"):
                children = data2["data"]["disease"].get("children", []) or []
                for child in children:
                    child_name = child.get("name", "")
                    if child_name:
                        cancer_terms.add(child_name.lower())
                print(f"  Cancer resolved: {disease_name} ({disease_id}) "
                      f"with {len(children)} subtypes")
            else:
                print(f"  Cancer resolved: {disease_name} ({disease_id}), no subtypes found")
        else:
            print(f"  Cancer type '{cancer_name}' not found in Open Targets. "
                  f"Using substring matching only.")

        # Also add the core keyword from the cancer name for substring matching
        # e.g., "neuroblastoma" from "Neuroblastoma", "lung" from "lung cancer"
        core_words = [w for w in cancer_name_lower.split()
                      if w not in {"cancer", "carcinoma", "tumor", "tumour",
                                   "neoplasm", "the", "of", "type"}]
        if core_words:
            cancer_terms.update(core_words)

        return {
            "disease_id": disease_id,
            "disease_name": disease_name,
            "cancer_terms": cancer_terms,
        }

    @staticmethod
    def matches_query_cancer(indications, cancer_terms):
        """
        Check if any drug indication matches the user's cancer type.
        Returns (matched: bool, max_phase_for_ind: float or None, matched_terms: list)
        """
        matched_terms = []
        max_phase = None
        for ind in indications:
            term = (ind.get("term", "") or "").lower()
            mesh = (ind.get("mesh", "") or "").lower()
            text = f"{term} {mesh}"
            for ct in cancer_terms:
                if ct in text:
                    matched_terms.append(ind.get("term", ""))
                    phase_str = ind.get("phase", "")
                    try:
                        phase = float(phase_str) if phase_str else None
                        if phase is not None:
                            if max_phase is None or phase > max_phase:
                                max_phase = phase
                    except (ValueError, TypeError):
                        pass
                    break
        return (len(matched_terms) > 0, max_phase, matched_terms)

    @staticmethod
    def get_cancer_indications(indications):
        """Return all cancer-related indications (for 'any cancer' check)."""
        results = []
        for ind in indications:
            term = (ind.get("term", "") or "").lower()
            if any(kw in term for kw in CANCER_KEYWORDS):
                results.append(ind.get("term", ""))
        return list(set(results))

    @staticmethod
    def get_non_cancer_indications(indications):
        """Return non-cancer indications."""
        results = []
        for ind in indications:
            term = (ind.get("term", "") or "").lower()
            if not any(kw in term for kw in CANCER_KEYWORDS):
                results.append(ind.get("term", ""))
        return list(set(results))


# ---------------------------------------------------------------------------
# Step 5: PRISM Client (optional)
# ---------------------------------------------------------------------------
class PRISMClient:
    """Optional: download and query PRISM Repurposing drug sensitivity data."""

    # Figshare download URLs for PRISM Repurposing 24Q2
    # We use the DepMap Figshare release
    DATA_MATRIX_URL = "https://figshare.com/ndownloader/files/47896686"  # Primary data matrix
    COMPOUND_META_URL = "https://figshare.com/ndownloader/files/47896680"  # Compound metadata
    CELL_LINE_META_URL = "https://figshare.com/ndownloader/files/47896677"  # Cell line metadata

    # Fallback: try direct DepMap portal filenames
    FALLBACK_MATRIX = "https://depmap.org/portal/download/allfiles/?release=PRISM+Repurposing+Public+24Q2&file=Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv"

    def __init__(self, cache_dir=None):
        self.cache_dir = Path(cache_dir or CACHE_DIR)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.data_matrix = None
        self.compound_map = None  # name -> compound_id
        self.cell_line_map = None  # cell_line_id -> lineage
        self._loaded = False

    def _download_file(self, url, dest):
        """Download a file with progress."""
        if dest.exists():
            return True
        try:
            print(f"  [PRISM] Downloading {dest.name}...")
            r = requests.get(url, timeout=300, stream=True)
            r.raise_for_status()
            with open(dest, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    f.write(chunk)
            print(f"  [PRISM] Downloaded {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
            return True
        except Exception as e:
            print(f"  [PRISM] Download failed for {dest.name}: {e}")
            return False

    def load(self):
        """Download (if needed) and load PRISM data into memory."""
        if self._loaded:
            return True

        matrix_path = self.cache_dir / "prism_primary_matrix.csv"
        compound_path = self.cache_dir / "prism_compound_meta.csv"
        cellline_path = self.cache_dir / "prism_cellline_meta.csv"

        # Try downloading from Figshare
        ok = True
        if not matrix_path.exists():
            ok = self._download_file(self.DATA_MATRIX_URL, matrix_path) or \
                 self._download_file(self.FALLBACK_MATRIX, matrix_path)
        if not compound_path.exists():
            ok = ok and self._download_file(self.COMPOUND_META_URL, compound_path)
        if not cellline_path.exists():
            ok = ok and self._download_file(self.CELL_LINE_META_URL, cellline_path)

        if not ok or not matrix_path.exists():
            print("  [PRISM] Could not load PRISM data. Skipping PRISM enrichment.")
            return False

        try:
            # Load compound metadata -> name to ID mapping
            comp_df = pd.read_csv(compound_path, low_memory=False)
            self.compound_map = {}
            name_col = "name" if "name" in comp_df.columns else comp_df.columns[1]
            id_col = "compound_id" if "compound_id" in comp_df.columns else comp_df.columns[0]
            for _, row in comp_df.iterrows():
                name = str(row[name_col]).strip().lower()
                cid = str(row[id_col])
                if name and name != "nan":
                    self.compound_map[name] = cid

            # Load cell line metadata -> ID to lineage mapping
            cl_df = pd.read_csv(cellline_path, low_memory=False)
            self.cell_line_map = {}
            cl_id_col = "cell_line_id" if "cell_line_id" in cl_df.columns else cl_df.columns[0]
            lineage_col = None
            for c in cl_df.columns:
                if "lineage" in c.lower():
                    lineage_col = c
                    break
            if lineage_col:
                for _, row in cl_df.iterrows():
                    self.cell_line_map[str(row[cl_id_col])] = str(row[lineage_col])

            # Load data matrix
            self.data_matrix = pd.read_csv(matrix_path, low_memory=False, index_col=0)
            self._loaded = True
            print(f"  [PRISM] Loaded: {len(self.data_matrix)} compounds x {len(self.data_matrix.columns)} cell lines")
            return True
        except Exception as e:
            print(f"  [PRISM] Error loading data: {e}")
            return False

    def get_sensitivity(self, drug_name):
        """Return sensitivity stats for a drug, or None."""
        if not self._loaded:
            return None
        name_lower = drug_name.strip().lower()
        cid = self.compound_map.get(name_lower)
        if not cid:
            # Try partial match
            for k, v in self.compound_map.items():
                if name_lower in k or k in name_lower:
                    cid = v
                    break
        if not cid:
            return None

        # Find row in data matrix
        if cid in self.data_matrix.index:
            row = self.data_matrix.loc[cid]
        else:
            # Try by name in index
            if name_lower in [str(x).lower() for x in self.data_matrix.index]:
                idx = [str(x).lower() for x in self.data_matrix.index].index(name_lower)
                row = self.data_matrix.iloc[idx]
            else:
                return None

        # AUC values: lower = more sensitive (PRISM uses log2 fold change)
        # Sensitivity threshold: AUC < 0 (negative = sensitive)
        values = pd.to_numeric(row, errors="coerce").dropna()
        if values.empty:
            return None

        sensitive = values[values < -0.5]  # threshold for sensitivity
        pct_sensitive = (len(sensitive) / len(values)) * 100

        # Top lineages
        lineage_sens = {}
        for cl_id, val in values.items():
            lineage = self.cell_line_map.get(str(cl_id), "Unknown")
            if val < -0.5:
                lineage_sens[lineage] = lineage_sens.get(lineage, 0) + 1

        top_lineages = sorted(lineage_sens.items(), key=lambda x: -x[1])[:5]
        top_lineage_str = "; ".join(f"{l}({n})" for l, n in top_lineages)

        return {
            "pct_sensitive": round(pct_sensitive, 1),
            "mean_auc": round(values.mean(), 3),
            "n_cell_lines": len(values),
            "top_sensitive_lineages": top_lineage_str,
        }


# ---------------------------------------------------------------------------
# Step 6: Scorer
# ---------------------------------------------------------------------------
class Scorer:
    """Scoring for gene-drug pairs (3-tier original and 5-tier cancer-specific)."""

    SCORE_LABELS_5TIER = {
        5: "Approved for your cancer",
        4: "Clinical trial for your cancer",
        3: "Approved for other cancers",
        2: "Approved (non-cancer) or trial (any cancer)",
        1: "Experimental",
    }

    @staticmethod
    def score(drug_entry, gene_cancers):
        """
        Assign score 1-3 based on FDA approval and cancer indication (original mode).

        Score 3: FDA-approved + cancer indication
        Score 2: FDA-approved (any disease) OR clinical trial + cancer evidence
        Score 1: experimental/preclinical, targets the gene
        """
        max_phase = drug_entry.get("max_phase", 0) or 0
        fda_approved = drug_entry.get("fda_approved", False)
        has_cancer_indication = drug_entry.get("cancer_indication", False)

        # Also check if gene itself is associated with cancer
        gene_has_cancer = len(gene_cancers) > 0

        # FDA approved = max_phase 4 or DGIdb approved flag
        is_fda = (max_phase >= 4) or fda_approved

        # Clinical trial = max_phase 1-3
        in_clinical = 1 <= max_phase <= 3

        if is_fda and has_cancer_indication:
            return 3
        elif is_fda and gene_has_cancer:
            # FDA-approved drug targeting a gene associated with cancer
            return 3
        elif is_fda:
            return 2
        elif in_clinical and (has_cancer_indication or gene_has_cancer):
            return 2
        else:
            return 1

    @staticmethod
    def score_cancer_specific(drug_entry, cancer_resolver, cancer_info):
        """
        Assign score 1-5 based on FDA approval and cancer-type-specific indication.

        Score 5: FDA-approved AND indicated for the user's cancer type
        Score 4: Clinical trial AND indicated for the user's cancer type
        Score 3: FDA-approved with cancer indications, but NOT for the user's cancer
        Score 2: FDA-approved (non-cancer) OR clinical trial for any cancer (not user's)
        Score 1: Experimental, no clinical cancer evidence
        """
        max_phase = drug_entry.get("max_phase", 0) or 0
        if max_phase < 0:
            max_phase = 0
        fda_approved = drug_entry.get("fda_approved", False)
        indications = drug_entry.get("indications", [])
        has_chembl_indications = len(indications) > 0

        is_fda = (max_phase >= 4) or fda_approved
        in_clinical = 1 <= max_phase <= 3

        # Check if drug is indicated for the user's cancer
        indicated_for_query = False
        query_cancer_phase = None
        other_cancer_inds = []

        if has_chembl_indications:
            indicated_for_query, query_cancer_phase, matched = \
                cancer_resolver.matches_query_cancer(indications, cancer_info["cancer_terms"])
            other_cancer_inds = cancer_resolver.get_cancer_indications(indications)
            # Remove query cancer matches from "other" list
            other_cancer_inds = [c for c in other_cancer_inds
                                 if c not in matched]
        else:
            # No ChEMBL indications available (DGIdb-only drug)
            # Can't determine cancer-specificity
            if is_fda:
                return 2, Scorer.SCORE_LABELS_5TIER[2], False, None, []
            else:
                return 1, Scorer.SCORE_LABELS_5TIER[1], False, None, []

        # Determine score
        if indicated_for_query:
            # Drug is indicated for the user's cancer
            # Use the indication-specific phase if available, else overall max_phase
            effective_phase = query_cancer_phase if query_cancer_phase else max_phase
            if is_fda or effective_phase >= 4:
                score = 5
            elif in_clinical or (1 <= effective_phase <= 3):
                score = 4
            else:
                # Indicated for the cancer but phase unknown/0 — still better than nothing
                score = 4 if is_fda else 3
        elif is_fda and other_cancer_inds:
            # FDA-approved for other cancers, not the user's
            score = 3
        elif is_fda:
            # FDA-approved for non-cancer disease
            score = 2
        elif in_clinical and other_cancer_inds:
            # Clinical trial for other cancers
            score = 2
        else:
            # Experimental
            score = 1

        label = Scorer.SCORE_LABELS_5TIER.get(score, "Unknown")
        return score, label, indicated_for_query, query_cancer_phase, other_cancer_inds


# ---------------------------------------------------------------------------
# Step 7: Heatmap Generator
# ---------------------------------------------------------------------------
class HeatmapGenerator:
    """Generate genes x drugs and genes x cancers heatmaps."""

    def __init__(self, output_dir):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def _setup_matplotlib(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        matplotlib.rcParams["font.family"] = ["Liberation Sans", "Arimo", "DejaVu Sans"]
        matplotlib.rcParams["svg.fonttype"] = "none"
        return plt

    def genes_drugs_heatmap(self, df, max_drugs=50, vmax=3):
        """Heatmap: rows=genes, columns=drugs, color=score."""
        plt = self._setup_matplotlib()
        import seaborn as sns

        if df.empty:
            print("  [WARN] No data for genes x drugs heatmap")
            return None

        # Pivot: gene x drug -> score
        pivot = df.pivot_table(index="Gene", columns="Drug_Name", values="Score",
                               aggfunc="max", fill_value=0)

        # For large drug lists, keep top-scoring drugs
        if pivot.shape[1] > max_drugs:
            # Keep drugs with highest total score, plus all score-3 drugs
            col_scores = pivot.sum(axis=0).sort_values(ascending=False)
            top_cols = list(col_scores.head(max_drugs).index)
            # Ensure all score-3 drugs are included
            score3_drugs = pivot.columns[(pivot == 3).any()].tolist()
            keep = list(set(top_cols + score3_drugs))[:max_drugs]
            pivot = pivot[keep]

        if pivot.empty or pivot.shape[1] == 0:
            print("  [WARN] No drug data after filtering for heatmap")
            return None

        # Sort genes by total score (descending)
        pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=False).index]

        fig, ax = plt.subplots(figsize=(max(8, pivot.shape[1] * 0.35 + 2),
                                        max(6, pivot.shape[0] * 0.4 + 2)))
        cmap = sns.color_palette("YlOrRd", as_cmap=True)
        sns.heatmap(pivot, annot=True, fmt="d", cmap=cmap, linewidths=0.5,
                    linecolor="white", cbar_kws={"label": "Score"},
                    vmin=0, vmax=vmax, ax=ax)
        ax.set_title("Gene-Drug Targeting Scores", fontsize=14, pad=12)
        ax.set_xlabel("Drug", fontsize=11)
        ax.set_ylabel("Gene", fontsize=11)
        plt.xticks(rotation=45, ha="right", fontsize=8)
        plt.yticks(rotation=0, fontsize=9)
        plt.tight_layout()

        png_path = self.output_dir / "heatmap_genes_drugs.png"
        svg_path = self.output_dir / "heatmap_genes_drugs.svg"
        fig.savefig(png_path, dpi=150, bbox_inches="tight")
        fig.savefig(svg_path, bbox_inches="tight")
        plt.close(fig)
        return png_path

    def genes_cancers_heatmap(self, df):
        """Heatmap: rows=genes, columns=cancer types, color=best score."""
        plt = self._setup_matplotlib()
        import seaborn as sns

        if df.empty:
            print("  [WARN] No data for genes x cancers heatmap")
            return None

        # Explode Cancer_Types into separate rows, then pivot
        records = []
        for _, row in df.iterrows():
            cancers = str(row.get("Cancer_Types", "")).split(";")
            cancers = [c.strip() for c in cancers if c.strip() and c.strip() != "nan"]
            if not cancers:
                # Use gene-level cancer associations from Open Targets
                cancers = str(row.get("Gene_Cancer_Associations", "")).split(";")
                cancers = [c.strip() for c in cancers if c.strip() and c.strip() != "nan"]
            for cancer in cancers:
                records.append({
                    "Gene": row["Gene"],
                    "Cancer": cancer,
                    "Score": row["Score"],
                })

        if not records:
            print("  [WARN] No cancer data for heatmap")
            return None

        cancer_df = pd.DataFrame(records)
        pivot = cancer_df.pivot_table(index="Gene", columns="Cancer",
                                      values="Score", aggfunc="max", fill_value=0)

        # Limit cancer columns if too many
        max_cancers = 40
        if pivot.shape[1] > max_cancers:
            col_scores = pivot.sum(axis=0).sort_values(ascending=False)
            pivot = pivot[col_scores.head(max_cancers).index]

        # Sort genes by total score
        pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=False).index]

        fig, ax = plt.subplots(figsize=(max(10, pivot.shape[1] * 0.4 + 2),
                                        max(6, pivot.shape[0] * 0.4 + 2)))
        cmap = sns.color_palette("YlOrRd", as_cmap=True)
        sns.heatmap(pivot, annot=True, fmt="d", cmap=cmap, linewidths=0.5,
                    linecolor="white", cbar_kws={"label": "Best Score"},
                    vmin=0, vmax=3, ax=ax)
        ax.set_title("Gene-Cancer Coverage Scores", fontsize=14, pad=12)
        ax.set_xlabel("Cancer Type", fontsize=11)
        ax.set_ylabel("Gene", fontsize=11)
        plt.xticks(rotation=45, ha="right", fontsize=8)
        plt.yticks(rotation=0, fontsize=9)
        plt.tight_layout()

        png_path = self.output_dir / "heatmap_genes_cancers.png"
        svg_path = self.output_dir / "heatmap_genes_cancers.svg"
        fig.savefig(png_path, dpi=150, bbox_inches="tight")
        fig.savefig(svg_path, bbox_inches="tight")
        plt.close(fig)
        return png_path

    def gene_scores_barchart(self, df, cancer_name, vmax=5):
        """Stacked bar chart: X=genes, Y=drug count, stacked by score tier."""
        plt = self._setup_matplotlib()

        if df.empty:
            print("  [WARN] No data for gene scores bar chart")
            return None

        # Count drugs per gene per score tier
        score_counts = df.groupby(["Gene", "Score"]).size().unstack(fill_value=0)

        # Ensure all score tiers 1..vmax are present as columns
        for s in range(1, vmax + 1):
            if s not in score_counts.columns:
                score_counts[s] = 0
        score_counts = score_counts[sorted(score_counts.columns, reverse=True)]

        # Sort genes by total drug count (descending)
        score_counts = score_counts.loc[score_counts.sum(axis=1).sort_values(ascending=False).index]

        # Color map: 5=dark red, 4=orange, 3=yellow, 2=light blue, 1=gray
        tier_colors = {
            5: "#B22222",   # firebrick (dark red)
            4: "#FF8C00",   # dark orange
            3: "#FFD700",   # gold (yellow)
            2: "#87CEEB",   # sky blue (light blue)
            1: "#A9A9A9",   # dark gray
        }
        colors = [tier_colors.get(s, "#A9A9A9") for s in score_counts.columns]

        fig, ax = plt.subplots(figsize=(max(8, len(score_counts) * 1.2 + 2), 6))
        x = range(len(score_counts))
        bottom = [0] * len(score_counts)
        for col_idx, score_tier in enumerate(score_counts.columns):
            vals = score_counts[score_tier].values
            label = Scorer.SCORE_LABELS_5TIER.get(score_tier, f"Score {score_tier}")
            ax.bar(x, vals, bottom=bottom, color=colors[col_idx],
                   label=f"{score_tier}: {label}", edgecolor="white", linewidth=0.5)
            bottom = [b + v for b, v in zip(bottom, vals)]

        ax.set_xticks(list(x))
        ax.set_xticklabels(score_counts.index, rotation=45, ha="right", fontsize=10)
        ax.set_ylabel("Number of Drugs", fontsize=11)
        ax.set_title(f"Drug Score Distribution for {cancer_name}", fontsize=14, pad=12)
        ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
        plt.tight_layout()

        png_path = self.output_dir / "barchart_gene_scores.png"
        svg_path = self.output_dir / "barchart_gene_scores.svg"
        fig.savefig(png_path, dpi=150, bbox_inches="tight")
        fig.savefig(svg_path, bbox_inches="tight")
        plt.close(fig)
        return png_path


# ---------------------------------------------------------------------------
# Main Pipeline
# ---------------------------------------------------------------------------
def run_pipeline(genes, output_dir, use_prism=False, min_score=2,
                 max_drugs_per_gene=30, cancer_type=None, skip_plots=False,
                 progress_cb=None):
    """Run the full drug-target scoring pipeline.

    If cancer_type is provided, uses 5-tier cancer-specific scoring.
    Otherwise, uses the original 3-tier cancer-agnostic scoring.

    skip_plots: if True, skip writing matplotlib PNGs (the SCRAP-AI UI renders
      its own interactive Plotly views from the returned DataFrame instead).
    progress_cb: optional callable(i, total, gene) invoked before each gene is
      queried, for UI progress reporting.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    resolver = GeneResolver()
    chembl = ChEMBLClient()
    dgidb = DGIdbClient()
    ot = OpenTargetsClient()
    scorer = Scorer()

    # Cancer-specific mode setup
    cancer_resolver = None
    cancer_info = None
    if cancer_type:
        cancer_resolver = CancerResolver()
        print(f"\nResolving cancer type: {cancer_type}...")
        cancer_info = cancer_resolver.resolve(cancer_type)
        print(f"  Cancer terms for matching: {len(cancer_info['cancer_terms'])} terms")
        terms_preview = ", ".join(sorted(cancer_info["cancer_terms"])[:10])
        print(f"  (includes subtypes: {terms_preview}...)\n")

    prism = None
    if use_prism:
        prism = PRISMClient()
        prism.load()

    all_results = []
    unresolved = []

    total = len(genes)
    for i, gene in enumerate(genes, 1):
        gene = gene.strip()
        if not gene:
            continue
        print(f"[{i}/{total}] Querying {gene}...")
        if progress_cb:
            try:
                progress_cb(i, total, gene)
            except Exception:
                pass

        # Step 1: Resolve gene
        gene_info = resolver.resolve(gene)
        if not gene_info:
            print(f"  [WARN] Could not resolve gene '{gene}'. Skipping.")
            unresolved.append(gene)
            continue

        print(f"  Resolved: {gene_info['symbol']} | UniProt: {gene_info.get('uniprot', 'N/A')} | "
              f"Ensembl: {gene_info.get('ensembl_id', 'N/A')}")

        # Step 4: Open Targets cancer associations
        # In cancer-specific mode, skip gene-disease queries (user already has this evidence)
        gene_cancers = []
        cancer_names = []
        if not cancer_type:
            print(f"  Querying Open Targets for cancer associations...")
            gene_cancers = ot.query_gene(gene_info.get("ensembl_id"))
            cancer_names = [c["cancer_name"] for c in gene_cancers]
            print(f"  Found {len(gene_cancers)} cancer associations")

        # Step 4b: Open Targets known-drugs-for-target (cancer-specific mode only).
        # ONE call covers every drug's disease indications + phase for this gene,
        # replacing the old per-drug ChEMBL molecule-by-name reverse lookup.
        ot_known_drugs = {}
        if cancer_type:
            print(f"  Querying Open Targets for known drugs + indications...")
            ot_known_drugs = ot.query_known_drugs(gene_info.get("ensembl_id"))
            print(f"  Found {len(ot_known_drugs)} known drugs with indication data")

        # Step 2: ChEMBL
        print(f"  Querying ChEMBL for drug-target bioactivity...")
        chembl_drugs = chembl.query_gene(gene_info)
        print(f"  Found {len(chembl_drugs)} compounds in ChEMBL")

        # Step 3: DGIdb
        print(f"  Querying DGIdb for drug-gene interactions...")
        dgidb_drugs = dgidb.query_gene(gene_info["symbol"])
        print(f"  Found {len(dgidb_drugs)} interactions in DGIdb")

        # Merge ChEMBL and DGIdb results
        # Build a merged dict keyed by lowercased drug name
        merged = {}

        for d in chembl_drugs:
            key = d["drug_name"].strip().lower()
            if key not in merged:
                merged[key] = {
                    "drug_name": d["drug_name"],
                    "chembl_id": d["chembl_id"],
                    "max_phase": d["max_phase"],
                    "molecule_type": d.get("molecule_type", ""),
                    "best_potency_nM": d.get("best_potency_nM"),
                    "potency_type": d.get("potency_type", ""),
                    "indications": d.get("indications", []),
                    "interaction_type": "",
                    "dgidb_approved": False,
                    "dgidb_sources": "",
                    "sources": ["ChEMBL"],
                }
            else:
                # Update with potentially better info
                entry = merged[key]
                if d["max_phase"] > entry["max_phase"]:
                    entry["max_phase"] = d["max_phase"]
                if d.get("best_potency_nM") and (entry["best_potency_nM"] is None or
                                                  d["best_potency_nM"] < entry["best_potency_nM"]):
                    entry["best_potency_nM"] = d["best_potency_nM"]
                    entry["potency_type"] = d.get("potency_type", "")

        for d in dgidb_drugs:
            key = d["drug_name"].strip().lower()
            if key not in merged:
                merged[key] = {
                    "drug_name": d["drug_name"],
                    "chembl_id": "",
                    "max_phase": -1,
                    "molecule_type": "",
                    "best_potency_nM": None,
                    "potency_type": "",
                    "indications": [],
                    "interaction_type": d.get("interaction_type", ""),
                    "dgidb_approved": d.get("approved", False),
                    "dgidb_sources": d.get("sources", ""),
                    "sources": ["DGIdb"],
                }
            else:
                entry = merged[key]
                entry["interaction_type"] = entry["interaction_type"] or d.get("interaction_type", "")
                entry["dgidb_approved"] = d.get("approved", False)
                entry["dgidb_sources"] = d.get("sources", "")
                if "DGIdb" not in entry["sources"]:
                    entry["sources"].append("DGIdb")

        if not merged:
            print(f"  [WARN] No drugs found for gene '{gene}'.")
            unresolved.append(f"{gene} (no drugs found)")
            continue

        # Score each drug
        for key, entry in merged.items():
            # Determine FDA approved
            fda = (entry["max_phase"] >= 4) or entry["dgidb_approved"]

            # Determine cancer indication from ChEMBL drug-specific indications
            chembl_cancer_inds = []
            all_diseases = []
            for ind in entry.get("indications", []):
                term = ind.get("term", "")
                all_diseases.append(term)
                term_lower = term.lower()
                if any(kw in term_lower for kw in CANCER_KEYWORDS):
                    chembl_cancer_inds.append(term)

            # Filter out overly generic terms
            GENERIC_TERMS = {"neoplasm", "cancer", "neoplasms", "tumor", "tumour",
                             "malignant neoplasm", "carcinoma", "neoplastic disease"}
            specific_cancer_inds = [c for c in set(chembl_cancer_inds)
                                    if c.lower() not in GENERIC_TERMS]

            # PRISM enrichment
            prism_sens = None
            if prism and prism._loaded:
                prism_sens = prism.get_sensitivity(entry["drug_name"])

            if cancer_type:
                # ---- 5-tier cancer-specific scoring ----
                # For DGIdb-only drugs without ChEMBL IDs, look up indications from the
                # pre-fetched Open Targets known-drugs-for-target dict (fast, reliable,
                # zero extra HTTP calls) instead of the old flaky per-drug ChEMBL lookup.
                if not entry.get("chembl_id") and entry.get("dgidb_approved") and not entry.get("indications"):
                    ot_entry = ot_known_drugs.get(entry["drug_name"].strip().lower())
                    extra_inds = ot_entry["indications"] if ot_entry else []
                    if extra_inds:
                        entry["indications"] = extra_inds
                        for ind in extra_inds:
                            term = ind.get("term", "")
                            all_diseases.append(term)
                            if any(kw in term.lower() for kw in CANCER_KEYWORDS):
                                chembl_cancer_inds.append(term)
                        specific_cancer_inds = [c for c in set(chembl_cancer_inds)
                                                if c.lower() not in GENERIC_TERMS]

                entry["fda_approved"] = fda
                entry["cancer_indication"] = len(specific_cancer_inds) > 0 or len(chembl_cancer_inds) > 0
                entry["diseases"] = "; ".join(all_diseases) if all_diseases else ""
                entry["cancer_types"] = "; ".join(specific_cancer_inds[:8]) if specific_cancer_inds else ""

                score, label, indicated_for_query, query_phase, other_cancer_inds = \
                    scorer.score_cancer_specific(entry, cancer_resolver, cancer_info)

                entry["score"] = score
                entry["score_label"] = label
                entry["indicated_for_query"] = indicated_for_query
                entry["query_cancer_phase"] = query_phase
                entry["other_cancer_inds"] = other_cancer_inds

                if prism_sens:
                    entry["prism_sensitivity"] = f"{prism_sens['pct_sensitive']}% sensitive"
                    entry["prism_mean_auc"] = prism_sens["mean_auc"]
                    entry["prism_top_lineages"] = prism_sens["top_sensitive_lineages"]
                else:
                    entry["prism_sensitivity"] = ""
                    entry["prism_mean_auc"] = None
                    entry["prism_top_lineages"] = ""

                entry["sources_str"] = "; ".join(entry["sources"])

            else:
                # ---- Original 3-tier scoring ----
                drug_cancer_types = specific_cancer_inds[:8]
                cancer_source = "ChEMBL_indications" if drug_cancer_types else ""

                if not drug_cancer_types and cancer_names and fda:
                    drug_cancer_types = [c for c in cancer_names[:5]
                                         if c.lower() not in GENERIC_TERMS][:5]
                    cancer_source = "gene_association"

                has_cancer_ind = len(drug_cancer_types) > 0 or len(cancer_names) > 0

                if prism_sens:
                    has_cancer_ind = True
                    if not drug_cancer_types:
                        drug_cancer_types = ["Cancer cell lines (PRISM)"]

                entry["fda_approved"] = fda
                entry["cancer_indication"] = has_cancer_ind
                entry["diseases"] = "; ".join(all_diseases) if all_diseases else ""
                entry["cancer_types"] = "; ".join(drug_cancer_types) if drug_cancer_types else ""
                entry["cancer_source"] = cancer_source
                entry["gene_cancer_associations"] = "; ".join(cancer_names[:10]) if cancer_names else ""
                entry["score"] = scorer.score(entry, gene_cancers)

                if prism_sens:
                    entry["prism_sensitivity"] = f"{prism_sens['pct_sensitive']}% sensitive"
                    entry["prism_mean_auc"] = prism_sens["mean_auc"]
                    entry["prism_top_lineages"] = prism_sens["top_sensitive_lineages"]
                else:
                    entry["prism_sensitivity"] = ""
                    entry["prism_mean_auc"] = None
                    entry["prism_top_lineages"] = ""

                entry["sources_str"] = "; ".join(entry["sources"])

        # Build result rows
        for key, entry in merged.items():
            if cancer_type:
                all_results.append({
                    "Gene": gene_info["symbol"],
                    "Drug_Name": entry["drug_name"],
                    "Drug_ChEMBL_ID": entry["chembl_id"],
                    "Max_Phase": str(entry["max_phase"]) if entry["max_phase"] >= 0 else ("N/A (DGIdb)" if entry["dgidb_approved"] else ""),
                    "FDA_Approved": "Yes" if entry["fda_approved"] else "No",
                    "Interaction_Type": entry["interaction_type"],
                    "Best_Potency_nM": entry["best_potency_nM"] if entry["best_potency_nM"] else None,
                    "Potency_Type": entry["potency_type"],
                    "Query_Cancer": cancer_info["disease_name"],
                    "Indicated_for_Query_Cancer": "Yes" if entry.get("indicated_for_query") else "No",
                    "Query_Cancer_Indication_Phase": entry.get("query_cancer_phase") or None,
                    "Other_Cancer_Indications": "; ".join(entry.get("other_cancer_inds", [])) if entry.get("other_cancer_inds") else "",
                    "Cancer_Types": entry["cancer_types"],
                    "Diseases": entry["diseases"],
                    "Score": entry["score"],
                    "Score_Label": entry.get("score_label", ""),
                    "PRISM_Sensitivity": entry["prism_sensitivity"],
                    "PRISM_Mean_AUC": entry["prism_mean_auc"],
                    "PRISM_Top_Lineages": entry["prism_top_lineages"],
                    "Sources": entry["sources_str"],
                })
            else:
                all_results.append({
                    "Gene": gene_info["symbol"],
                    "Drug_Name": entry["drug_name"],
                    "Drug_ChEMBL_ID": entry["chembl_id"],
                    "Max_Phase": str(entry["max_phase"]) if entry["max_phase"] >= 0 else ("N/A (DGIdb)" if entry["dgidb_approved"] else ""),
                    "FDA_Approved": "Yes" if entry["fda_approved"] else "No",
                    "Interaction_Type": entry["interaction_type"],
                    "Best_Potency_nM": entry["best_potency_nM"] if entry["best_potency_nM"] else None,
                    "Potency_Type": entry["potency_type"],
                    "Cancer_Indication": "Yes" if entry["cancer_indication"] else "No",
                    "Cancer_Source": entry.get("cancer_source", ""),
                    "Diseases": entry["diseases"],
                    "Cancer_Types": entry["cancer_types"],
                    "Gene_Cancer_Associations": entry["gene_cancer_associations"],
                    "PRISM_Sensitivity": entry["prism_sensitivity"],
                    "PRISM_Mean_AUC": entry["prism_mean_auc"],
                    "PRISM_Top_Lineages": entry["prism_top_lineages"],
                    "Score": entry["score"],
                    "Sources": entry["sources_str"],
                })

        print(f"  Scored {len(merged)} drugs for {gene}")

    # Create DataFrame
    df = pd.DataFrame(all_results)

    if df.empty:
        print("\n[ERROR] No results found for any gene. Check gene symbols and API connectivity.")
        return

    # Sort by gene, then score descending
    df = df.sort_values(["Gene", "Score", "Drug_Name"], ascending=[True, False, True])

    # Save full results (all scores)
    csv_full_path = output_dir / "drug_target_results_full.csv"
    df.to_csv(csv_full_path, index=False)

    # Save filtered results (min_score and above) — the main deliverable
    df_filtered = df[df["Score"] >= min_score].copy()
    csv_path = output_dir / "drug_target_results.csv"
    df_filtered.to_csv(csv_path, index=False)

    print(f"\nMain results (Score >= {min_score}) saved to {csv_path}")
    print(f"  {len(df_filtered)} gene-drug pairs")
    print(f"Full results (all scores) saved to {csv_full_path}")
    print(f"  {len(df)} gene-drug pairs total")

    if cancer_type:
        print(f"\nScore 5 (FDA-approved for {cancer_type}): {(df['Score'] == 5).sum()}")
        print(f"Score 4 (Clinical trial for {cancer_type}): {(df['Score'] == 4).sum()}")
        print(f"Score 3 (FDA-approved for other cancers): {(df['Score'] == 3).sum()}")
        print(f"Score 2 (FDA non-cancer / trial any cancer): {(df['Score'] == 2).sum()}")
        print(f"Score 1 (Experimental): {(df['Score'] == 1).sum()}")
    else:
        print(f"\nScore 3 (FDA + cancer): {(df['Score'] == 3).sum()}")
        print(f"Score 2 (FDA any / clinical + cancer): {(df['Score'] == 2).sum()}")
        print(f"Score 1 (experimental): {(df['Score'] == 1).sum()}")

    # Save unresolved genes
    if unresolved:
        unres_path = output_dir / "unresolved_genes.txt"
        unres_path.write_text("\n".join(unresolved))
        print(f"Unresolved genes saved to {unres_path}: {unresolved}")

    # Generate visualizations from filtered data
    if not skip_plots:
        print("\nGenerating visualizations...")
        hm = HeatmapGenerator(output_dir)

        if cancer_type:
            # Cancer-specific mode: heatmap (vmax=5) + stacked bar chart
            hm1_path = hm.genes_drugs_heatmap(
                df_filtered,
                max_drugs=max_drugs_per_gene * len(df_filtered["Gene"].unique()),
                vmax=5)
            if hm1_path:
                print(f"  Genes x Drugs heatmap: {hm1_path}")

            bc_path = hm.gene_scores_barchart(df_filtered, cancer_info["disease_name"], vmax=5)
            if bc_path:
                print(f"  Gene scores bar chart: {bc_path}")
        else:
            # Original mode: two heatmaps
            hm1_path = hm.genes_drugs_heatmap(
                df_filtered,
                max_drugs=max_drugs_per_gene * len(df_filtered["Gene"].unique()))
            if hm1_path:
                print(f"  Genes x Drugs heatmap: {hm1_path}")

            hm2_path = hm.genes_cancers_heatmap(df_filtered)
            if hm2_path:
                print(f"  Genes x Cancers heatmap: {hm2_path}")

    print("\nPipeline complete.")
    return df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Drug-Target Scoring Pipeline: find drugs targeting a gene list, "
                    "annotate FDA approval and cancer indications, score, and visualize."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--genes", type=str,
                       help="Comma-separated gene symbols (e.g., 'TP53,BRCA1,EGFR')")
    group.add_argument("--file", type=str,
                       help="Path to a text file with one gene symbol per line")
    parser.add_argument("--output", type=str, default="/mnt/results/",
                        help="Output directory for CSV and heatmaps (default: /mnt/results/)")
    parser.add_argument("--prism", action="store_true",
                        help="Enable PRISM/DepMap drug sensitivity enrichment (requires ~100MB download)")
    parser.add_argument("--min-score", type=int, default=2,
                        help="Minimum score to include in the main CSV (default: 2, i.e. FDA-approved or clinical trial). "
                             "Full results with all scores saved separately.")
    parser.add_argument("--max-drugs-per-gene", type=int, default=30,
                        help="Maximum drugs per gene to show in heatmaps (default: 30, keeps top-scoring)")
    parser.add_argument("--cancer", type=str, default=None,
                        help="Cancer type for cancer-specific 5-tier scoring (e.g., 'Neuroblastoma', 'Lung cancer'). "
                             "When provided, scores drugs by whether they are approved/trialed for YOUR cancer type. "
                             "When omitted, uses original 3-tier cancer-agnostic scoring.")

    args = parser.parse_args()

    # Parse gene list
    if args.genes:
        genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    else:
        gene_file = Path(args.file)
        if not gene_file.exists():
            print(f"[ERROR] Gene file not found: {args.file}")
            sys.exit(1)
        genes = [line.strip() for line in gene_file.read_text().splitlines()
                 if line.strip() and not line.startswith("#")]

    if not genes:
        print("[ERROR] No genes provided.")
        sys.exit(1)

    print(f"Drug-Target Scoring Pipeline")
    print(f"=" * 50)
    print(f"Genes: {len(genes)}")
    print(f"Output: {args.output}")
    print(f"PRISM: {'enabled' if args.prism else 'disabled'}")
    print(f"Min score for main CSV: {args.min_score}")
    if args.cancer:
        print(f"Cancer-specific mode: {args.cancer}")
    print(f"=" * 50)

    run_pipeline(genes, args.output, use_prism=args.prism,
                 min_score=args.min_score, max_drugs_per_gene=args.max_drugs_per_gene,
                 cancer_type=args.cancer)


if __name__ == "__main__":
    main()
