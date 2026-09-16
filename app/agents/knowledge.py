"""Literature + protein/drug-target knowledge for candidate genes.

Keyless public APIs, all read-only:
  - Europe PMC   : literature (gene x cancer type)
  - UniProt      : protein function, domains, sites, subcellular location, PDB xrefs
  - ChEMBL       : known drugs / mechanisms (mapped via UniProt accession)

Everything is cached in-process. Never writes user_data/. Network failures
degrade gracefully to empty results (the caller shows what it has).
"""
from __future__ import annotations

import json
import urllib.parse
import urllib.request
from functools import lru_cache

_UA = {"User-Agent": "INO80-DepMap-Explorer/1.0", "Accept": "application/json"}


def _get(url: str, timeout: float = 25.0) -> dict:
    req = urllib.request.Request(url, headers=_UA)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


# ---------------- Europe PMC ----------------
@lru_cache(maxsize=512)
def literature(gene: str, disease: str, n: int = 6) -> list[dict]:
    """Recent literature for gene x disease. Returns list of paper dicts."""
    query = f'("{gene}") AND ("{disease}")'
    q = urllib.parse.quote(query)
    sort = urllib.parse.quote("P_PDATE_D desc")
    url = (f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}"
           f"&format=json&pageSize={n}&resultType=lite&sort={sort}")
    try:
        data = _get(url)
    except Exception:
        return []
    out = []
    for r in data.get("resultList", {}).get("result", []):
        out.append({
            "title": r.get("title", "").strip(),
            "authors": r.get("authorString", ""),
            "journal": r.get("journalTitle", "") or r.get("bookOrReportDetails", ""),
            "year": r.get("pubYear", ""),
            "pmid": r.get("pmid", ""),
            "pmcid": r.get("pmcid", ""),
            "doi": r.get("doi", ""),
            "cited_by": r.get("citedByCount", 0),
            "url": (f"https://pubmed.ncbi.nlm.nih.gov/{r.get('pmid')}/" if r.get("pmid")
                    else (f"https://doi.org/{r.get('doi')}" if r.get("doi") else "")),
        })
    return out


# ---------------- UniProt ----------------
def _uniprot_feature_texts(feats: list, ftype: str, limit: int = 6) -> list[str]:
    out = []
    for f in feats or []:
        if f.get("type") == ftype:
            desc = f.get("description") or ftype
            loc = f.get("location", {})
            pos = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            span = f"{pos}" if pos == end else f"{pos}-{end}"
            out.append(f"{desc} ({span})" if pos else desc)
        if len(out) >= limit:
            break
    return out


@lru_cache(maxsize=512)
def uniprot(gene: str) -> dict | None:
    """Reviewed human UniProt entry for a gene symbol, distilled for drug design."""
    fields = ("accession,protein_name,cc_function,ft_domain,ft_act_site,ft_binding,"
              "cc_subcellular_location,xref_pdb,keyword,protein_existence")
    q = urllib.parse.quote(f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true")
    url = f"https://rest.uniprot.org/uniprotkb/search?query={q}&format=json&fields={fields}&size=1"
    try:
        data = _get(url)
    except Exception:
        return None
    res = data.get("results") or []
    if not res:
        return None
    e = res[0]
    acc = e.get("primaryAccession", "")
    name = ((e.get("proteinDescription", {}).get("recommendedName", {})
             .get("fullName", {}) or {}).get("value", ""))
    # function comment
    func = ""
    for c in e.get("comments", []):
        if c.get("commentType") == "FUNCTION":
            texts = c.get("texts", [])
            if texts:
                func = texts[0].get("value", "")
                break
    subcell = []
    for c in e.get("comments", []):
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            for loc in c.get("subcellularLocations", []):
                v = loc.get("location", {}).get("value")
                if v:
                    subcell.append(v)
    feats = e.get("features", [])
    pdb = [x["id"] for x in e.get("uniProtKBCrossReferences", []) if x.get("database") == "PDB"]
    kws = [k.get("name") for k in e.get("keywords", []) if k.get("name")]
    return {
        "accession": acc,
        "protein_name": name,
        "function": func,
        "domains": _uniprot_feature_texts(feats, "Domain"),
        "active_sites": _uniprot_feature_texts(feats, "Active site"),
        "binding_sites": _uniprot_feature_texts(feats, "Binding site"),
        "subcellular_location": subcell,
        "pdb_ids": pdb[:12],
        "n_pdb": len(pdb),
        "keywords": kws[:12],
        "uniprot_url": f"https://www.uniprot.org/uniprotkb/{acc}/entry" if acc else "",
    }


# ---------------- ChEMBL (via UniProt accession) ----------------
@lru_cache(maxsize=512)
def chembl_drugs(accession: str) -> dict:
    """Known drugs/mechanisms for a target identified by UniProt accession."""
    if not accession:
        return {"target_chembl_id": "", "n_mechanisms": 0, "mechanisms": []}
    try:
        t = _get("https://www.ebi.ac.uk/chembl/api/data/target.json"
                 f"?target_components__accession={accession}&limit=20")
    except Exception:
        return {"target_chembl_id": "", "n_mechanisms": 0, "mechanisms": []}
    single = [x for x in t.get("targets", []) if x.get("target_type") == "SINGLE PROTEIN"]
    if not single:
        return {"target_chembl_id": "", "n_mechanisms": 0, "mechanisms": []}
    tid = single[0]["target_chembl_id"]
    try:
        m = _get("https://www.ebi.ac.uk/chembl/api/data/mechanism.json"
                 f"?target_chembl_id={tid}&limit=1000")
    except Exception:
        return {"target_chembl_id": tid, "n_mechanisms": 0, "mechanisms": [],
                "chembl_url": f"https://www.ebi.ac.uk/chembl/target_report_card/{tid}/"}
    mechs = m.get("mechanisms", [])
    distilled, seen = [], set()
    for x in mechs:
        moa = x.get("mechanism_of_action")
        mol = x.get("molecule_chembl_id")
        key = (moa, mol)
        if moa and key not in seen:
            seen.add(key)
            distilled.append({"moa": moa, "molecule_chembl_id": mol,
                              "action_type": x.get("action_type", "")})
    return {
        "target_chembl_id": tid,
        "n_mechanisms": m.get("page_meta", {}).get("total_count", len(mechs)),
        "mechanisms": distilled[:12],
        "chembl_url": f"https://www.ebi.ac.uk/chembl/target_report_card/{tid}/",
    }


# ---------------- combined dossier ----------------
def _druggability(up: dict | None, ch: dict) -> str:
    bits = []
    if ch.get("n_mechanisms", 0) > 0:
        bits.append(f"{ch['n_mechanisms']} known drug mechanism(s) in ChEMBL")
    else:
        bits.append("no approved-drug mechanism in ChEMBL")
    if up:
        if up.get("n_pdb", 0) > 0:
            bits.append(f"{up['n_pdb']} PDB structure(s)")
        else:
            bits.append("no PDB structure")
        kw = " ".join(up.get("keywords", [])).lower()
        cls = [c for c in ("kinase", "protease", "transferase", "hydrolase",
                           "receptor", "transcription", "ion channel", "nuclease")
               if c in kw or c in (up.get("protein_name", "").lower())]
        if cls:
            bits.append("class: " + ", ".join(sorted(set(cls))))
    return "; ".join(bits)


def dossier(gene: str, disease: str, n_papers: int = 6) -> dict:
    """Full knowledge dossier for one gene: literature + protein + drugs + summary."""
    up = uniprot(gene)
    ch = chembl_drugs(up["accession"]) if up else {"n_mechanisms": 0, "mechanisms": []}
    return {
        "gene": gene,
        "literature": literature(gene, disease, n_papers),
        "uniprot": up,
        "chembl": ch,
        "druggability": _druggability(up, ch),
    }
