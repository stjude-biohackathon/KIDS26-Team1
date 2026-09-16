"""Minimal UCSC Xena client (plain urllib — no `xenaPython` dependency).

Xena hubs accept a small Lisp/Clojure query language POSTed as text/plain to
`{hub}/data/`. This module hand-implements the handful of query primitives the
tissue-specificity pipeline needs, reverse-engineered from `xenaPython`'s own
`.xq` templates and `xenaQuery.dataset_gene_str`. Verified live against
https://toil.xenahubs.net to return identical results to the real package
(field_codes index 89 == "Neuroblastoma", matching the upstream pipeline).

Kept dependency-free to match this app's existing keyless-evidence-source
philosophy (Europe PMC / UniProt / ChEMBL / STRING all use bare urllib too).
"""
from __future__ import annotations

import json
import urllib.error
import urllib.request

_UA = {"User-Agent": "SCRAP-AI/1.0", "Content-Type": "text/plain"}

_FIELD_CODES_XQ = """; fieldCodes
(fn [dataset fields]
	(query
	  {:select [:P.name [#sql/call [:group_concat :value :order :ordering :separator #sql/call [:chr 9]] :code]]
	   :from [[{:select [:field.id :field.name]
				:from [:field]
				:join [{:table [[[:name :varchar fields]] :T]} [:= :T.name :field.name]]
				:where [:= :dataset_id {:select [:id]
								 :from [:dataset]
								 :where [:= :name dataset]}]} :P]]
	   :left-join [:code [:= :P.id :field_id]]
	   :group-by [:P.id]}))
"""

_DATASET_SAMPLES_XQ = """; datasetSamples
(fn [dataset limit]
    (map :value
      (query
        {:select [:value]
         :from [:dataset]
         :join [:field [:= :dataset.id :dataset_id]
                :code [:= :field.id :field_id]]
         :limit limit
         :where [:and
                 [:= :dataset.name dataset]
                 [:= :field.name "sampleID"]]})))
"""

_DATASET_FETCH_XQ = """; datasetFetch
(fn [dataset samples probes]
  (fetch [{:table dataset
           :columns probes
           :samples samples}]))
"""

# Not wrapped via the generic call() mechanism upstream -- already a full
# standalone `(let [...] ...)` expression.
_DATASET_GENE_STR = """
(let [probemap (:probemap (car (query {:select [:probemap]
                                      :from [:dataset]
                                      :where [:= :name %s]})))
     probes-for-gene (fn [gene] ((xena-query {:select ["name"] :from [probemap] :where [:in :any "genes" [gene]]}) "name"))
     avg (fn [scores] (mean scores 0))
     scores-for-gene (fn [gene]
         (let [probes (probes-for-gene gene)
               scores (fetch [{:table %s
                               :samples %s
                               :columns (probes-for-gene gene)}])]
           {:gene gene
            :scores (if (car probes) (avg scores) [[]])}))]
 (map scores-for-gene %s))
"""

GENE_BATCH = 20  # genes per API call (matches upstream pipeline)


def _quote(s: str) -> str:
    return '"' + s + '"'


def _array(items: list[str]) -> str:
    return "[" + ", ".join(_quote(s) for s in items) + "]"


def _marshall(p):
    if p is None:
        return "nil"
    if isinstance(p, str):
        return _quote(p)
    return str(p)


def _post(hub: str, query: str, timeout: float = 40.0):
    req = urllib.request.Request(hub + "/data/", query.encode(), headers=_UA)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def field_codes(hub: str, dataset: str, fields: list[str]) -> dict[str, list[str]]:
    """Return {field_name: [category labels in code order]}."""
    q = "(%s %s %s)" % (_FIELD_CODES_XQ, _quote(dataset), _array(fields))
    rows = _post(hub, q)
    return {r["name"]: r["code"].split("\t") for r in rows}


def dataset_samples(hub: str, dataset: str, limit: int | None = None) -> list[str]:
    q = "(%s %s %s)" % (_DATASET_SAMPLES_XQ, _quote(dataset), _marshall(limit))
    return _post(hub, q)


def dataset_fetch(hub: str, dataset: str, samples: list[str], columns: list[str]):
    """Returns Xena's raw {'probes': [...], 'sampleIDs': [...], 'nan': v}-shaped payload."""
    q = "(%s %s %s %s)" % (_DATASET_FETCH_XQ, _quote(dataset), _array(samples), _array(columns))
    return _post(hub, q)


def dataset_gene_values(hub: str, dataset: str, samples: list[str],
                        genes: list[str]) -> dict[str, list[float | None]]:
    """Batched (GENE_BATCH genes/call) per-gene average-probe values for samples.

    Returns {gene: [value_per_sample...]} in the same order as `samples`.
    """
    out: dict[str, list] = {}
    for i in range(0, len(genes), GENE_BATCH):
        batch = genes[i:i + GENE_BATCH]
        q = _DATASET_GENE_STR % (_quote(dataset), _quote(dataset), _array(samples), _array(batch))
        try:
            res = _post(hub, q, timeout=90.0)
        except urllib.error.HTTPError:
            for g in batch:
                out[g] = [None] * len(samples)
            continue
        for entry in res:
            row = entry["scores"][0] if entry.get("scores") else None
            # Xena returns scores:[[]] (empty) when the gene isn't in this dataset's
            # probemap at all -- pad/replace with all-None so every column stays
            # aligned to `samples`, regardless of which genes resolved.
            if not row or len(row) != len(samples):
                row = [None] * len(samples)
            out[entry["gene"]] = row
        # genes Xena didn't return an entry for at all (dropped from `res`)
        for g in batch:
            if g not in out:
                out[g] = [None] * len(samples)
    return out
