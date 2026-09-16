"""Tissue-specificity triage: GTEx normal-tissue silence + TCGA pan-cancer +
target-cohort expression, generalized from the collaborator's Neuroblastoma-
specific `triage.py` to work for any selected cancer type / gene list.

Three gates + a composite selectivity score:
  pass_normal     -- TPM < cutoff across the chosen GTEx normal tissue(s)
  pass_pancancer  -- TPM < cutoff in >= `pancancer_fraction` of OTHER TCGA types
  pass_target     -- median TPM >= target_median OR >= target_fraction of
                      target-cohort samples with TPM >= cutoff
  selectivity_score = target_median * pancancer_low_fraction / (1 + max(normal tissues))
  pass_all        -- pass_normal AND pass_pancancer AND pass_target (strict; unchanged,
                      kept exactly as originally defined -- this is what the
                      collaborator-validated PHOX2B/MYCN reference result checks)
  n_gates_passed  -- 0-3, how many of the three gates above a gene actually passes
  specificity_tier -- n_gates_passed mapped to fail(0) / low(1) / medium(2) / high(3);
                      a graded alternative to pass_all that distinguishes a gene passing
                      2/3 gates from one passing 0/3, which pass_all alone collapses into
                      the same "False" bucket. Table rows are ranked by n_gates_passed
                      first, selectivity_score second (a strict superset of the old
                      pass_all-first ordering: pass_all=True genes are always tier=high,
                      i.e. n_gates_passed=3, so they still sort first).

Network calls (Xena) are cached; safe to call repeatedly for the same genes.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd

from agents import gtex_agent, xena_client

HUB = "https://toil.xenahubs.net"
PHENO = "TcgaTargetGTEX_phenotype.txt"
EXPR = "TcgaTargetGtex_rsem_gene_tpm"
TCGA_TUMOR_SAMPLE_TYPES = {0, 1, 10}  # Primary/Recurrent/Primary Solid Tumor


def _log2_to_tpm(v):
    if v is None:
        return np.nan
    try:
        return 2.0 ** float(v) - 0.001
    except (TypeError, ValueError):
        return np.nan


@lru_cache(maxsize=1)
def target_label_options() -> list[str]:
    """All selectable target-cohort / cancer-type labels (Xena detailed_category)."""
    codes = xena_client.field_codes(HUB, PHENO, ["detailed_category"])
    return codes.get("detailed_category", [])


def suggest_target_labels(disease_or_lineage: str, top_n: int = 5) -> list[str]:
    """Best-effort keyword-overlap suggestion from the live Xena label list.

    Not authoritative -- always shown alongside the full list for manual pick.
    """
    q = {w for w in disease_or_lineage.lower().replace("-", " ").replace("/", " ").split() if len(w) > 2}
    scored = []
    for label in target_label_options():
        lw = {w for w in label.lower().replace("-", " ").split() if len(w) > 2}
        overlap = len(q & lw)
        if overlap:
            scored.append((overlap, label))
    scored.sort(key=lambda t: -t[0])
    return [lbl for _, lbl in scored[:top_n]]


# best-effort OncotreeLineage -> GTEx tissue-name substrings (not authoritative)
LINEAGE_TISSUE_HINTS: dict[str, list[str]] = {
    "CNS/Brain": ["Brain_"], "Lung": ["Lung"], "Skin": ["Skin_"],
    "Breast": ["Breast"], "Kidney": ["Kidney_"], "Pancreas": ["Pancreas"],
    "Bowel": ["Colon_", "Small_Intestine"], "Liver": ["Liver"],
    "Esophagus/Stomach": ["Esophagus_", "Stomach"], "Ovary/Fallopian Tube": ["Ovary"],
    "Uterus": ["Uterus"], "Testis": ["Testis"], "Prostate": ["Prostate"],
    "Myeloid": ["Whole_Blood", "Spleen"], "Lymphoid": ["Whole_Blood", "Spleen"],
    "Peripheral Nervous System": ["Adrenal_Gland", "Brain_"],
    "Bladder/Urinary Tract": ["Bladder"], "Thyroid": ["Thyroid"],
    "Cervix": ["Cervix_"], "Biliary Tract": ["Liver"], "Head and Neck": ["Salivary_Gland"],
}


def suggest_gtex_tissues(lineage: str) -> list[str]:
    """Best-effort GTEx tissue-name substrings for a lineage; expand to real columns."""
    hints = LINEAGE_TISSUE_HINTS.get(lineage, [])
    all_tissues = gtex_agent.gtex_tissues()
    out = []
    for h in hints:
        out += [t for t in all_tissues if t.startswith(h)]
    return sorted(set(out))


@lru_cache(maxsize=1)
def _sample_metadata() -> pd.DataFrame:
    """study/detailed_category/sample_type per sample in the TCGA+TARGET+GTEx cohort."""
    samples = xena_client.dataset_samples(HUB, PHENO, None)
    fields = ["_study", "detailed_category", "_sample_type"]
    data = xena_client.dataset_fetch(HUB, PHENO, samples, fields)
    # data[i] = list of values for fields[i], one per sample (same order as `samples`)
    labels = target_label_options()
    meta = pd.DataFrame({f: data[i] for i, f in enumerate(fields)}, index=samples)
    for c in fields:
        meta[c] = pd.to_numeric(meta[c], errors="coerce")

    def _label(code):
        if pd.isna(code):
            return None
        code = int(code)
        return labels[code] if 0 <= code < len(labels) else None

    meta["cancer_type"] = meta["detailed_category"].map(_label)
    meta["is_tcga_tumor"] = (meta["_study"] == 0) & (meta["_sample_type"].isin(TCGA_TUMOR_SAMPLE_TYPES))
    return meta


def _tcga_target_expression(genes: tuple[str, ...], target_labels: tuple[str, ...]) -> dict:
    meta = _sample_metadata()
    target_samples = meta.index[meta["cancer_type"].isin(target_labels)].tolist()
    tcga_samples = meta.index[meta["is_tcga_tumor"]].tolist()
    needed = list(dict.fromkeys(target_samples + tcga_samples))
    if not needed:
        return {"target": pd.DataFrame(), "tcga_type_median": pd.DataFrame(), "meta": meta}

    raw = xena_client.dataset_gene_values(HUB, EXPR, needed, list(genes))
    n = len(needed)
    # Defensive: force every column to len(needed) regardless of what the client
    # returned, so an unresolved gene (or any upstream API quirk) degrades to all-NaN
    # for that gene rather than crashing the whole batch.
    cols = {}
    for g in genes:
        vals = raw.get(g)
        if not vals or len(vals) != n:
            vals = [None] * n
        cols[g] = [_log2_to_tpm(v) for v in vals]
    expr = pd.DataFrame(cols, index=needed).T
    expr.index = expr.index.str.upper()

    target_mat = expr[target_samples] if target_samples else pd.DataFrame(index=expr.index)
    tcga_mat = expr[tcga_samples] if tcga_samples else pd.DataFrame(index=expr.index)
    tcga_types = meta.loc[tcga_samples, "cancer_type"] if tcga_samples else pd.Series(dtype=object)
    type_median = (tcga_mat.T.groupby(tcga_types.values).median().T
                  if tcga_samples else pd.DataFrame(index=expr.index))
    return {"target": target_mat, "tcga_type_median": type_median, "meta": meta}


def triage(gene_symbols: list[str], gtex_tissues: list[str], target_labels: list[str],
          normal_cutoff: float = 5.0, pancancer_fraction: float = 0.80,
          target_median_tpm: float = 10.0, target_fraction_expressed: float = 0.50,
          normal_gate_mode: str = "median"
          ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the generalized triage. Returns (ranked_table, heatmap_matrix_raw_tpm)."""
    genes = tuple(g.upper() for g in gene_symbols)
    gtex_tissues = tuple(gtex_tissues)
    target_labels = tuple(target_labels)

    gtex = gtex_agent.load_gtex(genes, gtex_tissues)
    pass_normal = gtex_agent.normal_gate(gtex, normal_cutoff, mode=normal_gate_mode)

    xd = _tcga_target_expression(genes, target_labels)
    tm = xd["tcga_type_median"]
    # exclude the target label(s) from the "other cancers must be low" universe
    other_cols = [c for c in tm.columns if c not in target_labels]
    tm_other = tm[other_cols] if other_cols else tm
    low_frac = (tm_other < normal_cutoff).mean(axis=1) if not tm_other.empty else pd.Series(dtype=float)
    pass_pancancer = low_frac >= pancancer_fraction

    target_mat = xd["target"]
    target_median = target_mat.median(axis=1) if not target_mat.empty else pd.Series(dtype=float)
    target_frac = ((target_mat >= normal_cutoff).mean(axis=1)
                   if not target_mat.empty else pd.Series(dtype=float))
    pass_target = (target_median >= target_median_tpm) | (target_frac >= target_fraction_expressed)

    table = pd.DataFrame(index=list(genes))
    for t in gtex_tissues:
        if t in gtex.columns:
            table[t] = gtex[t]
    table["normal_max"] = gtex.max(axis=1) if not gtex.empty else np.nan
    table["pass_normal"] = pass_normal.reindex(table.index)
    table["pancancer_low_fraction"] = low_frac.reindex(table.index)
    table["pass_pancancer"] = pass_pancancer.reindex(table.index).fillna(False)
    table["target_median_tpm"] = target_median.reindex(table.index)
    table["target_fraction_expressed"] = target_frac.reindex(table.index)
    table["pass_target"] = pass_target.reindex(table.index).fillna(False)
    table["pass_all"] = (table["pass_normal"].fillna(False) & table["pass_pancancer"]
                         & table["pass_target"])
    # Graded alternative to the strict pass_all AND: count how many of the 3 gates each
    # gene actually passes (0-3) and map it to a coarse fail/low/medium/high label. This
    # is purely additive -- pass_all keeps its exact original boolean meaning (needed for
    # the collaborator-validated PHOX2B/MYCN reference check below) -- but lets a gene
    # that passes 2/3 gates be distinguished from one that passes 0/3, which pass_all
    # alone collapses into the same "False" bucket.
    table["n_gates_passed"] = (
        table["pass_normal"].fillna(False).astype(int)
        + table["pass_pancancer"].fillna(False).astype(int)
        + table["pass_target"].fillna(False).astype(int)
    )
    table["specificity_tier"] = table["n_gates_passed"].map(
        {0: "fail", 1: "low", 2: "medium", 3: "high"})
    with np.errstate(divide="ignore", invalid="ignore"):
        table["selectivity_score"] = (
            table["target_median_tpm"].fillna(0) * table["pancancer_low_fraction"].fillna(0)
            / (1.0 + table["normal_max"].fillna(0))
        )
    table = table.sort_values(["n_gates_passed", "selectivity_score"], ascending=[False, False])
    table.insert(0, "rank", range(1, len(table) + 1))
    table.index.name = "gene"

    heat_cols = list(gtex_tissues) + other_cols + list(target_labels)
    heat = pd.DataFrame(index=table.index, columns=heat_cols, dtype=float)
    for t in gtex_tissues:
        if t in gtex.columns:
            heat[t] = gtex[t].reindex(table.index)
    for c in other_cols:
        heat[c] = tm[c].reindex(table.index)
    for tl in target_labels:
        heat[tl] = target_median.reindex(table.index) if len(target_labels) == 1 else np.nan
    return table, heat
