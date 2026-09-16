"""Tests for GTEx/Xena tissue-specificity pipeline.

Offline: gate/score logic on synthetic data (no network).
Network-gated: live Xena/GTEx calls, skip gracefully if unreachable.
Run: uv run pytest app/tests -q
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

APP = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(APP))

from agents import gtex_agent  # noqa: E402


def _online() -> bool:
    try:
        from agents import xena_client
        xena_client.field_codes("https://toil.xenahubs.net", "TcgaTargetGTEX_phenotype.txt",
                                ["detailed_category"])
        return True
    except Exception:
        return False


# ---------------- offline: gate logic on synthetic data ----------------
def test_normal_gate_median_tolerant_of_one_outlier():
    # one atypical elevated subregion should not fail the median-mode gate,
    # but a gene elevated across most tissues should fail under both modes
    df = pd.DataFrame({
        "Adrenal_Gland": [0.3, 20.0],
        "Brain_A": [2.0, 22.0],
        "Brain_B": [2.5, 18.0],
        "Brain_C": [6.7, 1.0],  # outlier only for row 0
    }, index=["GENE_A", "GENE_B"])
    median_gate = gtex_agent.normal_gate(df, cutoff=5.0, mode="median")
    all_gate = gtex_agent.normal_gate(df, cutoff=5.0, mode="all")
    assert median_gate["GENE_A"] == True  # noqa: E712  (median ~2.25 < 5)
    assert all_gate["GENE_A"] == False  # noqa: E712    (Brain_C=6.7 fails "all")
    assert median_gate["GENE_B"] == False and all_gate["GENE_B"] == False  # noqa: E712


def test_normal_gate_empty_is_empty():
    out = gtex_agent.normal_gate(pd.DataFrame(), cutoff=5.0)
    assert out.empty


# ---------------- offline: suggestion helpers (no network needed for lineage hints) ----------------
def test_lineage_tissue_hints_present_for_common_lineages():
    from agents import tissue_specificity as ts
    assert "CNS/Brain" in ts.LINEAGE_TISSUE_HINTS
    assert "Lung" in ts.LINEAGE_TISSUE_HINTS
    assert ts.LINEAGE_TISSUE_HINTS["CNS/Brain"] == ["Brain_"]


# ---------------- network-gated: live wire-protocol + full pipeline ----------------
def test_xena_field_codes_matches_known_index():
    if not _online():
        pytest.skip("offline")
    from agents import xena_client
    codes = xena_client.field_codes("https://toil.xenahubs.net", "TcgaTargetGTEX_phenotype.txt",
                                    ["detailed_category"])
    labels = codes["detailed_category"]
    assert len(labels) > 50
    assert labels[89] == "Neuroblastoma"  # matches the collaborator's hardcoded NBL_CODE=89


def test_xena_gene_values_shape():
    if not _online():
        pytest.skip("offline")
    from agents import xena_client
    samples = xena_client.dataset_samples("https://toil.xenahubs.net",
                                          "TcgaTargetGtex_rsem_gene_tpm", 3)
    assert len(samples) == 3
    vals = xena_client.dataset_gene_values("https://toil.xenahubs.net",
                                           "TcgaTargetGtex_rsem_gene_tpm", samples, ["MYCN"])
    assert "MYCN" in vals and len(vals["MYCN"]) == 3


def test_full_triage_reproduces_collaborator_validated_result():
    """The exact example from user_data/gtex_other_expression/README.md:
    'Expect PHOX2B and MYCN at rank 1-2 passing all gates.'"""
    if not _online():
        pytest.skip("offline")
    from agents import tissue_specificity as ts
    genes = ["INO80", "MYCN", "PHOX2B", "ACTL6A", "ACTR8", "MCRS1", "TP53"]
    tissues = ["Adrenal_Gland"] + [t for t in gtex_agent.gtex_tissues() if t.startswith("Brain_")]
    table, heat = ts.triage(genes, tissues, ["Neuroblastoma"], normal_cutoff=5.0)
    top2 = table.sort_values("selectivity_score", ascending=False).head(2).index.tolist()
    assert set(top2) == {"PHOX2B", "MYCN"}
    assert table.loc["PHOX2B", "pass_all"] == True  # noqa: E712
    assert table.loc["MYCN", "pass_all"] == True  # noqa: E712
    assert not heat.empty


def test_suggest_target_labels_finds_neuroblastoma():
    if not _online():
        pytest.skip("offline")
    from agents import tissue_specificity as ts
    hits = ts.suggest_target_labels("Neuroblastoma")
    assert "Neuroblastoma" in hits


def test_pancancer_gate_excludes_target_from_other_cancers():
    """If the chosen target happens to itself be a TCGA type, it must be excluded
    from the 'other cancers must be low' universe (no self-comparison, no duplicate
    column in the heatmap)."""
    if not _online():
        pytest.skip("offline")
    from agents import tissue_specificity as ts
    genes = ["EGFR"]
    tissues = ["Lung"]
    table, heat = ts.triage(genes, tissues, ["Lung Adenocarcinoma"], normal_cutoff=5.0)
    assert list(heat.columns).count("Lung Adenocarcinoma") == 1
    assert len(table) == 1


def test_tissue_heatmap_plot_shape():
    from viz import plots
    import plotly.graph_objects as go
    heat = pd.DataFrame([[1.2, 0.1, 30.0], [0.0, 0.0, 12.0]],
                        index=["GENE_A", "GENE_B"], columns=["Adrenal_Gland", "Lung Adenocarcinoma", "Neuroblastoma"])
    fig = plots.tissue_heatmap(heat, "t")
    assert isinstance(fig, go.Figure)
    assert isinstance(plots.tissue_heatmap(pd.DataFrame(), "t"), go.Figure)  # empty-safe


def test_dataset_gene_values_pads_empty_scores(monkeypatch):
    """Regression: Xena returns scores:[[]] (empty) for a gene not in this
    dataset's probemap -- must be padded to len(samples), not left length-0
    (that mismatch crashed pd.DataFrame construction downstream)."""
    from agents import xena_client

    samples = ["s1", "s2", "s3"]

    def fake_post(hub, query, timeout=90.0):
        # simulate: GENE_A resolves fine, GENE_B is unresolved (empty scores),
        # GENE_C is missing from the response entirely
        return [
            {"gene": "GENE_A", "scores": [[1.0, 2.0, 3.0]]},
            {"gene": "GENE_B", "scores": [[]]},
        ]

    monkeypatch.setattr(xena_client, "_post", fake_post)
    out = xena_client.dataset_gene_values("http://fake", "ds", samples, ["GENE_A", "GENE_B", "GENE_C"])

    assert len(out["GENE_A"]) == 3 and out["GENE_A"] == [1.0, 2.0, 3.0]
    assert len(out["GENE_B"]) == 3 and out["GENE_B"] == [None, None, None]
    assert len(out["GENE_C"]) == 3 and out["GENE_C"] == [None, None, None]

    # the exact downstream construction that used to crash with
    # "Length of values (0) does not match length of index (N)"
    import pandas as pd
    df = pd.DataFrame({g: vals for g, vals in out.items()}, index=samples).T
    assert df.shape == (3, 3)


def test_tcga_target_expression_handles_unresolved_gene(monkeypatch):
    """End-to-end through tissue_specificity: an unresolved gene must not crash
    the batch, just produce NaN/False for that gene's row."""
    from agents import tissue_specificity as ts, xena_client

    samples = ["GTEX-1", "TCGA-1", "TCGA-2"]

    def fake_dataset_samples(hub, dataset, limit=None):
        return samples

    def fake_dataset_fetch(hub, dataset, smp, fields):
        # _study, detailed_category, _sample_type
        return [[2, 0, 0], [None, 5, 10], [None, 5, 10]]

    def fake_dataset_gene_values(hub, dataset, smp, genes):
        out = {}
        for g in genes:
            if g == "REALGENE":
                out[g] = [1.0, 5.0, 6.0]
            else:
                out[g] = []  # simulate an unresolved gene from the raw client too
        return out

    monkeypatch.setattr(xena_client, "dataset_samples", fake_dataset_samples)
    monkeypatch.setattr(xena_client, "dataset_fetch", fake_dataset_fetch)
    monkeypatch.setattr(xena_client, "dataset_gene_values", fake_dataset_gene_values)
    ts.target_label_options.cache_clear()
    ts._sample_metadata.cache_clear()
    monkeypatch.setattr(ts, "target_label_options", lambda: ["FakeDisease"])

    xd = ts._tcga_target_expression(("REALGENE", "FAKEGENE"), ("FakeDisease",))
    assert xd["target"].shape[0] == 2  # both genes get a row
    assert not xd["target"].loc["FAKEGENE"].notna().any()  # unresolved -> all-NaN, no crash


def test_triage_graded_tier_and_sort_order_offline(monkeypatch):
    """Fully offline (mocked Xena/GTEx calls): verifies n_gates_passed/specificity_tier
    for genes passing 0, 1, 2, and 3 of the 3 gates, and that the table is ranked by
    n_gates_passed first (a graded alternative to the strict pass_all AND that used to
    lump every non-perfect gene into the same 'False' bucket)."""
    from agents import tissue_specificity as ts

    genes = ["GENE0", "GENE1", "GENE2", "GENE3"]

    def fake_load_gtex(gene_symbols, tissues):
        # GENE3 alone is silent in the normal tissue; the rest are not.
        return pd.DataFrame(
            {"Tissue_A": {"GENE0": 10.0, "GENE1": 10.0, "GENE2": 10.0, "GENE3": 1.0}}
        )

    def fake_tcga_target_expression(g, target_labels):
        # GENE2 and GENE3 are low in both other cancer types; GENE0/GENE1 are not.
        tcga_type_median = pd.DataFrame(
            {"OtherA": {"GENE0": 10.0, "GENE1": 10.0, "GENE2": 1.0, "GENE3": 1.0},
             "OtherB": {"GENE0": 10.0, "GENE1": 10.0, "GENE2": 1.0, "GENE3": 1.0}}
        )
        # GENE1/GENE2/GENE3 are all high in the target cohort; GENE0 is not.
        target = pd.DataFrame(
            {"s1": {"GENE0": 1.0, "GENE1": 20.0, "GENE2": 20.0, "GENE3": 20.0},
             "s2": {"GENE0": 1.0, "GENE1": 20.0, "GENE2": 20.0, "GENE3": 20.0}}
        )
        return {"target": target, "tcga_type_median": tcga_type_median}

    monkeypatch.setattr(gtex_agent, "load_gtex", fake_load_gtex)
    monkeypatch.setattr(ts, "_tcga_target_expression", fake_tcga_target_expression)

    table, _heat = ts.triage(genes, ["Tissue_A"], ["MyTumor"],
                             normal_cutoff=5.0, pancancer_fraction=1.0,
                             target_median_tpm=10.0, target_fraction_expressed=0.50)

    assert table.loc["GENE0", "n_gates_passed"] == 0
    assert table.loc["GENE0", "specificity_tier"] == "fail"
    assert table.loc["GENE0", "pass_all"] == False  # noqa: E712

    assert table.loc["GENE1", "n_gates_passed"] == 1
    assert table.loc["GENE1", "specificity_tier"] == "low"

    assert table.loc["GENE2", "n_gates_passed"] == 2
    assert table.loc["GENE2", "specificity_tier"] == "medium"
    assert table.loc["GENE2", "pass_all"] == False  # noqa: E712

    assert table.loc["GENE3", "n_gates_passed"] == 3
    assert table.loc["GENE3", "specificity_tier"] == "high"
    assert table.loc["GENE3", "pass_all"] == True  # noqa: E712

    # ranked strictly by n_gates_passed (descending): 3 > 2 > 1 > 0
    assert table.index.tolist() == ["GENE3", "GENE2", "GENE1", "GENE0"]
