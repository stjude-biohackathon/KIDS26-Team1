"""Interactive Plotly figures for the DepMap explorer."""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

_NEGLOG_CAP = 50.0


def _neglog10(q: pd.Series) -> pd.Series:
    q = q.clip(lower=1e-300)
    return (-np.log10(q)).clip(upper=_NEGLOG_CAP)


def _add_gene_highlights(fig: go.Figure, df: pd.DataFrame, genes: list[str] | None,
                         x_col: str, y_col: str) -> go.Figure:
    """Overlay a distinct marker+label for specific genes on an existing scatter
    figure -- additive only (a new trace), never mutates the original data/trace,
    so existing callers that omit `genes` see byte-identical figures."""
    if not genes:
        return fig
    sub = df[df.index.isin(genes)] if df.index.name else df[df["gene"].isin(genes)]
    if sub.empty:
        return fig
    labels = list(sub.index) if df.index.name else list(sub["gene"])
    fig.add_trace(go.Scatter(
        x=sub[x_col], y=sub[y_col], mode="markers+text", text=labels,
        textposition="top center", name="requested gene(s)",
        marker=dict(size=16, color="rgba(0,0,0,0)", symbol="circle",
                   line=dict(width=3, color="#D11947")),
    ))
    return fig


def volcano(agg: pd.DataFrame, r_min: float, q_max: float, title: str,
            positive_only: bool = False, highlight_genes: list[str] | None = None) -> go.Figure:
    """x = mean_r, y = -log10(min_q); color = passes threshold."""
    if agg is None or agg.empty:
        return _empty("No candidates")
    df = agg.reset_index().copy()
    df["neglog_q"] = _neglog10(df["min_q"])
    score = df["mean_r"] if positive_only else df["mean_abs_r"]
    df["hit"] = (score.abs() >= r_min) & (df["min_q"] <= q_max) & (df["n_seeds_signif"] >= 1)
    df["status"] = np.where(df["hit"], "hit", "background")
    fig = px.scatter(
        df, x="mean_r", y="neglog_q", color="status",
        color_discrete_map={"hit": "#d62728", "background": "#c7c7c7"},
        custom_data=["gene", "best_seed", "best_r", "n_seeds_signif", "max_abs_r"],
        title=title, opacity=0.7,
    )
    fig.update_traces(
        hovertemplate=(
            "<b>%{customdata[0]}</b><br>mean_r=%{x:.3f}<br>-log10 q=%{y:.2f}"
            "<br>best seed=%{customdata[1]} (r=%{customdata[2]:.3f})"
            "<br>signif seeds=%{customdata[3]} · max|r|=%{customdata[4]:.3f}<extra></extra>"
        ),
        marker=dict(size=6),
    )
    fig.add_hline(y=-np.log10(q_max), line_dash="dash", line_color="#888",
                  annotation_text=f"q={q_max}")
    fig.add_vline(x=r_min, line_dash="dot", line_color="#888")
    fig.add_vline(x=-r_min, line_dash="dot", line_color="#888")
    fig.update_layout(xaxis_title="mean r vs seed set", yaxis_title="-log10(min q)",
                      legend_title="", height=520)
    return _add_gene_highlights(fig, df.set_index("gene"), highlight_genes, "mean_r", "neglog_q")


def heatmap(r_matrix: pd.DataFrame, top_k: int, agg: pd.DataFrame,
            title: str, order: list | None = None) -> go.Figure:
    """Candidate × seed correlation heatmap.

    If `order` is given, rows follow that gene order (already filtered/sorted by the
    caller, e.g. selective genes by weighted joint %). Otherwise top-K by mean_abs_r.
    """
    if r_matrix is None or r_matrix.empty:
        return _empty("No matrix")
    if order is not None:
        top = [g for g in order if g in r_matrix.index][:top_k]
    else:
        if agg is None or agg.empty:
            return _empty("No matrix")
        top = list(agg.head(top_k).index)
    sub = r_matrix.loc[[g for g in top if g in r_matrix.index]]
    if sub.empty:
        return _empty("No genes to show")
    fig = px.imshow(
        sub, color_continuous_scale="RdBu_r", zmin=-1, zmax=1, aspect="auto",
        labels=dict(x="seed gene", y="candidate", color="r"), title=title,
    )
    fig.update_traces(hovertemplate="candidate %{y}<br>seed %{x}<br>r=%{z:.3f}<extra></extra>")
    fig.update_layout(height=max(400, 18 * len(sub) + 120))
    return fig


def consensus_scatter(cons: pd.DataFrame, title: str) -> go.Figure:
    """co-expr mean_r (x) vs co-ess mean_r (y); size = combined |r|, color = common-essential."""
    if cons is None or cons.empty:
        return _empty("No consensus hits")
    df = cons.reset_index().copy()
    df["ce"] = np.where(df["common_essential"], "common-essential", "selective")
    fig = px.scatter(
        df, x="coexpr_mean_r", y="coess_mean_r", size="combined_abs_r",
        color="ce", color_discrete_map={"common-essential": "#ff7f0e", "selective": "#1f77b4"},
        custom_data=["gene", "coexpr_min_q", "coess_min_q", "combined_abs_r"],
        title=title, size_max=18,
    )
    fig.update_traces(hovertemplate=(
        "<b>%{customdata[0]}</b><br>co-expr r=%{x:.3f} (q=%{customdata[1]:.1e})"
        "<br>co-ess r=%{y:.3f} (q=%{customdata[2]:.1e})"
        "<br>combined |r|=%{customdata[3]:.3f}<extra></extra>"))
    fig.add_hline(y=0, line_color="#ddd"); fig.add_vline(x=0, line_color="#ddd")
    fig.update_layout(xaxis_title="co-expression mean r", yaxis_title="co-essentiality mean r",
                      legend_title="", height=560)
    return fig


def drilldown_scatter(x: pd.Series, y: pd.Series, labels: pd.Series,
                      xname: str, yname: str, modality: str) -> go.Figure:
    """Per-cell-line scatter of a candidate vs a seed gene."""
    df = pd.DataFrame({"x": x, "y": y, "cell_line": labels})
    fig = px.scatter(df, x="x", y="y", custom_data=["cell_line"], trendline="ols",
                     title=f"{yname} vs {xname}  ({modality})")
    fig.update_traces(hovertemplate="%{customdata[0]}<br>"+xname+"=%{x:.2f}<br>"+yname+"=%{y:.2f}<extra></extra>",
                      marker=dict(size=8, opacity=0.75))
    fig.update_layout(xaxis_title=xname, yaxis_title=yname, height=460)
    return fig


def essentiality_scatter(ranked: pd.DataFrame, title: str,
                         highlight_genes: list[str] | None = None) -> go.Figure:
    """x = % co-expressed, y = % essential; size = joint %, color = common-essential."""
    if ranked is None or ranked.empty:
        return _empty("No ranked candidates")
    df = ranked.reset_index().copy()
    df["ce"] = np.where(df.get("common_essential", False), "common-essential", "selective")
    fig = px.scatter(
        df, x="pct_coexpressed", y="pct_essential", size="pct_joint",
        color="ce", color_discrete_map={"common-essential": "#ff7f0e", "selective": "#1f77b4"},
        custom_data=["gene", "pct_joint", "pct_essential_given_coexpr", "n_joint", "n_lines"],
        title=title, size_max=22,
    )
    fig.update_traces(hovertemplate=(
        "<b>%{customdata[0]}</b><br>co-expressed %{x:.0f}%<br>essential %{y:.0f}%"
        "<br>joint %{customdata[1]:.0f}% (%{customdata[3]}/%{customdata[4]} lines)"
        "<br>essential|co-expr %{customdata[2]:.0f}%<extra></extra>"))
    fig.update_layout(xaxis_title="% cohort lines co-expressed",
                      yaxis_title="% cohort lines essential", legend_title="", height=560)
    return _add_gene_highlights(fig, ranked, highlight_genes, "pct_coexpressed", "pct_essential")


def ranked_bar(ranked: pd.DataFrame, top_k: int, title: str,
               highlight_genes: list[str] | None = None) -> go.Figure:
    """Horizontal bar of top-K genes by selectivity score, colored by common-essential.
    `highlight_genes` (if given) recolors those specific bars crimson, on top of the
    normal common-essential/selective palette -- additive only, existing callers that
    omit it see byte-identical figures."""
    if ranked is None or ranked.empty:
        return _empty("No ranked candidates")
    sub = ranked.head(top_k).iloc[::-1]
    color = np.where(sub.get("common_essential", False), "#ff7f0e", "#2ca02c")
    if highlight_genes:
        color = np.where(sub.index.isin(highlight_genes), "#D11947", color)
    fig = go.Figure()
    fig.add_bar(y=sub.index, x=sub["selectivity_score"], orientation="h",
                marker_color=color,
                customdata=np.stack([sub["pct_joint"], sub["selectivity_delta"],
                                     sub["pct_essential"], sub["pct_essential_pan"]], axis=-1),
                hovertemplate=("<b>%{y}</b><br>selectivity score %{x:.0f}"
                               "<br>joint %{customdata[0]:.0f}%"
                               "<br>selectivity Δ %{customdata[1]:+.0f} pts"
                               "<br>essential here %{customdata[2]:.0f}% vs pan %{customdata[3]:.0f}%"
                               "<extra></extra>"))
    fig.update_layout(title=title, xaxis_title="selectivity-weighted score",
                      height=max(400, 20 * len(sub) + 120))
    return fig


def classifier_bar(lookup_df: pd.DataFrame, title: str) -> go.Figure:
    """Horizontal bar of precomputed protein-classifier probabilities per gene.
    Not-cached genes get a 0-length grey bar (visually distinct, never fabricated)."""
    if lookup_df is None or lookup_df.empty:
        return _empty("No genes to show")
    sub = lookup_df.iloc[::-1]
    colors = np.where(~sub["cached"], "#bbbbbb",
                      np.where(sub["Prediction"] == "INTERACTOR", "#2ca02c", "#d62728"))
    x = sub["Probability"].fillna(0.0)
    labels = np.where(sub["cached"], sub["Prediction"].astype(str), "not generated/cached")
    fig = go.Figure()
    fig.add_bar(y=sub.index, x=x, orientation="h", marker_color=colors,
                customdata=np.stack([labels, sub["Confidence"].fillna("")], axis=-1),
                hovertemplate=("<b>%{y}</b><br>P(interactor) %{x:.3f}"
                              "<br>%{customdata[0]} (confidence: %{customdata[1]})"
                              "<extra></extra>"))
    fig.update_layout(title=title, xaxis_title="P(interactor) — precomputed, one fixed batch",
                      xaxis_range=[0, 1], height=max(360, 26 * len(sub) + 120))
    return fig


def _empty(msg: str) -> go.Figure:
    fig = go.Figure()
    fig.add_annotation(text=msg, showarrow=False, font=dict(size=16, color="#888"))
    fig.update_layout(height=300, xaxis_visible=False, yaxis_visible=False)
    return fig


def string_heatmap(matrix: pd.DataFrame, title: str) -> go.Figure:
    """Candidate x seed STRING score heatmap (0-1)."""
    if matrix is None or matrix.empty:
        return _empty("No STRING scores")
    fig = px.imshow(matrix, color_continuous_scale="Teal", zmin=0, zmax=1, aspect="auto",
                    labels=dict(x="seed", y="candidate", color="STRING"), title=title)
    fig.update_traces(hovertemplate="candidate %{y}<br>seed %{x}<br>score %{z:.3f}<extra></extra>")
    fig.update_layout(height=max(360, 20 * len(matrix) + 140))
    return fig


def string_network(edges: list[dict], node_size: dict, kind: str, title: str) -> go.Figure:
    """Bipartite-ish spring layout: candidates + seeds, edges weighted by STRING.

    kind: 'physical' | 'functional'. node_size maps candidate gene -> selectivity score.
    """
    if not edges:
        return _empty("No STRING interactions above cutoff")
    weighted = [(e["source"], e["target"], e[kind]) for e in edges if e.get(kind, 0) > 0]
    if not weighted:
        return _empty(f"No {kind} interactions above cutoff")
    cands = sorted({s for s, _, _ in weighted})
    seeds = sorted({t for _, t, _ in weighted})
    # two columns: seeds left, candidates right
    pos = {}
    for i, s in enumerate(seeds):
        pos[s] = (0.0, 1.0 - (i / max(1, len(seeds) - 1)) if len(seeds) > 1 else 0.5)
    for i, c in enumerate(cands):
        pos[c] = (1.0, 1.0 - (i / max(1, len(cands) - 1)) if len(cands) > 1 else 0.5)

    edge_traces = []
    for s, t, w in weighted:
        x0, y0 = pos[s]; x1, y1 = pos[t]
        edge_traces.append(go.Scatter(
            x=[x0, x1, None], y=[y0, y1, None], mode="lines",
            line=dict(width=0.5 + 4.5 * w, color=f"rgba(28,114,147,{0.25 + 0.6 * w:.2f})"),
            hoverinfo="text", text=f"{t}–{s}: {w:.3f}", showlegend=False))

    def _node_trace(names, color, sizes, tag):
        return go.Scatter(
            x=[pos[n][0] for n in names], y=[pos[n][1] for n in names],
            mode="markers+text", text=names, textposition="middle right" if tag == "cand" else "middle left",
            marker=dict(size=sizes, color=color, line=dict(width=1, color="white")),
            hovertext=[f"{n}" for n in names], hoverinfo="text", showlegend=False)

    seed_tr = _node_trace(seeds, "#0B2E4F", 16, "seed")
    # candidate size scaled by selectivity score
    svals = [node_size.get(c, 1.0) for c in cands]
    smax = max(svals) or 1.0
    csizes = [14 + 26 * (v / smax) for v in svals]
    cand_tr = _node_trace(cands, "#00A896", csizes, "cand")

    fig = go.Figure(edge_traces + [seed_tr, cand_tr])
    fig.update_layout(title=f"{title} ({kind})", height=max(420, 26 * max(len(cands), len(seeds)) + 120),
                      xaxis=dict(visible=False, range=[-0.3, 1.4]), yaxis=dict(visible=False),
                      margin=dict(l=10, r=10, t=48, b=10))
    return fig


def tissue_heatmap(heat: pd.DataFrame, title: str) -> go.Figure:
    """Genes x (normal tissues + TCGA types + target) log2(TPM+1) heatmap."""
    if heat is None or heat.empty:
        return _empty("No tissue-specificity data")
    log2 = np.log2(heat.astype(float).fillna(0) + 1)
    fig = px.imshow(log2, color_continuous_scale="Viridis", aspect="auto",
                    labels=dict(x="tissue / cancer type", y="gene", color="log2(TPM+1)"),
                    title=title)
    fig.update_traces(hovertemplate="gene %{y}<br>%{x}<br>log2(TPM+1)=%{z:.2f}<extra></extra>")
    fig.update_layout(height=max(380, 20 * len(log2) + 160), xaxis_tickangle=-45)
    return fig


def drug_target_heatmap(df: pd.DataFrame, max_drugs_per_gene: int = 15,
                        vmax: float = 5, title: str = "") -> go.Figure:
    """Genes (rows) x Drugs (columns) heatmap, color = score.

    df must have columns Gene, Drug_Name, Score. Keeps the top-scoring
    max_drugs_per_gene drugs per gene (by max Score across genes) to stay readable.
    """
    if df is None or df.empty:
        return _empty("No drug-target results")
    top_drugs = (df.groupby("Drug_Name")["Score"].max()
                .sort_values(ascending=False).head(max_drugs_per_gene * df["Gene"].nunique()).index)
    sub = df[df["Drug_Name"].isin(top_drugs)]
    mat = sub.pivot_table(index="Gene", columns="Drug_Name", values="Score", aggfunc="max")
    fig = px.imshow(mat, color_continuous_scale="RdYlGn", zmin=0, zmax=vmax, aspect="auto",
                    labels=dict(x="drug", y="gene", color="score"), title=title)
    fig.update_traces(hovertemplate="gene %{y}<br>drug %{x}<br>score %{z}<extra></extra>")
    fig.update_layout(height=max(380, 24 * len(mat) + 160), xaxis_tickangle=-45)
    return fig


def drug_target_score_bar(df: pd.DataFrame, title: str = "") -> go.Figure:
    """Stacked bar: per-gene drug count, stacked by score tier (1-5)."""
    if df is None or df.empty:
        return _empty("No drug-target results")
    counts = df.groupby(["Gene", "Score"]).size().unstack(fill_value=0)
    colors = {5: "#00A896", 4: "#4FB0C6", 3: "#F2C14E", 2: "#E29578", 1: "#C7C7C7"}
    fig = go.Figure()
    for score in sorted(counts.columns, reverse=True):
        fig.add_bar(y=counts.index, x=counts[score], name=f"score {score}", orientation="h",
                   marker_color=colors.get(score, "#888"),
                   hovertemplate=f"score {score}: " + "%{x} drugs<br>%{y}<extra></extra>")
    fig.update_layout(barmode="stack", title=title, xaxis_title="drug count",
                      height=max(380, 24 * len(counts) + 160), legend_title="")
    return fig
