from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import plotly.graph_objects as go

from ..data import Study
from ._common import empty_figure

try:
    import networkx as nx
    _HAS_NX = True
except ImportError:  # pragma: no cover
    _HAS_NX = False


def network_plot(
    study: Study,
    gene: str,
    network_type: str = "TF",
    active_clusters: Optional[Iterable[str]] = None,
    max_edges: int = 60,
) -> go.Figure:
    edges_df = (
        study.network_tf if network_type == "TF" else study.network_sig
    )
    if edges_df is None or len(edges_df) == 0:
        return empty_figure(f"No {network_type} network in this study")

    mask = (edges_df["source"] == gene) | (edges_df["target"] == gene)
    if active_clusters is not None:
        mask &= edges_df["cellType"].isin(list(active_clusters))
    sub = edges_df.loc[mask]
    if len(sub) == 0:
        return empty_figure(f"No {network_type} edges for {gene}")

    sub = sub.reindex(sub["mi"].abs().sort_values(ascending=False).index)
    sub = sub.head(max_edges)

    nodes = list({*sub["source"], *sub["target"]})
    if not _HAS_NX:
        # Fall back to a deterministic circular layout
        n = len(nodes)
        angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
        pos = {n_: (float(np.cos(a)), float(np.sin(a)))
               for n_, a in zip(nodes, angles)}
    else:
        g = nx.Graph()
        g.add_nodes_from(nodes)
        for _, r in sub.iterrows():
            g.add_edge(r["source"], r["target"], weight=abs(r["mi"]))
        pos = nx.spring_layout(g, seed=42, weight="weight")

    edge_x, edge_y, edge_colors_pos, edge_colors_neg = [], [], [], []
    for _, r in sub.iterrows():
        x0, y0 = pos[r["source"]]
        x1, y1 = pos[r["target"]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y, mode="lines",
        line=dict(width=0.6, color="#bbbbbb"),
        hoverinfo="skip", showlegend=False,
    )
    node_trace = go.Scatter(
        x=[pos[n_][0] for n_ in nodes],
        y=[pos[n_][1] for n_ in nodes],
        mode="markers+text",
        text=nodes,
        textposition="top center",
        marker=dict(
            size=[32 if n_ == gene else 18 for n_ in nodes],
            color=["#e15759" if n_ == gene else "#4e79a7" for n_ in nodes],
            line=dict(width=1, color="#fff"),
        ),
        hovertemplate="%{text}<extra></extra>",
        showlegend=False,
    )
    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        xaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
        yaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
        margin=dict(l=10, r=10, t=30, b=10),
        title=f"{network_type} network for {gene} ({len(sub)} edges)",
    )
    return fig
