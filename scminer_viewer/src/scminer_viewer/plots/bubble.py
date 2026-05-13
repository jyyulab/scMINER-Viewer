from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import plotly.graph_objects as go

from ..data import Study
from ._common import aggregate_by_cluster, empty_figure


def bubble_plot(
    study: Study,
    genes: list[str],
    relationship: str = "Express_normalized",
    active_clusters: Optional[Iterable[str]] = None,
) -> go.Figure:
    if not genes:
        return empty_figure("Add gene(s) to build bubble plot")
    agg = aggregate_by_cluster(
        study, list(genes), relationship,
        list(active_clusters) if active_clusters is not None else None,
    )
    if agg is None:
        return empty_figure("No data available for the selected genes")
    means, pcts = agg

    xs, ys, sizes, colors, hover = [], [], [], [], []
    for gene in means.index:
        for cluster in means.columns:
            xs.append(cluster)
            ys.append(gene)
            sizes.append(float(pcts.at[gene, cluster]))
            colors.append(float(means.at[gene, cluster]))
            hover.append(
                f"<b>{gene}</b> in <b>{cluster}</b><br>"
                f"mean: {means.at[gene, cluster]:.3f}<br>"
                f"pct expressing: {pcts.at[gene, cluster]:.1%}"
            )

    sizes_arr = np.array(sizes)
    max_size = float(sizes_arr.max()) if sizes_arr.size else 1.0
    if max_size <= 0:
        max_size = 1.0
    scaled = 2 + 28 * (sizes_arr / max_size)

    fig = go.Figure(
        data=[
            go.Scatter(
                x=xs, y=ys, mode="markers",
                marker=dict(
                    size=scaled,
                    color=colors,
                    colorscale=[(0, "#f7f7f7"), (0.5, "#7c9fd1"), (1, "#1f3a72")],
                    colorbar=dict(title=dict(text="mean")),
                    line=dict(width=0),
                ),
                text=hover,
                hovertemplate="%{text}<extra></extra>",
            )
        ]
    )
    fig.update_layout(
        xaxis=dict(title="Cluster"),
        yaxis=dict(title="Gene", autorange="reversed"),
        margin=dict(l=80, r=20, t=30, b=80),
    )
    return fig
