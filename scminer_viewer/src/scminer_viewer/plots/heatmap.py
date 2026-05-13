from __future__ import annotations

from typing import Iterable, Optional

import plotly.graph_objects as go

from ..data import Study
from ._common import aggregate_by_cluster, empty_figure


def heatmap_plot(
    study: Study,
    genes: list[str],
    relationship: str = "Express_normalized",
    active_clusters: Optional[Iterable[str]] = None,
) -> go.Figure:
    if not genes:
        return empty_figure("Add gene(s) to build heatmap")
    agg = aggregate_by_cluster(
        study, list(genes), relationship,
        list(active_clusters) if active_clusters is not None else None,
    )
    if agg is None:
        return empty_figure("No data available for the selected genes")
    means, _ = agg
    fig = go.Figure(
        data=go.Heatmap(
            z=means.values,
            x=list(means.columns),
            y=list(means.index),
            colorscale=[(0, "#f7f7f7"), (0.5, "#7c9fd1"), (1, "#1f3a72")],
            colorbar=dict(title=dict(text="mean")),
        )
    )
    fig.update_layout(
        xaxis=dict(title="Cluster"),
        yaxis=dict(title="Gene"),
        margin=dict(l=80, r=20, t=30, b=80),
    )
    return fig
