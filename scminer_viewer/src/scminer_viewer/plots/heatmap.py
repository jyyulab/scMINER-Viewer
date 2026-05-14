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
    cell_mask: Optional["np.ndarray"] = None,
) -> go.Figure:
    if not genes:
        return empty_figure("Add gene(s) to build heatmap")
    agg = aggregate_by_cluster(
        study, list(genes), relationship,
        list(active_clusters) if active_clusters is not None else None,
        cell_mask=cell_mask,
    )
    if agg is None:
        return empty_figure("No data available for the selected genes")
    means, _ = agg
    # Blue-white-red diverging palette (ColorBrewer RdBu reversed). zmid=0
    # centres white at zero so positive values lean red, negative blue.
    fig = go.Figure(
        data=go.Heatmap(
            z=means.values,
            x=list(means.columns),
            y=list(means.index),
            colorscale=[(0, "#2166ac"), (0.5, "#ffffff"), (1, "#b2182b")],
            zmid=0,
            colorbar=dict(title=dict(text="mean")),
        )
    )
    fig.update_layout(
        xaxis=dict(title="Cluster"),
        yaxis=dict(title="Gene"),
        margin=dict(l=80, r=20, t=30, b=80),
    )
    return fig
