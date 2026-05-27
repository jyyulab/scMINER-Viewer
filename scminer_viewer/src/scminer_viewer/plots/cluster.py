from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import plotly.graph_objects as go

from ..data import Study
from ._common import cells_mask, cluster_color_map


def cluster_plot(
    study: Study,
    active_clusters: Optional[Iterable[str]] = None,
    dot_size: float = 4,
    show_labels: bool = True,
    cell_mask: Optional[np.ndarray] = None,
) -> go.Figure:
    mask = cells_mask(study, active_clusters, cell_mask)
    cells = study.cells.loc[mask]
    colors = cluster_color_map(study)

    fig = go.Figure()
    for cluster_id in cells["cellType"].unique():
        sub = cells.loc[cells["cellType"] == cluster_id]
        fig.add_trace(
            go.Scattergl(
                x=sub["coord1"],
                y=sub["coord2"],
                mode="markers",
                name=cluster_id,
                marker=dict(
                    size=dot_size,
                    color=colors.get(cluster_id, "#888888"),
                    opacity=0.7,
                ),
                hovertemplate=(
                    f"<b>%{{customdata}}</b><br>"
                    f"{study.meta.coordinate}_1: %{{x:.3f}}<br>"
                    f"{study.meta.coordinate}_2: %{{y:.3f}}<extra></extra>"
                ),
                customdata=sub.index.to_numpy(),
            )
        )

    if show_labels and "label_1" in study.clusters.columns:
        visible = (
            study.clusters
            if active_clusters is None
            else study.clusters.loc[study.clusters.index.isin(active_clusters)]
        )
        for cluster_id, row in visible.iterrows():
            fig.add_annotation(
                x=row["label_1"], y=row["label_2"],
                text=cluster_id, showarrow=False,
                font=dict(size=14, color="#222"),
                bgcolor="rgba(255,255,255,0.7)",
                bordercolor="#666", borderpad=3,
            )

    fig.update_layout(
        xaxis=dict(title=f"{study.meta.coordinate}_1", zeroline=False),
        yaxis=dict(title=f"{study.meta.coordinate}_2", zeroline=False),
        legend=dict(title=dict(text="Cluster")),
        margin=dict(l=50, r=20, t=30, b=50),
    )
    return fig
