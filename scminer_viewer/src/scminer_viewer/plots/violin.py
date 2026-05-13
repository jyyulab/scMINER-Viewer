from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from ..data import Study
from ._common import cells_mask, cluster_color_map, empty_figure


def violin_plot(
    study: Study,
    gene: str,
    relationship: str = "Express_normalized",
    active_clusters: Optional[Iterable[str]] = None,
) -> go.Figure:
    vals = study.gene_values(gene, relationship)
    if vals is None:
        return empty_figure(f"No {relationship} data for {gene}")

    mask = cells_mask(study, active_clusters)
    cell_types = study.cells["cellType"].to_numpy()[mask]
    values = vals[mask]
    colors = cluster_color_map(study)

    df = pd.DataFrame({"cellType": cell_types, "value": values})
    fig = go.Figure()
    for ct in df["cellType"].unique():
        sub = df.loc[df["cellType"] == ct, "value"]
        fig.add_trace(
            go.Violin(
                y=sub,
                name=ct,
                box_visible=True,
                meanline_visible=True,
                points="outliers",
                line_color=colors.get(ct, "#888888"),
                fillcolor=colors.get(ct, "#888888"),
                opacity=0.6,
            )
        )
    fig.update_layout(
        yaxis=dict(title=gene),
        xaxis=dict(title=""),
        showlegend=False,
        margin=dict(l=50, r=20, t=30, b=80),
    )
    return fig
