from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from ..data import Study
from ._common import cells_mask, cluster_color_map, empty_figure


def _darken(hex_color: str, factor: float = 0.55) -> str:
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return "#{:02x}{:02x}{:02x}".format(
        int(r * factor), int(g * factor), int(b * factor)
    )


def violin_plot(
    study: Study,
    gene: str,
    relationship: str = "Express_normalized",
    active_clusters: Optional[Iterable[str]] = None,
    cell_mask: Optional[np.ndarray] = None,
) -> go.Figure:
    vals = study.gene_values(gene, relationship)
    if vals is None:
        return empty_figure(f"No {relationship} data for {gene}")

    mask = cells_mask(study, active_clusters, cell_mask)
    cell_types = study.cells["cellType"].to_numpy()[mask]
    values = vals[mask]
    colors = cluster_color_map(study)

    df = pd.DataFrame({"cellType": cell_types, "value": values})
    fig = go.Figure()
    for ct in df["cellType"].unique():
        sub = df.loc[df["cellType"] == ct, "value"]
        color = colors.get(ct, "#888888")
        fig.add_trace(
            go.Violin(
                x=[ct] * len(sub),
                y=sub,
                name=ct,
                box=dict(
                    visible=True,
                    width=0.25,
                    fillcolor=color,
                    line=dict(color=_darken(color), width=1),
                ),
                meanline=dict(visible=True, color=_darken(color), width=1),
                points="outliers",
                line=dict(color=_darken(color)),
                fillcolor=color,
                opacity=0.6,
                showlegend=False,
            )
        )
    fig.update_layout(
        yaxis=dict(title=gene),
        xaxis=dict(title=""),
        showlegend=False,
        margin=dict(l=50, r=20, t=30, b=80),
    )
    return fig
