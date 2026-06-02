from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import plotly.graph_objects as go

from ..data import Study
from ._common import cells_mask, empty_figure


def feature_plot(
    study: Study,
    gene: str,
    relationship: str = "Express_normalized",
    active_clusters: Optional[Iterable[str]] = None,
    dot_size: float = 4,
    cell_mask: Optional[np.ndarray] = None,
) -> go.Figure:
    vals = study.gene_values(gene, relationship)
    if vals is None:
        return empty_figure(f"No {relationship} data for {gene}")

    mask = cells_mask(study, active_clusters, cell_mask)
    cells = study.cells.loc[mask]
    v = vals[mask]

    fig = go.Figure(
        data=[
            go.Scattergl(
                x=cells["coord1"],
                y=cells["coord2"],
                mode="markers",
                marker=dict(
                    size=dot_size,
                    color=v,
                    colorscale=[(0, "#dddddd"), (0.5, "#7c9fd1"), (1, "#1f3a72")],
                    cmin=float(np.nanmin(v)) if v.size else 0,
                    cmax=float(np.nanmax(v)) if v.size else 1,
                    colorbar=dict(title=dict(text=gene)),
                    opacity=0.85,
                ),
                hovertemplate=(
                    "<b>%{customdata}</b><br>"
                    f"{gene}: %{{marker.color:.3f}}<extra></extra>"
                ),
                customdata=cells.index.to_numpy(),
            )
        ]
    )
    fig.update_layout(
        xaxis=dict(title=f"{study.meta.coordinate}_1", zeroline=False),
        yaxis=dict(title=f"{study.meta.coordinate}_2", zeroline=False),
        margin=dict(l=50, r=20, t=30, b=50),
        autosize=True,
        # Opt out of the equal-aspect lock (see ASPECT_CORRECT_JS) so the
        # feature plot stretches to fill the full panel, like the cluster
        # plot.
        meta=dict(scvFill=True),
    )
    return fig
