from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ..data import Study
from ._common import cells_mask, cluster_color_map, empty_figure


def heatmap_plot(
    study: Study,
    genes: list[str],
    relationship: str = "Express_normalized",
    active_clusters: Optional[Iterable[str]] = None,
    cell_mask: Optional[np.ndarray] = None,
) -> go.Figure:
    if not genes:
        return empty_figure("Add gene(s) to build heatmap")

    base_mask = cells_mask(study, active_clusters, cell_mask)
    if not base_mask.any():
        return empty_figure("No cells in current selection")

    cell_types = study.cells["cellType"].to_numpy()
    cell_ids = study.cells.index.to_numpy()

    if active_clusters is not None:
        cluster_order = list(active_clusters)
    else:
        cluster_order = list(study.clusters.index)
    present_in_selection = set(cell_types[base_mask])
    ordered_clusters = [c for c in cluster_order if c in present_in_selection]

    keep_idx = np.where(base_mask)[0]
    ct_keep = cell_types[keep_idx]
    order = np.concatenate(
        [keep_idx[ct_keep == c] for c in ordered_clusters]
    ) if ordered_clusters else keep_idx

    present_genes: list[str] = []
    rows: list[np.ndarray] = []
    for g in genes:
        vals = study.gene_values(g, relationship)
        if vals is None:
            continue
        present_genes.append(g)
        rows.append(vals[order])
    if not present_genes:
        return empty_figure("No data available for the selected genes")

    z = np.vstack(rows)
    sorted_cell_ids = cell_ids[order]
    sorted_cell_types = cell_types[order]

    color_map = cluster_color_map(study)
    cluster_to_idx = {c: i for i, c in enumerate(ordered_clusters)}
    ann_z = np.array([[cluster_to_idx[c] for c in sorted_cell_types]])

    n_clusters = len(ordered_clusters)
    if n_clusters <= 1:
        only = ordered_clusters[0] if ordered_clusters else "#cccccc"
        color = color_map.get(only, "#cccccc") if ordered_clusters else "#cccccc"
        ann_colorscale = [(0.0, color), (1.0, color)]
        zmin_ann, zmax_ann = -0.5, 0.5
    else:
        ann_colorscale = []
        for i, c in enumerate(ordered_clusters):
            lo = i / n_clusters
            hi = (i + 1) / n_clusters
            ann_colorscale.append((lo, color_map[c]))
            ann_colorscale.append((hi, color_map[c]))
        zmin_ann, zmax_ann = -0.5, n_clusters - 0.5

    fig = make_subplots(
        rows=2, cols=1,
        row_heights=[0.94, 0.06],
        vertical_spacing=0.02,
        shared_xaxes=True,
    )

    fig.add_trace(
        go.Heatmap(
            z=z,
            x=sorted_cell_ids,
            y=present_genes,
            colorscale=[(0, "#2166ac"), (0.5, "#ffffff"), (1, "#b2182b")],
            zmid=0,
            colorbar=dict(title=dict(text="value")),
            customdata=np.tile(sorted_cell_types, (len(present_genes), 1)),
            hovertemplate=(
                "gene: %{y}<br>cell: %{x}<br>"
                "cellType: %{customdata}<br>value: %{z:.3f}<extra></extra>"
            ),
        ),
        row=1, col=1,
    )

    fig.add_trace(
        go.Heatmap(
            z=ann_z,
            x=sorted_cell_ids,
            y=["cellType"],
            colorscale=ann_colorscale,
            zmin=zmin_ann,
            zmax=zmax_ann,
            showscale=False,
            customdata=np.array([sorted_cell_types]),
            hovertemplate=(
                "cell: %{x}<br>cellType: %{customdata}<extra></extra>"
            ),
        ),
        row=2, col=1,
    )

    for c in ordered_clusters:
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None],
                mode="markers",
                marker=dict(size=10, color=color_map[c],
                            line=dict(width=0)),
                name=c, showlegend=True, hoverinfo="skip",
            ),
            row=1, col=1,
        )

    fig.update_xaxes(showticklabels=False, ticks="", row=1, col=1)
    fig.update_xaxes(showticklabels=False, ticks="",
                     title=dict(text="Cell"), row=2, col=1)
    fig.update_yaxes(title=dict(text="Gene"), autorange="reversed",
                     row=1, col=1)
    fig.update_yaxes(showticklabels=False, ticks="", row=2, col=1)

    fig.update_layout(
        margin=dict(l=80, r=20, t=30, b=100),
        legend=dict(
            title=dict(text="Cell type"),
            orientation="h",
            x=0.5, xanchor="center",
            y=-0.18, yanchor="top",
        ),
    )
    return fig
