"""Helpers shared across plot modules."""

from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import pandas as pd

from ..data import Study

DEFAULT_COLORS = [
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
]


def cluster_color_map(study: Study) -> dict[str, str]:
    colors = list(study.clusters["color"])
    if any(c is None or (isinstance(c, str) and not c) for c in colors):
        colors = [DEFAULT_COLORS[i % len(DEFAULT_COLORS)]
                  for i in range(len(colors))]
    return dict(zip(study.clusters.index, colors))


def cells_mask(study: Study,
               active_clusters: Optional[Iterable[str]]) -> np.ndarray:
    if active_clusters is None:
        return np.ones(study.n_cells, dtype=bool)
    active = set(active_clusters)
    return study.cells["cellType"].isin(active).to_numpy()


def aggregate_by_cluster(
    study: Study,
    genes: list[str],
    relationship: str,
    active_clusters: Optional[list[str]] = None,
) -> Optional[tuple[pd.DataFrame, pd.DataFrame]]:
    """Return (mean, pct_expressing) per gene per cluster.

    Both DataFrames are indexed by gene with cluster columns.
    """
    mat = {
        "Express_normalized": study.expression,
        "Activity_tf": study.activity_tf,
        "Activity_sig": study.activity_sig,
    }.get(relationship)
    if mat is None or not genes:
        return None
    rows = []
    for g in genes:
        idx = study._gene_to_row.get(g)
        if idx is not None:
            rows.append((g, idx))
    if not rows:
        return None
    gene_names = [g for g, _ in rows]
    row_idx = np.array([i for _, i in rows], dtype=np.int64)

    clusters = active_clusters or list(study.clusters.index)
    cell_types = study.cells["cellType"].to_numpy()

    means = np.zeros((len(rows), len(clusters)), dtype=np.float64)
    pcts = np.zeros_like(means)
    for j, cluster in enumerate(clusters):
        mask = cell_types == cluster
        n = int(mask.sum())
        if n == 0:
            continue
        sub = mat[row_idx][:, mask]
        # sub is a CSR matrix
        means[:, j] = np.asarray(sub.mean(axis=1)).ravel()
        # pct expressing = fraction of cells with non-zero value
        nnz_per_row = np.asarray((sub != 0).sum(axis=1)).ravel()
        pcts[:, j] = nnz_per_row / n
    return (
        pd.DataFrame(means, index=gene_names, columns=clusters),
        pd.DataFrame(pcts, index=gene_names, columns=clusters),
    )


def empty_figure(title: str = ""):
    import plotly.graph_objects as go
    fig = go.Figure()
    if title:
        fig.update_layout(title=title)
    fig.update_layout(margin=dict(l=20, r=20, t=40, b=20))
    return fig
