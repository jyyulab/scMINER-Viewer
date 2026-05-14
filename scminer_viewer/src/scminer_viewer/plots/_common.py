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


def cells_mask(
    study: Study,
    active_clusters: Optional[Iterable[str]],
    cell_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Intersect the cluster filter with an optional per-cell mask.

    `cell_mask` (boolean ndarray of length n_cells) is typically supplied
    by the downsampling reactive.
    """
    if active_clusters is None:
        base = np.ones(study.n_cells, dtype=bool)
    else:
        active = set(active_clusters)
        base = study.cells["cellType"].isin(active).to_numpy()
    if cell_mask is not None:
        base = base & np.asarray(cell_mask, dtype=bool)
    return base


def aggregate_by_cluster(
    study: Study,
    genes: list[str],
    relationship: str,
    active_clusters: Optional[list[str]] = None,
    cell_mask: Optional[np.ndarray] = None,
) -> Optional[tuple[pd.DataFrame, pd.DataFrame]]:
    """Return (mean, pct_expressing) per gene per cluster (lazy reads).

    Each gene's row is loaded from disk via `study.gene_values(...)`.
    Optional `cell_mask` restricts the per-cluster aggregation to the
    sampled subset of cells (used by the Sampling % control).
    """
    if not genes:
        return None
    present_genes: list[str] = []
    rows: list[np.ndarray] = []
    for g in genes:
        vals = study.gene_values(g, relationship)
        if vals is None:
            continue
        present_genes.append(g)
        rows.append(vals)
    if not present_genes:
        return None

    clusters = active_clusters or list(study.clusters.index)
    cell_types = study.cells["cellType"].to_numpy()
    sample = (np.asarray(cell_mask, dtype=bool)
              if cell_mask is not None
              else np.ones(len(cell_types), dtype=bool))
    means = np.zeros((len(present_genes), len(clusters)), dtype=np.float64)
    pcts = np.zeros_like(means)
    for i, vals in enumerate(rows):
        for j, cluster in enumerate(clusters):
            mask = (cell_types == cluster) & sample
            sub = vals[mask]
            sub = sub[~np.isnan(sub)]
            n = sub.size
            if n == 0:
                continue
            means[i, j] = float(sub.mean())
            pcts[i, j] = float((sub != 0).sum()) / n
    return (
        pd.DataFrame(means, index=present_genes, columns=clusters),
        pd.DataFrame(pcts, index=present_genes, columns=clusters),
    )


def empty_figure(title: str = ""):
    import plotly.graph_objects as go
    fig = go.Figure()
    if title:
        fig.update_layout(title=title)
    fig.update_layout(margin=dict(l=20, r=20, t=40, b=20))
    return fig
