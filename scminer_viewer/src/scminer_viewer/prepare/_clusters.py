"""Cluster auto-fill: per-cluster cell counts, palette colors, centroid
labels. Python port of scminerViewer/R/cluster_fill.R.

Used by `prepare_study_data()` so the auto-generated `study_meta.csv`
always carries real counts + colors + label positions, and exposed
publicly so callers can pre-fill clusters before `write_bundle()`.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from ._palettes import cluster_colors, DEFAULT_PALETTE


def fill_clusters(
    cells: pd.DataFrame,
    clusters: Optional[pd.DataFrame] = None,
    palette: str = DEFAULT_PALETTE,
) -> pd.DataFrame:
    """Auto-fill a clusters DataFrame from a cells DataFrame.

    Missing columns are computed from the cells:

    * `count` — cells per `cellType` (from `value_counts`)
    * `color` — colors from a named palette (default `npg`)
    * `label_1`, `label_2` — per-cluster `(coord1, coord2)` centroid

    Existing values are preserved; only NA / empty cells are filled in.

    Args:
        cells: DataFrame with `cellID`, `cellType`, `coord1`, `coord2`.
        clusters: Existing clusters DataFrame, or `None` to build from
            scratch using `unique(cells.cellType)`.
        palette: Palette name (see `_palettes.PALETTES`).

    Returns:
        A DataFrame with columns `cellType`, `count`, `color`,
        `label_1`, `label_2`.
    """
    required = ("cellType", "coord1", "coord2")
    miss = [c for c in required if c not in cells.columns]
    if miss:
        raise ValueError(f"`cells` missing columns: {', '.join(miss)}")

    counts = (
        cells["cellType"]
        .astype(str)
        .value_counts(sort=False)
        .rename_axis("cellType")
        .reset_index(name="count")
    )

    if clusters is None or len(clusters) == 0:
        clusters = pd.DataFrame(
            {"cellType": counts["cellType"], "count": counts["count"]}
        )
    else:
        clusters = clusters.copy()
        if "count" not in clusters.columns or clusters["count"].isna().all():
            mapping = dict(zip(counts["cellType"], counts["count"]))
            clusters["count"] = (
                clusters["cellType"].map(mapping).fillna(0).astype(int)
            )

    # Colors — fill any missing/empty.
    if "color" not in clusters.columns:
        clusters["color"] = cluster_colors(len(clusters), palette=palette)
    else:
        auto = cluster_colors(len(clusters), palette=palette)
        col = clusters["color"].astype("object")
        for i in range(len(clusters)):
            val = col.iloc[i]
            if pd.isna(val) or not str(val):
                col.iloc[i] = auto[i]
        clusters["color"] = col

    # Centroids — fill any missing/NA label_1 / label_2.
    need_labels = (
        "label_1" not in clusters.columns
        or "label_2" not in clusters.columns
        or clusters.get("label_1", pd.Series(dtype=float)).isna().any()
        or clusters.get("label_2", pd.Series(dtype=float)).isna().any()
    )
    if need_labels:
        labels = _compute_cluster_labels(cells, clusters["cellType"].astype(str).tolist())
        for col_name, vals in zip(("label_1", "label_2"), (labels[0], labels[1])):
            if col_name not in clusters.columns:
                clusters[col_name] = vals
            else:
                current = clusters[col_name].to_numpy(dtype=float)
                na = np.isnan(current)
                current[na] = np.asarray(vals)[na]
                clusters[col_name] = current

    return clusters.reset_index(drop=True)


def _compute_cluster_labels(
    cells: pd.DataFrame, cell_types: list[str]
) -> tuple[list[float], list[float]]:
    """Mean coord1/coord2 per cellType, ordered to match `cell_types`."""
    grouped = (
        cells[["cellType", "coord1", "coord2"]]
        .assign(cellType=lambda df: df["cellType"].astype(str))
        .groupby("cellType", sort=False)
        .mean()
    )
    label_1: list[float] = []
    label_2: list[float] = []
    for ct in cell_types:
        if ct in grouped.index:
            label_1.append(float(grouped.loc[ct, "coord1"]))
            label_2.append(float(grouped.loc[ct, "coord2"]))
        else:
            label_1.append(0.0)
            label_2.append(0.0)
    return label_1, label_2
