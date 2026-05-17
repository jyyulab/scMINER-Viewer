"""Write a self-contained scMINER study metadata bundle (HDF5).

Python port of scminerViewer/R/write_bundle.R. The bundle is a small
`.scminer.h5` file containing study metadata, cell + cluster info, the
master gene list, networks, and per-matrix **gene indexes** (which
genes are available in `expression` / `activity_tf` / `activity_sig`).

Expression and activity matrix values are *not* stored in the bundle —
they are read on demand from the on-disk shard tree by
`Study.gene_values()` at runtime.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Optional, Sequence

import h5py
import numpy as np
import pandas as pd


_BUNDLE_VERSION = 1

_META_REQUIRED = (
    "studyID", "studyAbbr", "longTitle", "shortTitle",
    "species", "coordinate",
)
_CELLS_REQUIRED = ("cellID", "cellType", "coord1", "coord2")
_CLUSTERS_REQUIRED = ("cellType", "count")
_NETWORK_REQUIRED = (
    "source", "target", "cellType",
    "mi", "pearson", "spearman", "rho", "pvalue",
)


def write_bundle(
    bundle_path: str | Path,
    meta: Mapping[str, object],
    cells: pd.DataFrame,
    clusters: pd.DataFrame,
    genes: Sequence[str],
    expression_genes: Optional[Sequence[str]] = None,
    activity_tf_genes: Optional[Sequence[str]] = None,
    activity_sig_genes: Optional[Sequence[str]] = None,
    default_genes: Optional[Sequence[str]] = None,
    network_tf: Optional[pd.DataFrame] = None,
    network_sig: Optional[pd.DataFrame] = None,
    overwrite: bool = False,
) -> Path:
    """Write a `.scminer.h5` bundle.

    Args:
        bundle_path: Output path; should end in `.scminer.h5`.
        meta: Mapping with keys `studyID`, `studyAbbr`, `longTitle`,
            `shortTitle`, `species`, `coordinate`.
        cells: DataFrame with columns `cellID`, `cellType`, `cellGroup`
            (optional, defaults to `cellType`), `coord1`, `coord2`.
        clusters: DataFrame with columns `cellType`, `count`, `color`
            (optional), `label_1` (optional), `label_2` (optional).
        genes: Master gene-symbol sequence.
        expression_genes, activity_tf_genes, activity_sig_genes: Optional
            sequences listing genes that have shards in each matrix.
        default_genes: Optional sequence; if present, the app pre-selects
            these on startup.
        network_tf, network_sig: Optional DataFrames with the canonical
            columns.
        overwrite: Allow overwriting `bundle_path`.

    Returns:
        Resolved `Path` to the written bundle.
    """
    bundle_path = Path(bundle_path)
    if bundle_path.exists():
        if not overwrite:
            raise FileExistsError(
                f"Bundle file already exists: {bundle_path} "
                "(pass overwrite=True to replace)"
            )
        bundle_path.unlink()
    bundle_path.parent.mkdir(parents=True, exist_ok=True)

    _validate_meta(meta)
    _validate_cells(cells)
    _validate_clusters(clusters)
    genes = list(genes)
    if len(genes) == 0:
        raise ValueError("`genes` must be a non-empty sequence")

    with h5py.File(bundle_path, "w") as f:
        _write_meta(f, meta)
        _write_cells(f, cells)
        _write_clusters(f, clusters)
        _write_genes(f, genes)
        _write_index(f, expression_genes, activity_tf_genes, activity_sig_genes)
        _write_defaults(f, default_genes)
        if network_tf is not None and len(network_tf) > 0:
            _write_network(f, "network_tf", network_tf)
        if network_sig is not None and len(network_sig) > 0:
            _write_network(f, "network_sig", network_sig)

    return bundle_path


# ---------------------------------------------------------------- helpers


def _utf8_array(values) -> np.ndarray:
    """Coerce to a numpy object array of bytes / strings h5py accepts."""
    arr = np.asarray(list(values), dtype=object)
    return np.array([str(v) for v in arr], dtype=object)


def _write_string_dset(grp: h5py.Group, name: str, values) -> None:
    """Variable-length UTF-8 string dataset — matches the R hdf5r writer."""
    arr = _utf8_array(values)
    dt = h5py.string_dtype(encoding="utf-8")
    grp.create_dataset(name, data=arr.astype(dt), dtype=dt)


def _write_meta(f: h5py.File, meta: Mapping[str, object]) -> None:
    grp = f.create_group("meta")
    for k in _META_REQUIRED:
        _write_string_dset(grp, k, [str(meta.get(k, "") or "")])
    grp.create_dataset("bundleVersion", data=np.array([_BUNDLE_VERSION], dtype=np.int64))


def _write_cells(f: h5py.File, cells: pd.DataFrame) -> None:
    grp = f.create_group("cells")
    cell_id = cells["cellID"].astype(str).tolist()
    cell_type = cells["cellType"].astype(str).tolist()
    cell_group = (
        cells["cellGroup"].astype(str).tolist()
        if "cellGroup" in cells.columns
        else cell_type
    )
    _write_string_dset(grp, "cellID", cell_id)
    _write_string_dset(grp, "cellType", cell_type)
    _write_string_dset(grp, "cellGroup", cell_group)
    grp.create_dataset(
        "coord1", data=np.asarray(cells["coord1"], dtype=np.float64)
    )
    grp.create_dataset(
        "coord2", data=np.asarray(cells["coord2"], dtype=np.float64)
    )


def _write_clusters(f: h5py.File, clusters: pd.DataFrame) -> None:
    grp = f.create_group("clusters")
    cell_type = clusters["cellType"].astype(str).tolist()
    _write_string_dset(grp, "cellType", cell_type)
    grp.create_dataset(
        "count", data=np.asarray(clusters["count"], dtype=np.int64)
    )
    color = (
        clusters["color"].astype(str).tolist()
        if "color" in clusters.columns
        else ["#888888"] * len(clusters)
    )
    _write_string_dset(grp, "color", color)
    if "label_1" in clusters.columns:
        grp.create_dataset(
            "label_1", data=np.asarray(clusters["label_1"], dtype=np.float64)
        )
    if "label_2" in clusters.columns:
        grp.create_dataset(
            "label_2", data=np.asarray(clusters["label_2"], dtype=np.float64)
        )


def _write_genes(f: h5py.File, genes: Sequence[str]) -> None:
    grp = f.create_group("genes")
    _write_string_dset(grp, "symbol", [str(g) for g in genes])


def _write_index(
    f: h5py.File,
    exp_g: Optional[Sequence[str]],
    tf_g: Optional[Sequence[str]],
    sig_g: Optional[Sequence[str]],
) -> None:
    if exp_g is None and tf_g is None and sig_g is None:
        return
    grp = f.create_group("index")
    if exp_g is not None:
        _write_string_dset(grp, "expression", [str(g) for g in exp_g])
    if tf_g is not None:
        _write_string_dset(grp, "activity_tf", [str(g) for g in tf_g])
    if sig_g is not None:
        _write_string_dset(grp, "activity_sig", [str(g) for g in sig_g])


def _write_defaults(f: h5py.File, default_genes: Optional[Sequence[str]]) -> None:
    if default_genes is None or len(list(default_genes)) == 0:
        return
    grp = f.create_group("defaults")
    _write_string_dset(grp, "genes", [str(g) for g in default_genes])


def _write_network(f: h5py.File, name: str, df: pd.DataFrame) -> None:
    miss = [c for c in _NETWORK_REQUIRED if c not in df.columns]
    if miss:
        raise ValueError(
            f"network `{name}` missing columns: {', '.join(miss)}"
        )
    grp = f.create_group(name)
    _write_string_dset(grp, "source",   df["source"].astype(str).tolist())
    _write_string_dset(grp, "target",   df["target"].astype(str).tolist())
    _write_string_dset(grp, "cellType", df["cellType"].astype(str).tolist())
    for col in ("mi", "pearson", "spearman", "rho", "pvalue"):
        grp.create_dataset(col, data=np.asarray(df[col], dtype=np.float64))


def _validate_meta(meta: Mapping[str, object]) -> None:
    if not isinstance(meta, Mapping):
        raise TypeError("`meta` must be a Mapping")
    miss = [k for k in _META_REQUIRED if k not in meta]
    if miss:
        raise ValueError(f"`meta` missing fields: {', '.join(miss)}")


def _validate_cells(cells: pd.DataFrame) -> None:
    if not isinstance(cells, pd.DataFrame):
        raise TypeError("`cells` must be a pandas DataFrame")
    miss = [c for c in _CELLS_REQUIRED if c not in cells.columns]
    if miss:
        raise ValueError(f"`cells` missing columns: {', '.join(miss)}")


def _validate_clusters(clusters: pd.DataFrame) -> None:
    if not isinstance(clusters, pd.DataFrame):
        raise TypeError("`clusters` must be a pandas DataFrame")
    miss = [c for c in _CLUSTERS_REQUIRED if c not in clusters.columns]
    if miss:
        raise ValueError(f"`clusters` missing columns: {', '.join(miss)}")
