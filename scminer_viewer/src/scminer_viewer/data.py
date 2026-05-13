"""Load a `.scminer.h5` bundle into a Study dataclass.

The bundle layout is defined by `scminerViewer::write_bundle()` in R and
documented in `IMPLEMENTATION.md`. This module mirrors that contract:

* String datasets are decoded utf-8.
* Sparse matrices use scipy CSR (zero-based indices) so we can
  reconstruct them with `scipy.sparse.csr_matrix((data, indices, indptr),
  shape=shape)` directly.
* Optional groups are returned as `None` when absent in the bundle.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


def _read_scalar(dset: h5py.Dataset):
    """Read a 1-elem (or 0-d) HDF5 dataset that R wrote as a scalar."""
    val = dset[()]
    if isinstance(val, np.ndarray):
        val = val.item() if val.shape == () else val[0]
    if isinstance(val, bytes):
        val = val.decode("utf-8")
    return val


def _read_string_array(dset: h5py.Dataset) -> np.ndarray:
    arr = dset[:]
    if arr.dtype == object or arr.dtype.kind in ("S", "U"):
        return np.array(
            [s.decode("utf-8") if isinstance(s, (bytes, bytearray)) else s
             for s in arr],
            dtype=object,
        )
    return arr


def _read_sparse(grp: h5py.Group) -> csr_matrix:
    data = np.asarray(grp["data"][:], dtype=np.float64)
    indices = np.asarray(grp["indices"][:], dtype=np.int64)
    indptr = np.asarray(grp["indptr"][:], dtype=np.int64)
    shape = tuple(int(x) for x in grp["shape"][:])
    if len(shape) != 2:
        raise ValueError(f"Expected shape of length 2, got {shape!r}")
    return csr_matrix((data, indices, indptr), shape=shape)


def _read_network(grp: h5py.Group) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source": _read_string_array(grp["source"]),
            "target": _read_string_array(grp["target"]),
            "cellType": _read_string_array(grp["cellType"]),
            "mi": np.asarray(grp["mi"][:], dtype=np.float64),
            "pearson": np.asarray(grp["pearson"][:], dtype=np.float64),
            "spearman": np.asarray(grp["spearman"][:], dtype=np.float64),
            "rho": np.asarray(grp["rho"][:], dtype=np.float64),
            "pvalue": np.asarray(grp["pvalue"][:], dtype=np.float64),
        }
    )


@dataclass
class Meta:
    studyID: str
    studyAbbr: str
    longTitle: str
    shortTitle: str
    species: str
    coordinate: str
    bundleVersion: int


@dataclass
class Study:
    """An scMINER study loaded from an HDF5 bundle.

    Attributes:
        meta:        Study metadata.
        cells:       DataFrame indexed by cellID with columns
                     cellType, cellGroup, coord1, coord2.
        clusters:    DataFrame indexed by cellType with columns
                     count, color, label_1?, label_2?.
        genes:       1-D ndarray of gene symbols (row order of matrices).
        expression:  Optional CSR sparse matrix, shape (G, N).
        activity_tf: Optional CSR sparse matrix, shape (G, N).
        activity_sig: Optional CSR sparse matrix, shape (G, N).
        network_tf:  Optional DataFrame of TF edges.
        network_sig: Optional DataFrame of SIG edges.
        bundle_path: Filesystem path the bundle was loaded from.
    """

    meta: Meta
    cells: pd.DataFrame
    clusters: pd.DataFrame
    genes: np.ndarray
    expression: Optional[csr_matrix] = None
    activity_tf: Optional[csr_matrix] = None
    activity_sig: Optional[csr_matrix] = None
    network_tf: Optional[pd.DataFrame] = None
    network_sig: Optional[pd.DataFrame] = None
    bundle_path: Optional[str] = None
    _gene_to_row: dict = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        self._gene_to_row = {g: i for i, g in enumerate(self.genes)}

    @property
    def n_cells(self) -> int:
        return len(self.cells)

    @property
    def n_genes(self) -> int:
        return int(len(self.genes))

    @property
    def cell_types(self) -> list[str]:
        return list(self.clusters.index)

    def gene_values(
        self, gene: str, relationship: str = "Express_normalized"
    ) -> Optional[np.ndarray]:
        """Return the dense row vector for `gene` from the requested matrix.

        relationship: one of 'Express_normalized', 'Activity_tf',
        'Activity_sig'. Returns None if the matrix is absent or the
        gene is unknown.
        """
        mat = {
            "Express_normalized": self.expression,
            "Activity_tf": self.activity_tf,
            "Activity_sig": self.activity_sig,
        }.get(relationship)
        if mat is None:
            return None
        idx = self._gene_to_row.get(gene)
        if idx is None:
            return None
        row = mat.getrow(idx).toarray().ravel()
        return row

    def __repr__(self) -> str:
        def mark(v) -> str:
            return "yes" if v is not None else "-"

        return (
            f"<Study {self.meta.studyAbbr} ({self.meta.studyID}) "
            f"cells={self.n_cells} genes={self.n_genes} "
            f"clusters={len(self.clusters)} "
            f"expression={mark(self.expression)} "
            f"activity_tf={mark(self.activity_tf)} "
            f"activity_sig={mark(self.activity_sig)} "
            f"network_tf={'yes' if self.network_tf is not None else '-'} "
            f"network_sig={'yes' if self.network_sig is not None else '-'}>"
        )


def load_study(bundle_path: str | Path) -> Study:
    """Load a `.scminer.h5` bundle from disk."""
    bundle_path = Path(bundle_path)
    if not bundle_path.exists():
        raise FileNotFoundError(f"Bundle not found: {bundle_path}")

    with h5py.File(bundle_path, "r") as f:
        meta = Meta(
            studyID=str(_read_scalar(f["meta/studyID"])),
            studyAbbr=str(_read_scalar(f["meta/studyAbbr"])),
            longTitle=str(_read_scalar(f["meta/longTitle"])),
            shortTitle=str(_read_scalar(f["meta/shortTitle"])),
            species=str(_read_scalar(f["meta/species"])),
            coordinate=str(_read_scalar(f["meta/coordinate"])),
            bundleVersion=int(_read_scalar(f["meta/bundleVersion"])),
        )

        cell_ids = _read_string_array(f["cells/cellID"])
        cells = pd.DataFrame(
            {
                "cellType": _read_string_array(f["cells/cellType"]),
                "cellGroup": _read_string_array(f["cells/cellGroup"]),
                "coord1": np.asarray(f["cells/coord1"][:], dtype=np.float64),
                "coord2": np.asarray(f["cells/coord2"][:], dtype=np.float64),
            },
            index=pd.Index(cell_ids, name="cellID"),
        )

        ct = _read_string_array(f["clusters/cellType"])
        cluster_data = {
            "count": np.asarray(f["clusters/count"][:], dtype=np.int64),
            "color": _read_string_array(f["clusters/color"]),
        }
        if "label_1" in f["clusters"]:
            cluster_data["label_1"] = np.asarray(
                f["clusters/label_1"][:], dtype=np.float64
            )
        if "label_2" in f["clusters"]:
            cluster_data["label_2"] = np.asarray(
                f["clusters/label_2"][:], dtype=np.float64
            )
        clusters = pd.DataFrame(cluster_data, index=pd.Index(ct, name="cellType"))

        genes = _read_string_array(f["genes/symbol"])

        expression = _read_sparse(f["expression"]) if "expression" in f else None
        activity_tf = (
            _read_sparse(f["activity_tf"]) if "activity_tf" in f else None
        )
        activity_sig = (
            _read_sparse(f["activity_sig"]) if "activity_sig" in f else None
        )
        network_tf = (
            _read_network(f["network_tf"]) if "network_tf" in f else None
        )
        network_sig = (
            _read_network(f["network_sig"]) if "network_sig" in f else None
        )

    return Study(
        meta=meta,
        cells=cells,
        clusters=clusters,
        genes=genes,
        expression=expression,
        activity_tf=activity_tf,
        activity_sig=activity_sig,
        network_tf=network_tf,
        network_sig=network_sig,
        bundle_path=str(bundle_path),
    )
