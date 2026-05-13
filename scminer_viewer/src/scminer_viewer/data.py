"""Load a `.scminer.h5` bundle into a Study dataclass.

The bundle (v2) stores only metadata + per-matrix gene **indexes** (not
the values). Expression / activity values are read lazily from the
on-disk shard tree via `Study.gene_values(...)`.

Bundle layout is documented in `IMPLEMENTATION.md`. The bundle is
expected to be co-located with the shard tree by default:
`<shard_dir>/<studyID>.scminer.h5` and
`<shard_dir>/{expression_files,activity_files}/<studyID>/...`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import gzip
import h5py
import numpy as np
import pandas as pd


_RELATIONSHIPS = {
    "Express_normalized": ("expression_files", "expression_index", "exp"),
    "Activity_tf":        ("activity_files",   "activity_tf_index", "tf"),
    "Activity_sig":       ("activity_files",   "activity_sig_index", "sig"),
}


def _shard_letter(gene: str) -> str:
    first = gene[:1].lower()
    return first if "a" <= first <= "z" else "nm"


def _read_scalar(dset: h5py.Dataset):
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


def _read_network(grp: h5py.Group) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source":   _read_string_array(grp["source"]),
            "target":   _read_string_array(grp["target"]),
            "cellType": _read_string_array(grp["cellType"]),
            "mi":       np.asarray(grp["mi"][:],       dtype=np.float64),
            "pearson":  np.asarray(grp["pearson"][:],  dtype=np.float64),
            "spearman": np.asarray(grp["spearman"][:], dtype=np.float64),
            "rho":      np.asarray(grp["rho"][:],      dtype=np.float64),
            "pvalue":   np.asarray(grp["pvalue"][:],   dtype=np.float64),
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
        meta:               Study metadata.
        cells:              DataFrame indexed by cellID with columns
                            cellType, cellGroup, coord1, coord2.
        clusters:           DataFrame indexed by cellType with columns
                            count, color, label_1?, label_2?.
        genes:              1-D ndarray of master gene symbols (the
                            full picker list).
        expression_index:   ndarray of genes with expression shards
                            (or None if the bundle has no expression).
        activity_tf_index:  ndarray of genes with TF activity shards.
        activity_sig_index: ndarray of genes with SIG activity shards.
        default_genes:      ndarray of genes the app should auto-load,
                            or None.
        network_tf:         Optional DataFrame of TF edges.
        network_sig:        Optional DataFrame of SIG edges.
        bundle_path:        Path the bundle was loaded from.
        shard_dir:          Directory containing the shard tree (defaults
                            to the bundle's parent dir).
    """

    meta: Meta
    cells: pd.DataFrame
    clusters: pd.DataFrame
    genes: np.ndarray
    expression_index: Optional[np.ndarray] = None
    activity_tf_index: Optional[np.ndarray] = None
    activity_sig_index: Optional[np.ndarray] = None
    default_genes: Optional[np.ndarray] = None
    network_tf: Optional[pd.DataFrame] = None
    network_sig: Optional[pd.DataFrame] = None
    bundle_path: Optional[str] = None
    shard_dir: Optional[str] = None
    _cache: dict = field(default_factory=dict, repr=False)
    _index_sets: dict = field(default_factory=dict, repr=False)
    _meta_perm: dict = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        self._index_sets = {
            "Express_normalized": (
                set(self.expression_index) if self.expression_index is not None else None
            ),
            "Activity_tf": (
                set(self.activity_tf_index) if self.activity_tf_index is not None else None
            ),
            "Activity_sig": (
                set(self.activity_sig_index) if self.activity_sig_index is not None else None
            ),
        }

    @property
    def n_cells(self) -> int:
        return len(self.cells)

    @property
    def n_genes(self) -> int:
        return int(len(self.genes))

    @property
    def cell_types(self) -> list[str]:
        return list(self.clusters.index)

    def has_gene(self, gene: str, relationship: str = "Express_normalized") -> bool:
        idx = self._index_sets.get(relationship)
        return idx is not None and gene in idx

    def gene_values(
        self, gene: str, relationship: str = "Express_normalized"
    ) -> Optional[np.ndarray]:
        """Lazily read one gene's values from the shard tree.

        Returns a numpy array of length `n_cells`, aligned to
        `study.cells.index`; cells absent from the shard meta become
        NaN. Returns None if the gene isn't in the requested index or
        the shard file is missing.
        """
        if relationship not in _RELATIONSHIPS:
            return None
        idx = self._index_sets.get(relationship)
        if idx is None or gene not in idx:
            return None

        cache_key = (relationship, gene)
        cached = self._cache.get(cache_key)
        if cached is not None:
            return cached

        perm = self._get_perm(relationship)
        if perm is None:
            return None
        perm_idx, n_shard = perm

        shard_path = self._shard_path(gene, relationship)
        if not shard_path.exists():
            return None

        try:
            with gzip.open(shard_path, "rt") as fh:
                line = fh.readline().strip()
        except OSError:
            return None
        if not line:
            return None
        try:
            shard_vals = np.fromstring(line, sep=",", dtype=np.float64)
        except ValueError:
            return None
        if shard_vals.size != n_shard:
            return None
        aligned = np.full(self.n_cells, np.nan, dtype=np.float64)
        valid = perm_idx >= 0
        aligned[valid] = shard_vals[perm_idx[valid]]
        self._cache[cache_key] = aligned
        return aligned

    def _shard_path(self, gene: str, relationship: str) -> Path:
        subdir, _, _ = _RELATIONSHIPS[relationship]
        kind = ""
        if relationship == "Activity_tf":
            kind = "TF"
        elif relationship == "Activity_sig":
            kind = "SIG"
        name = gene.replace("/", "_")
        letter = _shard_letter(gene)
        base = Path(self.shard_dir) / subdir / self.meta.studyID
        if kind:
            base = base / kind
        return base / letter / f"{name}.csv.gz"

    def _get_perm(self, relationship: str):
        """Return (perm_idx, n_shard_cells) for the relationship's meta.

        perm_idx[i] is the column index in the shard for cell `i` in
        `study.cells`, or -1 if that cell isn't represented.
        """
        tag = _RELATIONSHIPS[relationship][2]
        if tag in self._meta_perm:
            return self._meta_perm[tag]
        subdir, _, _ = _RELATIONSHIPS[relationship]
        meta_path = Path(self.shard_dir) / subdir / self.meta.studyID / "meta.csv"
        if not meta_path.exists():
            self._meta_perm[tag] = None
            return None
        with open(meta_path) as fh:
            header = fh.readline().strip()
        if not header:
            self._meta_perm[tag] = None
            return None
        shard_cells = [s.strip() for s in header.split(",")]
        # Build index from shard cell to position
        shard_pos = {c: i for i, c in enumerate(shard_cells)}
        perm_idx = np.fromiter(
            (shard_pos.get(c, -1) for c in self.cells.index),
            dtype=np.int64,
            count=self.n_cells,
        )
        result = (perm_idx, len(shard_cells))
        self._meta_perm[tag] = result
        return result

    def __repr__(self) -> str:
        def count(v) -> str:
            return "-" if v is None else f"{len(v)}"

        return (
            f"<Study {self.meta.studyAbbr} ({self.meta.studyID}) "
            f"cells={self.n_cells} genes={self.n_genes} "
            f"clusters={len(self.clusters)} "
            f"expression_index={count(self.expression_index)} "
            f"activity_tf_index={count(self.activity_tf_index)} "
            f"activity_sig_index={count(self.activity_sig_index)} "
            f"defaults={count(self.default_genes)} "
            f"network_tf={'yes' if self.network_tf is not None else '-'} "
            f"network_sig={'yes' if self.network_sig is not None else '-'} "
            f"shard_dir={self.shard_dir}>"
        )


def load_study(bundle_path: str | Path, shard_dir: str | Path | None = None) -> Study:
    """Load a `.scminer.h5` bundle from disk.

    Args:
        bundle_path: Path to the bundle file.
        shard_dir:   Directory containing `expression_files/<studyID>/`
                     and `activity_files/<studyID>/`. Defaults to the
                     bundle's parent directory.
    """
    bundle_path = Path(bundle_path)
    if not bundle_path.exists():
        raise FileNotFoundError(f"Bundle not found: {bundle_path}")
    shard_dir = Path(shard_dir) if shard_dir is not None else bundle_path.parent

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
                "cellType":  _read_string_array(f["cells/cellType"]),
                "cellGroup": _read_string_array(f["cells/cellGroup"]),
                "coord1":    np.asarray(f["cells/coord1"][:], dtype=np.float64),
                "coord2":    np.asarray(f["cells/coord2"][:], dtype=np.float64),
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

        def opt(path):
            if path in f:
                return _read_string_array(f[path])
            return None

        expression_index   = opt("index/expression")
        activity_tf_index  = opt("index/activity_tf")
        activity_sig_index = opt("index/activity_sig")
        default_genes      = opt("defaults/genes")

        network_tf  = _read_network(f["network_tf"])  if "network_tf"  in f else None
        network_sig = _read_network(f["network_sig"]) if "network_sig" in f else None

    return Study(
        meta=meta,
        cells=cells,
        clusters=clusters,
        genes=genes,
        expression_index=expression_index,
        activity_tf_index=activity_tf_index,
        activity_sig_index=activity_sig_index,
        default_genes=default_genes,
        network_tf=network_tf,
        network_sig=network_sig,
        bundle_path=str(bundle_path),
        shard_dir=str(shard_dir),
    )
