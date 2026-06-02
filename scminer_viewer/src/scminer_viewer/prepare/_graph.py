"""Write the graph-import directory layout + per-gene shard tree.

Python port of scminerViewer/R/write_graph.R. Translates the original
scMINER-portal-datapre-R `h_metadata.R` / `h_activity.R` /
`h_networks.R` writers; per-row write loops vectorised; no Biobase /
yaml dependency at this layer (those live in `orchestrator.py`).
"""

from __future__ import annotations

import csv
import gzip
import io
from pathlib import Path
from string import ascii_lowercase
from typing import Iterable, Mapping, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp


_GRAPH_SUBDIRS = (
    "Header", "Study", "Gene", "Cell",
    "Study_Contains_Gene", "Study_Contains_Cell",
    "Network_TF_Activity", "Network_SIG_Activity",
    "study_gene_expression", "study_gene_tf",
    "study_gene_sig", "study_meta",
)

_SHARD_LETTERS = tuple(ascii_lowercase) + ("nm",)


def shard_letter(gene: str) -> str:
    """Return the shard-directory letter for `gene` (matches R)."""
    first = gene[:1].lower()
    return first if "a" <= first <= "z" else "nm"


def ensure_graph_tree(out_dir: Path) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for sub in _GRAPH_SUBDIRS:
        (out_dir / sub).mkdir(parents=True, exist_ok=True)


def ensure_shard_tree(out_dir: Path, kind: str) -> Path:
    """`kind` is e.g. ``expression_files/<sid>`` or ``activity_files/<sid>/TF``."""
    base = Path(out_dir) / kind
    base.mkdir(parents=True, exist_ok=True)
    for letter in _SHARD_LETTERS:
        (base / letter).mkdir(parents=True, exist_ok=True)
    return base


# --- Study / clusters -------------------------------------------------------


def write_graph_study(out_dir: Path, meta: Mapping[str, object]) -> None:
    out_dir = Path(out_dir)
    study_id = str(meta["studyID"])
    header = "\t".join([
        "studyID:ID(Study)", "studyAbbr", "longTitle", "shortTitle",
    ])
    (out_dir / "Header" / f"{study_id}_study.header.tsv").write_text(
        header + "\n", encoding="utf-8"
    )
    row = "\t".join([
        str(meta["studyID"]), str(meta["studyAbbr"]),
        str(meta["longTitle"]), str(meta["shortTitle"]),
    ])
    (out_dir / "Study" / f"{study_id}_study.tsv").write_text(
        row + "\n", encoding="utf-8"
    )


def write_graph_clusters(
    out_dir: Path,
    meta: Mapping[str, object],
    clusters: Optional[pd.DataFrame],
) -> None:
    if clusters is None or len(clusters) == 0:
        return
    out_dir = Path(out_dir)
    study_id = str(meta["studyID"])
    study_abbr = str(meta["studyAbbr"])

    df = clusters.copy()
    if "NetworkCellType" not in df.columns:
        df["NetworkCellType"] = df["cellType"]
    if "cellGroup" not in df.columns:
        df["cellGroup"] = df["cellType"]
    if "color" not in df.columns:
        df["color"] = "#888888"
    if "label_1" not in df.columns:
        df["label_1"] = 0.0
    if "label_2" not in df.columns:
        df["label_2"] = 0.0

    out = pd.DataFrame({
        "StudyID":         study_id,
        "StudyAbbr":       study_abbr,
        "CellType":        df["cellType"].astype(str),
        "CellGroup":       df["cellGroup"].astype(str),
        "Color":           df["color"].astype(str),
        "Label_1":         pd.to_numeric(df["label_1"], errors="coerce"),
        "Label_2":         pd.to_numeric(df["label_2"], errors="coerce"),
        "NetworkCellType": df["NetworkCellType"].astype(str),
    })
    path = out_dir / "study_meta" / f"{study_id}_study_meta.csv"
    out.to_csv(path, index=False)


# --- Genes ------------------------------------------------------------------


def write_graph_genes(
    out_dir: Path,
    meta: Mapping[str, object],
    genes: Iterable[str],
) -> None:
    out_dir = Path(out_dir)
    study_id = str(meta["studyID"])
    genes_list = [str(g) for g in genes]

    (out_dir / "Header" / f"{study_id}_n_gene.header.tsv").write_text(
        "geneSymbol:ID(Gene)\n", encoding="utf-8"
    )
    (out_dir / "Gene" / f"{study_id}_n_gene.tsv").write_text(
        "\n".join(genes_list) + "\n", encoding="utf-8"
    )
    (out_dir / "Header" / f"{study_id}_r_study_gene.header.tsv").write_text(
        ":START_ID(Study)\t:END_ID(Gene)\n", encoding="utf-8"
    )
    lines = "\n".join(f"{study_id}\t{g}" for g in genes_list) + "\n"
    (out_dir / "Study_Contains_Gene" / f"{study_id}_r_study_gene.tsv").write_text(
        lines, encoding="utf-8"
    )


# --- Cells ------------------------------------------------------------------


def write_graph_cells(
    out_dir: Path,
    meta: Mapping[str, object],
    cells: pd.DataFrame,
) -> None:
    out_dir = Path(out_dir)
    study_id = str(meta["studyID"])
    coordinate = str(meta.get("coordinate") or "UMAP")

    header_cols = [
        "cellID:ID(Cell)", "cellType", "cellGroup",
        f"{coordinate}_1:float", f"{coordinate}_2:float",
        "coordinateName",
    ]
    (out_dir / "Header" / f"{study_id}_n_cell.header.tsv").write_text(
        "\t".join(header_cols) + "\n", encoding="utf-8"
    )

    cell_id = cells["cellID"].astype(str).tolist()
    cell_type = cells["cellType"].astype(str).tolist()
    cell_group = (
        cells["cellGroup"].astype(str).tolist()
        if "cellGroup" in cells.columns
        else cell_type
    )
    coord1 = cells["coord1"].astype(str).tolist()
    coord2 = cells["coord2"].astype(str).tolist()

    # Preserve the original 7-column form: cellID is written twice — the
    # original sapply-over-row layout, and downstream Java backend expects it.
    lines = [
        "\t".join((cid, cid, ct, cg, c1, c2, coordinate))
        for cid, ct, cg, c1, c2 in zip(cell_id, cell_type, cell_group, coord1, coord2)
    ]
    (out_dir / "Cell" / f"{study_id}_n_cell.tsv").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )

    (out_dir / "Header" / f"{study_id}_r_study_cell.header.tsv").write_text(
        ":START_ID(Study)\t:END_ID(Cell)\n", encoding="utf-8"
    )
    rs = "\n".join(f"{study_id}\t{cid}" for cid in cell_id) + "\n"
    (out_dir / "Study_Contains_Cell" / f"{study_id}_r_study_cell.tsv").write_text(
        rs, encoding="utf-8"
    )


# --- Networks ---------------------------------------------------------------


_NETWORK_HEADER = "\t".join([
    ":START_ID(Gene)", ":END_ID(Gene)",
    "relationshipType:string", "studyID:string", "cellType:string",
    "mi:float", "pearson:float", "spearman:float",
    "rho:float", "pvalue:float",
])


def write_graph_networks(
    out_dir: Path,
    meta: Mapping[str, object],
    network_tf: Optional[pd.DataFrame],
    network_sig: Optional[pd.DataFrame],
) -> None:
    out_dir = Path(out_dir)
    study_id = str(meta["studyID"])
    for kind, df in (("TF", network_tf), ("SIG", network_sig)):
        if df is None or len(df) == 0:
            continue
        sub = "Network_TF_Activity" if kind == "TF" else "Network_SIG_Activity"
        (out_dir / "Header" / f"{study_id}_{kind}.header.tsv").write_text(
            _NETWORK_HEADER + "\n", encoding="utf-8"
        )
        body = "\n".join(
            "\t".join((
                str(r.source), str(r.target),
                kind, study_id, str(r.cellType),
                str(r.mi), str(r.pearson), str(r.spearman),
                str(r.rho), str(r.pvalue),
            ))
            for r in df.itertuples(index=False)
        )
        (out_dir / sub / f"{study_id}_{kind}.tsv").write_text(
            body + "\n", encoding="utf-8"
        )


# --- Per-gene shards (expression + activity) --------------------------------


def write_graph_shards(
    out_dir: Path,
    meta: Mapping[str, object],
    mat,
    *,
    kind: str,
    meta_kind: str,
    manifest_dir: str,
    manifest_name: str,
    cell_ids: Iterable[str],
    type_label: str,
    genes: Iterable[str],
    verbose: bool = False,
    progress: bool = False,
) -> None:
    """Write the per-gene shard tree + manifest CSV + cell-header.

    Args:
        out_dir: Study output root.
        meta: Study metadata (uses studyID, studyAbbr, species).
        mat: 2-D matrix (genes x cells). Accepts numpy ndarray,
            scipy.sparse matrix, or pandas DataFrame.
        kind: Shard subdir, e.g. ``expression_files/<sid>`` or
            ``activity_files/<sid>/TF``.
        meta_kind: Where the cell-header meta.csv lives — for activity,
            one level above `kind` so TF and SIG share one meta.csv.
        manifest_dir: e.g. ``study_gene_expression``.
        manifest_name: e.g. ``expression``.
        cell_ids: Cell IDs in column order (written as meta.csv).
        type_label: Manifest's Type column value (e.g. "Expression").
        genes: Row names matching `mat.shape[0]`.
        verbose: Emit progress every 1000 shards.
    """
    if mat is None:
        return
    out_dir = Path(out_dir)
    study_id = str(meta["studyID"])
    abbr = str(meta["studyAbbr"])
    species = str(meta.get("species") or "")

    shard_base = ensure_shard_tree(out_dir, kind)
    meta_dir = out_dir / meta_kind
    meta_dir.mkdir(parents=True, exist_ok=True)
    meta_csv_rel = f"{meta_kind}/meta.csv"
    (meta_dir / "meta.csv").write_text(
        ",".join(str(c) for c in cell_ids) + "\n", encoding="utf-8"
    )

    genes_list = [str(g) for g in genes]
    if len(genes_list) != _nrows(mat):
        raise ValueError(
            f"genes length ({len(genes_list)}) != mat rows ({_nrows(mat)})"
        )

    manifest_path = out_dir / manifest_dir / f"{study_id}_{manifest_name}.csv"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    with open(manifest_path, "w", encoding="utf-8", newline="") as manifest_fh:
        manifest_writer = csv.writer(manifest_fh)
        manifest_writer.writerow([
            "GeneSymbol", "Species", "StudyID", "StudyAbbr",
            "Type", "FileHeader", "File",
        ])

        # A live progress bar over the per-gene shard loop — the slow part
        # of prepare_study (one gzipped file per gene). Falls back to the
        # every-1000 verbose message when `progress` is off or tqdm is
        # unavailable.
        rows = _progress_rows(
            _iter_rows(mat), len(genes_list),
            desc=f"  {manifest_name}", enabled=progress,
        )
        for i, row_vals in enumerate(rows):
            gene = genes_list[i]
            name = gene.replace("/", "_")
            letter = shard_letter(gene)
            rel = f"{kind}/{letter}/{name}.csv.gz"
            manifest_writer.writerow([
                gene, species, study_id, abbr, type_label,
                meta_csv_rel, rel,
            ])

            csv_gz_path = shard_base / letter / f"{name}.csv.gz"
            # Write a single-row CSV directly to gzip to skip the
            # intermediate uncompressed file the R version writes.
            _write_shard_csv_gz(csv_gz_path, row_vals)

            if not progress and verbose and (i + 1) % 1000 == 0:
                print(
                    f"  [{manifest_name}] {i + 1}/{len(genes_list)} shards written",
                    flush=True,
                )


def _progress_rows(iterable, total: int, desc: str, enabled: bool):
    """Wrap a row iterator in a tqdm progress bar when `enabled`.

    Falls back to the bare iterator if `enabled` is False or tqdm isn't
    installed (the per-1000 verbose print then provides coarse progress).
    """
    if not enabled:
        return iterable
    try:
        from tqdm import tqdm
    except ImportError:
        return iterable
    return tqdm(iterable, total=total, desc=desc, unit="gene", leave=True)


def _nrows(mat) -> int:
    if hasattr(mat, "shape"):
        return int(mat.shape[0])
    if isinstance(mat, pd.DataFrame):
        return mat.shape[0]
    raise TypeError(f"Unsupported matrix type: {type(mat).__name__}")


def _iter_rows(mat):
    """Yield each row of `mat` as a 1-D numpy array of floats."""
    if isinstance(mat, pd.DataFrame):
        for i in range(mat.shape[0]):
            yield np.asarray(mat.iloc[i, :].to_numpy(), dtype=np.float64)
        return
    if sp.issparse(mat):
        # Convert to CSR once for cheap row slicing.
        csr = mat.tocsr()
        n = csr.shape[0]
        for i in range(n):
            row = csr.getrow(i).toarray().ravel()
            yield np.asarray(row, dtype=np.float64)
        return
    arr = np.asarray(mat)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2-D matrix; got shape {arr.shape}")
    for i in range(arr.shape[0]):
        yield np.asarray(arr[i, :], dtype=np.float64)


def _write_shard_csv_gz(path: Path, values: np.ndarray) -> None:
    """Write one row of comma-separated values directly to gzip.

    Matches what R's data.table::fwrite + R.utils::gzip would produce:
    a single line of comma-separated values terminated by `\\n`. We use
    np.savetxt to a StringIO buffer then write through gzip.GzipFile,
    which beats per-element str()/format() concatenation on large rows.
    """
    buf = io.BytesIO()
    np.savetxt(buf, values.reshape(1, -1), delimiter=",", fmt="%.15g")
    payload = buf.getvalue()
    with gzip.open(path, "wb") as gz:
        gz.write(payload)
