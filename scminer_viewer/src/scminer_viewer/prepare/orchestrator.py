"""Top-level prepare-study orchestrators.

Three entry points, mirroring the R package's `prepare_study` family:

* `prepare_study_data` — lowest-level; takes already-extracted
  structures (meta dict, cells DataFrame, matrices). No YAML / AnnData
  dependency.
* `prepare_study_from_anndata` — accepts one or two `anndata.AnnData`
  objects (expression + optional activity) and pulls the canonical
  cells / genes / matrices out internally. Replaces the R
  `prepare_study_from_eset` (Biobase ExpressionSet) entry point.
* `prepare_study` — drop-in YAML driver: reads a config file, loads
  the referenced `.h5ad` / `.tsv` files, then calls
  `prepare_study_from_anndata`.

All three accept `emit=("graph", "bundle")`; either or both. The
bundle always lands at `<out_dir>/<studyID>/<studyID>.scminer.h5`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence, Union

import numpy as np
import pandas as pd

from ._bundle import write_bundle
from ._clusters import fill_clusters
from ._config import load_study_config, validate_default_genes
from ._eset import (
    extract_activity,
    extract_cells,
    extract_expression,
    extract_genes,
)
from ._graph import (
    ensure_graph_tree,
    write_graph_cells,
    write_graph_clusters,
    write_graph_genes,
    write_graph_networks,
    write_graph_shards,
    write_graph_study,
)
from ._networks import read_networks


_EMIT_VALUES = ("graph", "bundle")


def prepare_study_data(
    out_dir: Union[str, Path],
    meta: Mapping[str, Any],
    cells: pd.DataFrame,
    genes: Sequence[str],
    *,
    clusters: Optional[pd.DataFrame] = None,
    expression=None,
    activity_tf=None,
    activity_sig=None,
    network_tf: Optional[pd.DataFrame] = None,
    network_sig: Optional[pd.DataFrame] = None,
    default_genes: Optional[Sequence[str]] = None,
    cluster_palette: str = "npg",
    emit: Sequence[str] = _EMIT_VALUES,
    verbose: bool = False,
) -> dict:
    """Lowest-level orchestrator — write graph layout and/or bundle.

    Args:
        out_dir: Root directory; each study lands under ``<out_dir>/<studyID>/``.
        meta: Mapping with ``studyID, studyAbbr, longTitle, shortTitle,
            species, coordinate``.
        cells: DataFrame with ``cellID, cellType, cellGroup?, coord1, coord2``.
        genes: Master gene-symbol sequence (in row order of the matrices).
        clusters: Optional clusters DataFrame; missing columns are
            auto-filled from `cells` via `fill_clusters`.
        expression / activity_tf / activity_sig: Optional 2-D matrices
            (genes × cells) — ndarray, scipy.sparse, or DataFrame.
        network_tf / network_sig: Optional DataFrames with the canonical
            columns.
        default_genes: Genes auto-loaded on app launch.
        cluster_palette: Name of the palette used to fill any missing
            cluster colors.
        emit: Subset of ``("graph", "bundle")``.
        verbose: Emit progress while writing shards.

    Returns:
        ``{"out_dir": ..., "root_dir": ..., "bundle_path": ...}``.
    """
    emit = tuple(emit)
    bad = [e for e in emit if e not in _EMIT_VALUES]
    if bad:
        raise ValueError(
            f"emit must be a subset of {_EMIT_VALUES}; got: {bad}"
        )

    if not isinstance(meta, Mapping):
        raise TypeError("`meta` must be a Mapping")
    if not isinstance(cells, pd.DataFrame):
        raise TypeError("`cells` must be a pandas DataFrame")
    genes = list(genes)

    clusters = fill_clusters(cells, clusters, palette=cluster_palette)

    # Validate matrix row counts against genes.
    for name, mat in (
        ("expression", expression),
        ("activity_tf", activity_tf),
        ("activity_sig", activity_sig),
    ):
        if mat is None:
            continue
        nrows = mat.shape[0] if hasattr(mat, "shape") else len(mat)
        if nrows != len(genes):
            raise ValueError(
                f"{name} row count ({nrows}) does not match "
                f"len(genes) ({len(genes)})"
            )

    study_id = str(meta["studyID"])
    out_dir = Path(out_dir)
    study_out = out_dir / study_id
    study_out.mkdir(parents=True, exist_ok=True)

    if "graph" in emit:
        if verbose:
            print(f"Writing graph-import layout to {study_out}")
        ensure_graph_tree(study_out)
        write_graph_study(study_out, meta)
        write_graph_clusters(study_out, meta, clusters)
        write_graph_genes(study_out, meta, genes)
        write_graph_cells(study_out, meta, cells)
        write_graph_networks(study_out, meta, network_tf, network_sig)

        exp_root = f"expression_files/{study_id}"
        act_root = f"activity_files/{study_id}"

        cell_ids = cells["cellID"].astype(str).tolist()

        if expression is not None:
            write_graph_shards(
                study_out, meta, expression,
                kind=exp_root,
                meta_kind=exp_root,
                manifest_dir="study_gene_expression",
                manifest_name="expression",
                cell_ids=cell_ids,
                type_label="Expression",
                genes=genes,
                verbose=verbose,
            )
        if activity_tf is not None:
            write_graph_shards(
                study_out, meta, activity_tf,
                kind=f"{act_root}/TF",
                meta_kind=act_root,
                manifest_dir="study_gene_tf",
                manifest_name="activity_tf",
                cell_ids=cell_ids,
                type_label="TF",
                genes=genes,
                verbose=verbose,
            )
        if activity_sig is not None:
            write_graph_shards(
                study_out, meta, activity_sig,
                kind=f"{act_root}/SIG",
                meta_kind=act_root,
                manifest_dir="study_gene_sig",
                manifest_name="activity_sig",
                cell_ids=cell_ids,
                type_label="SIG",
                genes=genes,
                verbose=verbose,
            )

    bundle_path = None
    if "bundle" in emit:
        bundle_path = study_out / f"{study_id}.scminer.h5"
        if verbose:
            print(f"Writing bundle to {bundle_path}")
        write_bundle(
            bundle_path=bundle_path,
            meta=meta,
            cells=cells,
            clusters=clusters,
            genes=genes,
            expression_genes=genes if expression is not None else None,
            activity_tf_genes=genes if activity_tf is not None else None,
            activity_sig_genes=genes if activity_sig is not None else None,
            default_genes=default_genes,
            network_tf=network_tf,
            network_sig=network_sig,
            overwrite=True,
        )

    return {
        "out_dir": str(study_out),
        "root_dir": str(out_dir),
        "bundle_path": str(bundle_path) if bundle_path else None,
    }


def prepare_study_from_anndata(
    out_dir: Union[str, Path],
    expression_adata,
    *,
    activity_adata=None,
    networks_path: Optional[Union[str, Path]] = None,
    meta: Mapping[str, Any],
    cell_id_col: str = "cellID",
    cell_type_col: str = "cellGroup",
    cell_group_col: Optional[str] = None,
    coordinate_col: str = "UMAP",
    gene_symbol_col: str = "geneSymbol",
    clusters: Optional[pd.DataFrame] = None,
    cluster_palette: str = "npg",
    default_genes: Optional[Sequence[str]] = None,
    emit: Sequence[str] = _EMIT_VALUES,
    verbose: bool = False,
) -> dict:
    """Mid-level orchestrator: accepts AnnData objects directly.

    Internally calls :func:`extract_cells`, :func:`extract_genes`,
    :func:`extract_expression`, optionally :func:`extract_activity`,
    optionally :func:`read_networks`, then :func:`prepare_study_data`.
    """
    cells = extract_cells(
        expression_adata,
        cell_id_col=cell_id_col,
        cell_type_col=cell_type_col,
        cell_group_col=cell_group_col,
        coordinate_col=coordinate_col,
    )
    genes = extract_genes(expression_adata, gene_symbol_col=gene_symbol_col)
    expression = extract_expression(expression_adata, genes=genes)

    if activity_adata is not None:
        act = extract_activity(activity_adata, master_genes=genes)
        activity_tf = act["tf"]
        activity_sig = act["sig"]
    else:
        activity_tf = activity_sig = None

    if networks_path is not None:
        nets = read_networks(networks_path)
        network_tf = nets["tf"]
        network_sig = nets["sig"]
    else:
        network_tf = network_sig = None

    meta = dict(meta)
    meta.setdefault("coordinate", coordinate_col)
    default_genes = validate_default_genes(default_genes, genes)

    return prepare_study_data(
        out_dir=out_dir,
        meta=meta,
        cells=cells,
        clusters=clusters,
        genes=genes,
        expression=expression,
        activity_tf=activity_tf,
        activity_sig=activity_sig,
        network_tf=network_tf,
        network_sig=network_sig,
        default_genes=default_genes,
        cluster_palette=cluster_palette,
        emit=emit,
        verbose=verbose,
    )


def prepare_study(
    config_path: Union[str, Path],
    *,
    emit: Sequence[str] = _EMIT_VALUES,
    verbose: bool = False,
) -> dict:
    """High-level YAML-driven entry point.

    Reads `config_path` via :func:`load_study_config`, opens the
    referenced AnnData files (`input.expression`, optional
    `input.activity`), then calls :func:`prepare_study_from_anndata`.

    Input file paths are resolved relative to the YAML's parent
    directory if not absolute.
    """
    try:
        import anndata as ad
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "anndata is required for prepare_study(); install with "
            "`pip install anndata` (or `pip install scminer-viewer[prepare]`)."
        ) from exc

    cfg = load_study_config(config_path)
    base = Path(config_path).resolve().parent

    def _resolve(p: str) -> Path:
        path = Path(p)
        return path if path.is_absolute() else (base / path).resolve()

    expression_path = _resolve(cfg["input"]["expression"])
    expression_adata = ad.read_h5ad(expression_path)

    activity_adata = None
    activity_in = cfg["input"].get("activity")
    if activity_in is not None:
        activity_adata = ad.read_h5ad(_resolve(activity_in))

    networks_path = cfg["input"].get("networks")
    if networks_path is not None:
        networks_path = _resolve(networks_path)

    meta = {
        "studyID":    str(cfg["study"]["ID"]),
        "studyAbbr":  str(cfg["study"]["studyAbbr"]),
        "longTitle":  str(cfg["study"]["longTitle"]),
        "shortTitle": str(cfg["study"]["shortTitle"]),
        "species":    cfg["species"],
        "coordinate": cfg["coordinate"],
    }

    return prepare_study_from_anndata(
        out_dir=cfg["output"],
        expression_adata=expression_adata,
        activity_adata=activity_adata,
        networks_path=networks_path,
        meta=meta,
        cell_id_col=cfg["cellID"],
        cell_type_col=cfg["cellType"],
        cell_group_col=cfg["cellGroup"],
        coordinate_col=cfg["coordinate"],
        gene_symbol_col=cfg["geneSymbol"],
        cluster_palette=cfg["cluster_palette"],
        default_genes=cfg["default_genes"],
        emit=emit,
        verbose=verbose,
    )
