"""Read a study from the existing graph-import directory layout.

Python port of scminerViewer/R/graph_read.R. Reconstructs the inputs
that `prepare_study_data()` / `write_bundle()` consume from the
on-disk layout. Matrix values are *not* loaded (those are read on
demand by `Study.gene_values()`); instead, per-matrix gene **indexes**
are read from the three manifest CSVs and returned for the bundle to
embed.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


def read_graph_study(data_dir: str | Path, study_id: str) -> dict:
    """Reconstruct a study from its on-disk graph-import layout.

    Supports both layouts:
        (A) flat:    ``<data_dir>/Study/<sid>_study.tsv``
        (B) wrapped: ``<data_dir>/<sid>/Study/<sid>_study.tsv``

    (B) is what `prepare_study_data()` writes since the multi-study
    refactor.

    Returns a dict with keys: ``meta, cells, clusters, genes,
    expression_genes, activity_tf_genes, activity_sig_genes,
    network_tf, network_sig``. Any slot is ``None`` when its source
    file isn't present.
    """
    study_id = str(study_id)
    data_dir = Path(data_dir)

    study_file = data_dir / "Study" / f"{study_id}_study.tsv"
    if not study_file.exists():
        wrapped = data_dir / study_id / "Study" / f"{study_id}_study.tsv"
        if wrapped.exists():
            data_dir = data_dir / study_id
            study_file = wrapped
        else:
            raise FileNotFoundError(
                f"Study file not found: {study_file} (also tried {wrapped})"
            )

    study_row = study_file.read_text(encoding="utf-8").splitlines()[0].split("\t")
    if len(study_row) < 4:
        raise ValueError(f"Study file {study_file} has fewer than 4 fields")

    meta_csv_path = data_dir / "study_meta" / f"{study_id}_study_meta.csv"
    if not meta_csv_path.exists():
        fallback = data_dir / "study_meta" / "study_meta.csv"
        if fallback.exists():
            meta_csv_path = fallback

    meta_csv: pd.DataFrame
    if meta_csv_path.exists():
        meta_csv = pd.read_csv(meta_csv_path, dtype=str)
        if "StudyID" in meta_csv.columns:
            meta_csv = meta_csv[meta_csv["StudyID"] == study_id].reset_index(
                drop=True
            )
    else:
        meta_csv = pd.DataFrame()

    cell_file = data_dir / "Cell" / f"{study_id}_n_cell.tsv"
    if not cell_file.exists():
        raise FileNotFoundError(f"Cell file not found: {cell_file}")

    cell_df = pd.read_csv(cell_file, sep="\t", header=None, dtype=str)
    if cell_df.shape[1] >= 7:
        cells = pd.DataFrame({
            "cellID":    cell_df[1].astype(str),
            "cellType":  cell_df[2].astype(str),
            "cellGroup": cell_df[3].astype(str),
            "coord1":    pd.to_numeric(cell_df[4], errors="coerce"),
            "coord2":    pd.to_numeric(cell_df[5], errors="coerce"),
        })
        coordinate_name = str(cell_df[6].dropna().iloc[0]) if len(cell_df) else "UMAP"
    elif cell_df.shape[1] == 6:
        cells = pd.DataFrame({
            "cellID":    cell_df[0].astype(str),
            "cellType":  cell_df[1].astype(str),
            "cellGroup": cell_df[2].astype(str),
            "coord1":    pd.to_numeric(cell_df[3], errors="coerce"),
            "coord2":    pd.to_numeric(cell_df[4], errors="coerce"),
        })
        coordinate_name = str(cell_df[5].dropna().iloc[0]) if len(cell_df) else "UMAP"
    else:
        raise ValueError(f"Unexpected Cell tsv column count: {cell_df.shape[1]}")

    species_from_meta = ""
    if not meta_csv.empty and "Species" in meta_csv.columns:
        species_from_meta = str(meta_csv["Species"].iloc[0])

    meta = {
        "studyID":    study_row[0],
        "studyAbbr":  study_row[1],
        "longTitle":  study_row[2],
        "shortTitle": study_row[3],
        "species":    species_from_meta if species_from_meta else "",
        "coordinate": coordinate_name or "UMAP",
    }

    counts = (
        cells["cellType"].value_counts(sort=False)
        .rename_axis("cellType").reset_index(name="count")
    )

    if not meta_csv.empty:
        cnt_map = dict(zip(counts["cellType"], counts["count"]))
        clusters = pd.DataFrame({
            "cellType": meta_csv["CellType"],
            "count":    meta_csv["CellType"].map(cnt_map).fillna(0).astype(int),
            "color":    meta_csv.get("Color", "#888888"),
        })
        if "Label_1" in meta_csv.columns:
            clusters["label_1"] = pd.to_numeric(meta_csv["Label_1"], errors="coerce")
        if "Label_2" in meta_csv.columns:
            clusters["label_2"] = pd.to_numeric(meta_csv["Label_2"], errors="coerce")
    else:
        clusters = pd.DataFrame({
            "cellType": counts["cellType"],
            "count":    counts["count"],
            "color":    ["#888888"] * len(counts),
        })
    clusters["count"] = clusters["count"].fillna(0).astype(int)

    gene_file = data_dir / "Gene" / f"{study_id}_n_gene.tsv"
    if not gene_file.exists():
        raise FileNotFoundError(f"Gene file not found: {gene_file}")
    genes = [
        line for line in gene_file.read_text(encoding="utf-8").splitlines()
        if line
    ]

    network_tf = _read_graph_network(
        data_dir / "Network_TF_Activity" / f"{study_id}_TF.tsv"
    )
    network_sig = _read_graph_network(
        data_dir / "Network_SIG_Activity" / f"{study_id}_SIG.tsv"
    )

    # Enrich species from a manifest if still empty.
    if not meta["species"]:
        for m in (
            data_dir / "study_gene_expression" / f"{study_id}_expression.csv",
            data_dir / "study_gene_tf" / f"{study_id}_activity_tf.csv",
            data_dir / "study_gene_sig" / f"{study_id}_activity_sig.csv",
        ):
            if m.exists():
                try:
                    first = pd.read_csv(m, nrows=1)
                except Exception:
                    continue
                if "Species" in first.columns and len(first) > 0:
                    meta["species"] = str(first["Species"].iloc[0])
                    break

    expression_genes = _read_manifest_genes(
        data_dir / "study_gene_expression" / f"{study_id}_expression.csv"
    )
    activity_tf_genes = _read_manifest_genes(
        data_dir / "study_gene_tf" / f"{study_id}_activity_tf.csv"
    )
    activity_sig_genes = _read_manifest_genes(
        data_dir / "study_gene_sig" / f"{study_id}_activity_sig.csv"
    )

    return {
        "meta":               meta,
        "cells":              cells,
        "clusters":           clusters,
        "genes":              genes,
        "expression_genes":   expression_genes,
        "activity_tf_genes":  activity_tf_genes,
        "activity_sig_genes": activity_sig_genes,
        "network_tf":         network_tf,
        "network_sig":        network_sig,
    }


def _read_manifest_genes(path: Path) -> Optional[list[str]]:
    if not path.exists():
        return None
    try:
        df = pd.read_csv(path, dtype=str)
    except Exception:
        return None
    if "GeneSymbol" not in df.columns or df.empty:
        return None
    seen: set[str] = set()
    out: list[str] = []
    for g in df["GeneSymbol"]:
        s = str(g)
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


def _read_graph_network(path: Path) -> Optional[pd.DataFrame]:
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if df.shape[1] < 10:
        return None
    return pd.DataFrame({
        "source":   df[0].astype(str),
        "target":   df[1].astype(str),
        "cellType": df[4].astype(str),
        "mi":       pd.to_numeric(df[5], errors="coerce"),
        "pearson":  pd.to_numeric(df[6], errors="coerce"),
        "spearman": pd.to_numeric(df[7], errors="coerce"),
        "rho":      pd.to_numeric(df[8], errors="coerce"),
        "pvalue":   pd.to_numeric(df[9], errors="coerce"),
    })
