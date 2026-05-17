"""Extract cells / genes / matrices from an AnnData (or AnnData-like) object.

Python equivalent of the R `extract_cells` / `extract_genes` /
`extract_expression` / `extract_activity` helpers, which originally
consumed a Biobase `ExpressionSet`. AnnData is the de-facto Python
single-cell standard and exposes the same three slots:

    * `adata.obs`   ↔ pData (cell-level metadata)
    * `adata.var`   ↔ fData (gene-level metadata)
    * `adata.X`     ↔ exprs (matrix, cells × genes by AnnData convention)

NOTE on orientation: AnnData stores `X` as **cells × genes**; R's
ExpressionSet stores `exprs` as **genes × cells**. The shard writer
expects genes × cells, so `extract_expression` returns the transposed
matrix.

The activity helper takes a second AnnData whose `var_names` end in
`_TF` / `.TF` or `_SIG` / `.SIG`, splits on suffix, and reindexes onto
the master gene list (rows not present become zero rows).
"""

from __future__ import annotations

import re
from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp


_ACTIVITY_SUFFIX = re.compile(r"[._](TF|SIG)$")


def extract_cells(
    adata,
    cell_id_col: str = "cellID",
    cell_type_col: str = "cellGroup",
    cell_group_col: Optional[str] = None,
    coordinate_col: str = "UMAP",
) -> pd.DataFrame:
    """Build the canonical cells DataFrame from `adata`.

    Reads `adata.obs` plus the UMAP coordinates from either `adata.obs`
    (`<coordinate_col>_1` / `_2`) or `adata.obsm` (e.g. `adata.obsm["X_umap"]`).

    Args:
        adata: An `anndata.AnnData` (or duck-typed object exposing
            ``.obs`` / ``.obs_names`` / ``.obsm``).
        cell_id_col: Column in `obs` (defaults to obs_names if absent).
        cell_type_col: Column in `obs`.
        cell_group_col: Column in `obs`; defaults to `cell_type_col`.
        coordinate_col: Either the stem of the coord columns
            (`<stem>_1` / `<stem>_2`) in `obs`, OR an `obsm` key for a
            ``(n_obs, 2)`` matrix.

    Returns:
        A DataFrame with columns `cellID, cellType, cellGroup, coord1, coord2`.
    """
    if cell_group_col is None:
        cell_group_col = cell_type_col

    obs = adata.obs.copy()
    # Coerce category to plain str for downstream consistency.
    for col in obs.columns:
        if isinstance(obs[col].dtype, pd.CategoricalDtype):
            obs[col] = obs[col].astype(str)

    # Fill cellID from obs_names if the column is missing.
    if cell_id_col not in obs.columns:
        obs[cell_id_col] = obs.index.astype(str)

    if cell_type_col not in obs.columns:
        raise ValueError(
            f"adata.obs is missing required cellType column: '{cell_type_col}'"
        )

    coord1, coord2 = _resolve_coords(adata, obs, coordinate_col)

    group_vals = (
        obs[cell_group_col] if cell_group_col in obs.columns
        else obs[cell_type_col]
    )

    return pd.DataFrame({
        "cellID":    obs[cell_id_col].astype(str).to_numpy(),
        "cellType":  obs[cell_type_col].astype(str).to_numpy(),
        "cellGroup": group_vals.astype(str).to_numpy(),
        "coord1":    np.asarray(coord1, dtype=np.float64),
        "coord2":    np.asarray(coord2, dtype=np.float64),
    })


def _resolve_coords(adata, obs: pd.DataFrame, coordinate_col: str):
    c1_name = f"{coordinate_col}_1"
    c2_name = f"{coordinate_col}_2"
    if c1_name in obs.columns and c2_name in obs.columns:
        return obs[c1_name].to_numpy(), obs[c2_name].to_numpy()

    # Try obsm
    obsm = getattr(adata, "obsm", None)
    if obsm is not None:
        for key in (coordinate_col, f"X_{coordinate_col.lower()}"):
            if key in obsm:
                arr = np.asarray(obsm[key])
                if arr.ndim == 2 and arr.shape[1] >= 2:
                    return arr[:, 0], arr[:, 1]
    raise ValueError(
        f"Coordinates not found: expected obs columns '{c1_name}'/'{c2_name}' "
        f"or obsm key '{coordinate_col}' / 'X_{coordinate_col.lower()}'"
    )


def extract_genes(adata, gene_symbol_col: str = "geneSymbol") -> list[str]:
    """Read the master gene-symbol list from `adata.var`."""
    var = adata.var
    if gene_symbol_col in var.columns:
        return [str(g) for g in var[gene_symbol_col]]
    # Fall back to var_names (anndata's typical layout).
    return [str(g) for g in adata.var_names]


def extract_expression(adata, genes: Optional[list[str]] = None):
    """Return the expression matrix as **genes × cells** (matches R).

    `adata.X` is cells × genes; we transpose. Sparse matrices stay
    sparse (CSR for cheap row slicing later).
    """
    X = adata.X
    if sp.issparse(X):
        mat = X.T.tocsr()
    else:
        mat = np.asarray(X).T
    if genes is not None and mat.shape[0] != len(genes):
        raise ValueError(
            f"Expression rows ({mat.shape[0]}) != len(genes) ({len(genes)})"
        )
    return mat


def extract_activity(activity_adata, master_genes: list[str]) -> dict:
    """Split an activity AnnData into TF + SIG matrices keyed on master_genes.

    `activity_adata.var_names` (or its `geneSymbol` column) must end in
    `_TF` / `.TF` or `_SIG` / `.SIG`. Rows are split on suffix, suffix
    stripped, and reindexed onto `master_genes` (rows not present become
    zero rows; rows present but not in master are dropped).

    Returns ``{"tf": mat_or_None, "sig": mat_or_None}``, each shaped
    ``(len(master_genes), n_cells)`` if non-empty.
    """
    # Activity rows are in var_names (cells × genes AnnData layout).
    var_names = [str(v) for v in activity_adata.var_names]
    tf_mask  = np.array([bool(re.search(r"[._]TF$",  v)) for v in var_names])
    sig_mask = np.array([bool(re.search(r"[._]SIG$", v)) for v in var_names])

    if not tf_mask.any() and not sig_mask.any():
        import warnings
        warnings.warn(
            "activity_adata has no rows ending in _TF/.TF or _SIG/.SIG; "
            "returning None for both kinds.",
            UserWarning,
            stacklevel=2,
        )
        return {"tf": None, "sig": None}

    X = activity_adata.X  # cells x activity-rows
    X_T = X.T.tocsr() if sp.issparse(X) else np.asarray(X).T

    def _build(mask: np.ndarray):
        if not mask.any():
            return None
        sub = X_T[mask, :]
        sub_names = [_ACTIVITY_SUFFIX.sub("", v) for v in
                     np.array(var_names)[mask].tolist()]
        return _reindex_rows(sub, sub_names, master_genes)

    return {"tf": _build(tf_mask), "sig": _build(sig_mask)}


def _reindex_rows(mat, src_names: list[str], master_genes: list[str]):
    """Reorder `mat`'s rows to match `master_genes`; missing → zero rows.

    Sparse path constructs the output CSR directly from the source
    matrix's `data` / `indices` / `indptr` buffers, avoiding the
    LIL intermediate. ~50–100× faster on 600 K-cell × 5 K-row activity
    matrices than `lil_matrix.__setitem__` per row.
    """
    src_index = {name: i for i, name in enumerate(src_names)}
    G = len(master_genes)
    N = mat.shape[1]
    src_idx = np.fromiter(
        (src_index.get(g, -1) for g in master_genes),
        dtype=np.int64,
        count=G,
    )
    present = src_idx >= 0

    if sp.issparse(mat):
        return _reindex_csr(mat, src_idx, present, G, N)

    out = np.zeros((G, N), dtype=np.float64)
    if present.any():
        out[present, :] = np.asarray(mat)[src_idx[present], :]
    return out


def _reindex_csr(mat, src_idx: np.ndarray, present: np.ndarray, G: int, N: int):
    """Fast CSR reindex via direct buffer construction (fully vectorized).

    For each target row r where ``src_idx[r] >= 0``, copy mat row
    ``src_idx[r]``; for the rest, emit an empty row. We compute the
    output's ``indptr`` from per-row lengths (zero at absent rows), then
    use vectorized ``np.repeat`` arithmetic to map every output nnz back
    to its position in the source ``data`` / ``indices`` arrays, and do
    one fancy-index per buffer.
    """
    mat_csr = mat.tocsr()
    dtype = mat_csr.dtype if mat_csr.dtype != object else np.float64
    indptr_dtype = mat_csr.indptr.dtype
    indices_dtype = mat_csr.indices.dtype

    if not present.any():
        return sp.csr_matrix(
            (np.empty(0, dtype=dtype),
             np.empty(0, dtype=indices_dtype),
             np.zeros(G + 1, dtype=indptr_dtype)),
            shape=(G, N),
        )

    present_target = np.where(present)[0]                # ascending
    present_source = src_idx[present_target]
    src_row_lens = (
        mat_csr.indptr[present_source + 1] - mat_csr.indptr[present_source]
    ).astype(indptr_dtype, copy=False)

    # out_indptr: zero-length rows everywhere except present positions.
    row_lens = np.zeros(G, dtype=indptr_dtype)
    row_lens[present_target] = src_row_lens
    out_indptr = np.empty(G + 1, dtype=indptr_dtype)
    out_indptr[0] = 0
    np.cumsum(row_lens, out=out_indptr[1:])

    nnz = int(out_indptr[-1])
    if nnz == 0:
        return sp.csr_matrix(
            (np.empty(0, dtype=dtype),
             np.empty(0, dtype=indices_dtype),
             out_indptr),
            shape=(G, N),
        )

    # Map each output nnz position back to its source position:
    #   block_of_pos[k] = which present-row block k belongs to
    #   offset[k]       = position within that block
    #   src_pos[k]      = mat.indptr[present_source[block_of_pos[k]]] + offset[k]
    cum_present = np.empty(len(present_source) + 1, dtype=indptr_dtype)
    cum_present[0] = 0
    np.cumsum(src_row_lens, out=cum_present[1:])

    block_of_pos = np.repeat(
        np.arange(len(present_source), dtype=indptr_dtype),
        src_row_lens,
    )
    offset_in_block = (
        np.arange(nnz, dtype=indptr_dtype)
        - np.repeat(cum_present[:-1], src_row_lens)
    )
    src_pos = (
        mat_csr.indptr[present_source[block_of_pos]] + offset_in_block
    )

    out_data = mat_csr.data[src_pos]
    out_indices = mat_csr.indices[src_pos]
    return sp.csr_matrix((out_data, out_indices, out_indptr), shape=(G, N))
