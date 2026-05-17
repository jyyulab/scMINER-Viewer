"""Parse a scMINER networks TSV into the canonical TF / SIG DataFrames.

Python port of `read_networks` / `.read_networks_file` from
scminerViewer/R/prepare_study.R.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd


_SCORE_COLS = ("mi", "pearson", "spearman", "rho", "pvalue")


def read_networks(
    path: str | Path,
) -> dict[str, Optional[pd.DataFrame]]:
    """Read a scMINER networks TSV.

    The file must have columns ``source, target, NetworkType, CellGroup``
    plus the score columns (``mi, pearson, spearman, rho, pvalue``).
    Rows are split by ``NetworkType`` into TF and SIG groups.

    Returns:
        ``{"tf": tf_df_or_None, "sig": sig_df_or_None}`` — each DataFrame
        has columns ``source, target, cellType, mi, pearson, spearman,
        rho, pvalue``.
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t", quoting=3, dtype=str)
    if "NetworkType" not in df.columns:
        raise ValueError(f"Networks file {path} missing NetworkType column")
    if "CellGroup" not in df.columns:
        raise ValueError(f"Networks file {path} missing CellGroup column")

    src_col, tgt_col = df.columns[0], df.columns[1]
    missing_scores = [c for c in _SCORE_COLS if c not in df.columns]
    if missing_scores:
        # Fall back to positional: cols 5..9 of the original h_networks layout.
        if df.shape[1] < 9:
            raise ValueError(
                f"Networks file {path} is missing required score columns: "
                f"{', '.join(missing_scores)}"
            )
        for i, col in enumerate(_SCORE_COLS, start=4):
            if col not in df.columns:
                df[col] = df.iloc[:, i]

    def _to_canon(sub: pd.DataFrame) -> Optional[pd.DataFrame]:
        if sub is None or sub.empty:
            return None
        return pd.DataFrame({
            "source":   sub[src_col].astype(str),
            "target":   sub[tgt_col].astype(str),
            "cellType": sub["CellGroup"].astype(str),
            "mi":       pd.to_numeric(sub["mi"],       errors="coerce"),
            "pearson":  pd.to_numeric(sub["pearson"],  errors="coerce"),
            "spearman": pd.to_numeric(sub["spearman"], errors="coerce"),
            "rho":      pd.to_numeric(sub["rho"],      errors="coerce"),
            "pvalue":   pd.to_numeric(sub["pvalue"],   errors="coerce"),
        }).reset_index(drop=True)

    tf_df  = _to_canon(df[df["NetworkType"] == "TF"])
    sig_df = _to_canon(df[df["NetworkType"] == "SIG"])
    return {"tf": tf_df, "sig": sig_df}
