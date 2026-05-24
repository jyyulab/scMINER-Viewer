"""YAML config parser + default-gene helpers.

Python port of the YAML / default-gene helpers in
scminerViewer/R/prepare_study.R: `load_study_config`,
`parse_default_genes`, `validate_default_genes`.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path
from typing import Any, Iterable, Optional


_DEFAULT_GENE_SPLIT = re.compile(r"[,;\s\n]+")


def load_study_config(config_path: str | Path) -> dict[str, Any]:
    """Parse + validate a YAML config (does not touch the RDS / TSV files).

    Required top-level keys: ``output``, ``study``, ``input``.
    Required ``study`` keys: ``ID``, ``studyAbbr``, ``longTitle``,
    ``shortTitle``. Required ``input`` key: ``expression``.

    Optional keys get sensible defaults (`cellID="cellID"`,
    `cellType="cellGroup"`, `coordinate="UMAP"`,
    `cluster_palette="npg"`, ...). Mirrors R's `load_study_config()`.
    """
    try:
        import yaml  # PyYAML
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "PyYAML is required for load_study_config(); "
            "install with `pip install pyyaml` (or `pip install "
            "scminer-viewer[prepare]`)."
        ) from exc

    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh) or {}

    for k in ("output", "study", "input"):
        if cfg.get(k) is None:
            raise ValueError(
                f"Config {config_path} is missing required top-level key: '{k}'"
            )
    for k in ("ID", "studyAbbr", "longTitle", "shortTitle"):
        if cfg["study"].get(k) is None:
            raise ValueError(
                f"Config {config_path} is missing required key: study.{k}"
            )
    if cfg["input"].get("expression") is None:
        raise ValueError(
            f"Config {config_path} is missing required key: input.expression"
        )

    # Defaults
    cfg["species"] = str(cfg.get("species") or "")
    # Coordinate label: stem of `<stem>_1` / `<stem>_2` AND the display
    # name shown in the viewer's axis labels. `coordinateName` is
    # accepted as an alias so YAMLs that already split column-names from
    # the display label keep working.
    cfg["coordinate"] = str(
        cfg.get("coordinate") or cfg.get("coordinateName") or "UMAP"
    )
    # Optional explicit per-axis column overrides. Use these when the
    # embedding columns don't follow the `<stem>_1` / `<stem>_2`
    # convention (e.g. spatial layouts with ``X`` / ``Y``). Stored as
    # `None` when absent so downstream callers can branch cleanly.
    cfg["coordinate_1"] = (
        str(cfg["coordinate_1"]) if cfg.get("coordinate_1") is not None
        else None
    )
    cfg["coordinate_2"] = (
        str(cfg["coordinate_2"]) if cfg.get("coordinate_2") is not None
        else None
    )
    cfg["cellID"] = str(cfg.get("cellID") or "cellID")
    cfg["cellType"] = str(cfg.get("cellType") or "cellGroup")
    cfg["cellGroup"] = str(cfg.get("cellGroup") or cfg["cellType"])
    cfg["geneSymbol"] = str(cfg.get("geneSymbol") or "geneSymbol")
    cfg["cluster_palette"] = str(cfg.get("cluster_palette") or "npg")
    cfg["output"] = str(cfg["output"])

    # default_genes can come from any of three keys (legacy aliases)
    raw_defaults = (
        cfg.get("default_genes")
        or (cfg.get("defaults") or {}).get("genes")
        or cfg.get("preGenes")
    )
    cfg["default_genes"] = parse_default_genes(raw_defaults)
    return cfg


def parse_default_genes(x: Any) -> Optional[list[str]]:
    """Normalize a YAML-derived default-genes value to a string list.

    Accepts:
        * `None` / empty → returns `None`
        * a list of strings → trimmed, deduplicated
        * a single string with comma / semicolon / whitespace / newline
          separators → split, trimmed, deduplicated

    First occurrence wins for duplicates.
    """
    if x is None:
        return None
    if isinstance(x, (list, tuple)):
        items = [str(v) for v in x]
    else:
        items = [str(x)]
    # Split anything that still has internal separators.
    parts: list[str] = []
    for item in items:
        for piece in _DEFAULT_GENE_SPLIT.split(item):
            piece = piece.strip()
            if piece:
                parts.append(piece)
    if not parts:
        return None
    # Dedupe while preserving order.
    seen: set[str] = set()
    out: list[str] = []
    for p in parts:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def validate_default_genes(
    default_genes: Optional[Iterable[str]],
    master_genes: Iterable[str],
    warn: bool = True,
) -> Optional[list[str]]:
    """Filter `default_genes` to those present in `master_genes`.

    Warns about drops if `warn=True`. Returns `None` if the input is
    empty or every gene was dropped.
    """
    if default_genes is None:
        return None
    defaults = [str(g) for g in default_genes]
    if not defaults:
        return None
    master = set(str(g) for g in master_genes)
    keep = [g for g in defaults if g in master]
    missing = [g for g in defaults if g not in master]
    if missing and warn:
        sample = ", ".join(missing[:6])
        warnings.warn(
            f"default_genes: {len(missing)}/{len(defaults)} not in "
            f"master gene list (dropping): {sample}",
            UserWarning,
            stacklevel=2,
        )
    return keep or None
