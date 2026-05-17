"""Prepare scMINER study bundles + shard trees from AnnData inputs.

This module is the Python equivalent of the R-side `prepare_study`
family in `scminerViewer/`. Public surface:

* :func:`prepare_study`              — YAML-driven (config + h5ad files)
* :func:`prepare_study_from_anndata` — AnnData input(s)
* :func:`prepare_study_data`         — already-extracted structures
* :func:`write_bundle`               — HDF5 writer alone
* :func:`read_graph_study`           — round-trip the on-disk layout
* :func:`fill_clusters`              — palette colors + centroid labels
* :func:`load_study_config`          — parse the YAML schema
* :func:`parse_default_genes`        — normalise the YAML value
* :func:`validate_default_genes`     — filter against a master gene list
* :func:`read_networks`              — TSV → ``{tf, sig}`` DataFrames
* :func:`extract_cells` / :func:`extract_genes` /
  :func:`extract_expression` / :func:`extract_activity` — AnnData
  extractors for callers that want granular control.

Optional deps: ``anndata`` (only for `prepare_study_from_anndata` and
`prepare_study`), ``pyyaml`` (only for `load_study_config`). Install
together via ``pip install scminer-viewer[prepare]``.
"""

from __future__ import annotations

from ._bundle import write_bundle
from ._clusters import fill_clusters
from ._config import (
    load_study_config,
    parse_default_genes,
    validate_default_genes,
)
from ._eset import (
    extract_activity,
    extract_cells,
    extract_expression,
    extract_genes,
)
from ._graph_read import read_graph_study
from ._networks import read_networks
from .orchestrator import (
    prepare_study,
    prepare_study_data,
    prepare_study_from_anndata,
)


__all__ = [
    "prepare_study",
    "prepare_study_from_anndata",
    "prepare_study_data",
    "write_bundle",
    "read_graph_study",
    "fill_clusters",
    "load_study_config",
    "parse_default_genes",
    "validate_default_genes",
    "read_networks",
    "extract_cells",
    "extract_genes",
    "extract_expression",
    "extract_activity",
]
