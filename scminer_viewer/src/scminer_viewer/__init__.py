"""scMINER Viewer — Python webui + data prep for scMINER study bundles."""

from .data import Study, load_study
from .app import build_app, run_app
from .browser import discover_studies, build_browser, run_browser

__all__ = [
    "Study", "load_study",
    "build_app", "run_app",
    "discover_studies", "build_browser", "run_browser",
    # Data prep — lazy-imported via attribute access so the heavy deps
    # (anndata, pyyaml, scipy.sparse write paths) only load when used.
    "prepare",
    "prepare_study", "prepare_study_from_anndata", "prepare_study_data",
    "write_bundle", "read_graph_study", "fill_clusters",
    "load_study_config", "parse_default_genes", "validate_default_genes",
    "read_networks",
]
__version__ = "0.3.0"


def __getattr__(name):  # PEP 562 — lazy prepare imports
    _prepare_exports = {
        "prepare",
        "prepare_study", "prepare_study_from_anndata", "prepare_study_data",
        "write_bundle", "read_graph_study", "fill_clusters",
        "load_study_config", "parse_default_genes", "validate_default_genes",
        "read_networks",
    }
    if name in _prepare_exports:
        from . import prepare as _prepare
        if name == "prepare":
            return _prepare
        return getattr(_prepare, name)
    raise AttributeError(f"module 'scminer_viewer' has no attribute {name!r}")
