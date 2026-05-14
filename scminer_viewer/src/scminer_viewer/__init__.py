"""scMINER Viewer — Python webui for scMINER study bundles."""

from .data import Study, load_study
from .app import build_app, run_app
from .browser import discover_studies, build_browser, run_browser

__all__ = [
    "Study", "load_study",
    "build_app", "run_app",
    "discover_studies", "build_browser", "run_browser",
]
__version__ = "0.2.0"
