"""Multi-study browser: discover bundles + serve a card-grid landing page.

Mirrors the R package's `discover_studies()` / `run_browser()` /
`build_browser()`. Per-study bundles live at
`<root_dir>/<studyID>/<studyID>.scminer.h5`; the browser scans the root,
loads each bundle's meta info, and serves an index page that drills into
the standard single-study viewer when a card is clicked.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import pandas as pd
from shiny import App, Inputs, Outputs, Session, reactive, render, ui

from .app import _server_factory as _study_server, _safe_id
from .app import _ui_factory as _study_ui_factory  # noqa: F401  (kept for symmetry)
from .data import Study, load_study


def discover_studies(root_dir: str | Path) -> pd.DataFrame:
    """Return one row per `<studyID>/<studyID>.scminer.h5` under `root_dir`.

    Columns: ``studyID, studyAbbr, shortTitle, longTitle, species,
    n_cells, n_genes, n_clusters, bundle_path, study_dir``.
    Empty DataFrame when no bundles are present.
    """
    root = Path(root_dir)
    if not root.exists():
        return pd.DataFrame(columns=[
            "studyID", "studyAbbr", "shortTitle", "longTitle",
            "species", "n_cells", "n_genes", "n_clusters",
            "bundle_path", "study_dir",
        ])

    rows: list[dict] = []
    for sub in sorted(p for p in root.iterdir() if p.is_dir()):
        sid = sub.name
        bundle = sub / f"{sid}.scminer.h5"
        if not bundle.exists():
            continue
        try:
            s = load_study(bundle)
        except Exception:  # noqa: BLE001 -- skip unreadable bundles
            continue
        rows.append({
            "studyID":     str(s.meta.studyID),
            "studyAbbr":   str(s.meta.studyAbbr),
            "shortTitle":  str(s.meta.shortTitle),
            "longTitle":   str(s.meta.longTitle),
            "species":     str(s.meta.species),
            "n_cells":     s.n_cells,
            "n_genes":     s.n_genes,
            "n_clusters":  len(s.clusters),
            "bundle_path": str(bundle.resolve()),
            "study_dir":   str(sub.resolve()),
        })
    if not rows:
        return pd.DataFrame(columns=[
            "studyID", "studyAbbr", "shortTitle", "longTitle",
            "species", "n_cells", "n_genes", "n_clusters",
            "bundle_path", "study_dir",
        ])
    return pd.DataFrame(rows)


_BROWSER_CSS = """
  .browser-header { padding: 16px 4px 8px 4px; }
  .browser-header h2 { margin: 0; }
  .back-link { padding: 8px 4px; font-size: 13px; }
  .back-link a { color: #6c757d; text-decoration: none; cursor: pointer; }
  .back-link a:hover { color: #2c3e50; }
  .study-card { cursor: pointer; transition: box-shadow 0.15s ease;
                height: 100%; }
  .study-card.is-loading { opacity: 0.55; pointer-events: none; }
  .study-card:hover { box-shadow: 0 4px 12px rgba(0,0,0,0.08); }
  .study-card a { color: inherit; text-decoration: none;
                   display: block; height: 100%; }
  .study-card-meta { color: #6c757d; font-size: 12px; }
  .no-data-msg { padding: 30px; text-align: center; color: #6c757d; }

  /* Full-page overlay shown while a study bundle is being loaded.
     Triggered from the card's onclick; cleared on Shiny's `shiny:idle`. */
  #study-loading-overlay {
      display: none;
      position: fixed; inset: 0;
      background: rgba(255, 255, 255, 0.72);
      z-index: 9999;
      align-items: center; justify-content: center;
      flex-direction: column; gap: 14px;
      backdrop-filter: blur(1px);
  }
  #study-loading-overlay.active { display: flex; }
  #study-loading-overlay .spinner {
      width: 48px; height: 48px;
      border: 4px solid #dee2e6;
      border-top-color: #3C5488;
      border-radius: 50%;
      animation: scm-spin 0.9s linear infinite;
  }
  #study-loading-overlay .label {
      font-size: 14px; color: #2c3e50;
      font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
  }
  #study-loading-overlay .sublabel {
      font-size: 12px; color: #6c757d;
  }
  @keyframes scm-spin { to { transform: rotate(360deg); } }
"""


_LOADING_OVERLAY_JS = """
  (function() {
    // Hide the overlay whenever Shiny goes idle (load finished or any
    // recompute settled). Also hide on shiny:error so a failed load
    // doesn't leave the page locked.
    function clear() {
        var el = document.getElementById('study-loading-overlay');
        if (el) {
            el.classList.remove('active');
            el.dataset.sid = '';
            document.querySelectorAll('.study-card.is-loading')
                .forEach(function(c) { c.classList.remove('is-loading'); });
        }
    }

    // shiny-for-python dispatches `shiny:idle`/`shiny:error` via
    // `jQuery(document).trigger(...)`. Those custom events are NOT
    // delivered to native `addEventListener` listeners -- they only
    // reach `jQuery(document).on(...)` handlers -- so register through
    // jQuery (bundled with Shiny) once it's available.
    function bindShinyEvents() {
        if (window.jQuery) {
            window.jQuery(document).on('shiny:idle shiny:error', clear);
            return true;
        }
        return false;
    }
    if (!bindShinyEvents()) {
        // jQuery not loaded yet -- retry after DOMContentLoaded.
        document.addEventListener('DOMContentLoaded', bindShinyEvents);
    }

    // Backup signal that doesn't rely on jQuery event dispatch at all:
    // shiny.js toggles a `shiny-busy` class on <html> while a flush is
    // in flight. Hide the overlay whenever that class transitions off.
    if (window.MutationObserver) {
        new MutationObserver(function() {
            if (!document.documentElement.classList.contains('shiny-busy')) {
                clear();
            }
        }).observe(document.documentElement, {
            attributes: true, attributeFilter: ['class'],
        });
    }
  })();
"""


def _loading_overlay() -> ui.Tag:
    return ui.tags.div(
        {"id": "study-loading-overlay", "role": "status",
         "aria-live": "polite"},
        ui.tags.div({"class": "spinner"}),
        ui.tags.div({"class": "label", "id": "study-loading-label"},
                     "Loading study…"),
        ui.tags.div({"class": "sublabel"},
                     "Reading the bundle and warming caches"),
    )


def _study_card(row) -> ui.Tag:
    safe = _safe_id(row["studyID"])
    short_title = str(row["shortTitle"]).replace("\\", "\\\\").replace("'", "\\'")
    return ui.div(
        {"class": "card study-card", "data-sid": row["studyID"]},
        ui.tags.a(
            {
                "id": f"open_study_{safe}",
                "onclick": (
                    # Show the overlay synchronously, then notify the
                    # server. Shiny:idle (or :error) clears both.
                    "var ov = document.getElementById('study-loading-overlay');"
                    f"var lab = document.getElementById('study-loading-label');"
                    f"if (lab) lab.textContent = 'Loading ‘{short_title}’…';"
                    "if (ov) { ov.dataset.sid = "
                    f"{row['studyID']!r}; ov.classList.add('active'); }}"
                    "this.closest('.study-card').classList.add('is-loading');"
                    f"Shiny.setInputValue('open_study', "
                    f"{{ sid: {row['studyID']!r}, t: Date.now() }}, "
                    f"{{priority: 'event'}});"
                ),
            },
            ui.div(
                {"class": "card-body"},
                ui.tags.h5({"class": "mb-1"}, row["shortTitle"]),
                ui.div(
                    {"class": "study-card-meta"},
                    f"{row['studyAbbr']} - {row['species']} - "
                    f"{row['n_cells']:,} cells x {row['n_genes']:,} genes x "
                    f"{row['n_clusters']} clusters",
                ),
                ui.tags.p(
                    {"class": "mt-2 mb-0 small text-muted",
                     "title": row["longTitle"]},
                    row["longTitle"],
                ),
            ),
        ),
    )


def _index_content(studies: pd.DataFrame) -> ui.TagChild:
    if studies.empty:
        return ui.tags.div(
            ui.div({"class": "browser-header"}, ui.tags.h2("Studies")),
            ui.div(
                {"class": "no-data-msg"},
                "No studies found. Run scminerViewer::prepare_study() or "
                "the scminer-viewer migration helper to create one.",
            ),
        )
    cards = [_study_card(row) for _, row in studies.iterrows()]
    grid = ui.layout_column_wrap(*cards, width="320px", gap="12px")
    return ui.tags.div(
        ui.div(
            {"class": "browser-header"},
            ui.tags.h2("Studies"),
            ui.tags.small(
                {"class": "text-muted"},
                f"{len(studies)} "
                f"{'study' if len(studies) == 1 else 'studies'} available",
            ),
        ),
        grid,
    )


def build_browser(
    root_dir: str | Path, shard_dir: str | Path | None = None
) -> App:
    """Build a multi-study Shiny App without launching it.

    Args:
        root_dir: Directory containing one `<studyID>/<studyID>.scminer.h5`
            per study.
        shard_dir: Where the per-gene shard tree lives. ``None`` (default)
            uses each bundle's parent directory.
    """
    root_dir = Path(root_dir).resolve()
    studies = discover_studies(root_dir)
    shard_dir_p = Path(shard_dir).resolve() if shard_dir is not None else None

    # Build a quick studyID → bundle_path map for routing
    bundle_by_id = {row["studyID"]: row["bundle_path"]
                    for _, row in studies.iterrows()}

    app_ui = ui.page_fluid(
        ui.tags.style(_BROWSER_CSS),
        ui.tags.style("""
          .panel-card { border: 1px solid #dee2e6; border-radius: 6px;
                        background: #fff; margin-bottom: 10px; }
          .panel-card-header { background: #f5f7fa; padding: 8px 14px;
                               font-weight: 600;
                               border-bottom: 1px solid #dee2e6; }
          .panel-card-body { padding: 12px 14px; }
        """),
        _loading_overlay(),
        ui.tags.script(_LOADING_OVERLAY_JS),
        ui.output_ui("page_content"),
        title="scMINER Viewer",
    )

    def server(input: Inputs, output: Outputs, session: Session) -> None:
        current_sid = reactive.value(None)

        @reactive.effect
        @reactive.event(input.open_study)
        def _route_to_study():
            msg = input.open_study()
            if not msg:
                return
            sid = msg.get("sid") if isinstance(msg, dict) else None
            if sid and sid in bundle_by_id:
                current_sid.set(sid)

        @reactive.effect
        @reactive.event(input.back_to_index, ignore_init=True)
        def _route_back():
            current_sid.set(None)

        # Lazy: a study is loaded only after the user clicks its card.
        loaded_studies: dict[str, Study] = {}
        wired_studies: set[str] = set()

        def load(sid: str) -> Optional[Study]:
            if sid not in loaded_studies:
                try:
                    loaded_studies[sid] = load_study(
                        bundle_by_id[sid], shard_dir=shard_dir_p
                    )
                except Exception:  # noqa: BLE001
                    return None
            return loaded_studies[sid]

        @render.ui
        def page_content():
            sid = current_sid()
            if sid is None:
                return _index_content(studies)
            study = load(sid)
            if study is None:
                return ui.div(
                    {"class": "no-data-msg"},
                    f"Could not load study {sid}",
                )
            # Standard single-study UI body — the same content the
            # `run_app()` viewer renders, with a "back" link prepended.
            return ui.tags.div(
                ui.div(
                    {"class": "back-link"},
                    ui.tags.a(
                        {
                            "id": "back_to_index_link",
                            "onclick": (
                                "Shiny.setInputValue('back_to_index', "
                                "Date.now(), {priority: 'event'});"
                            ),
                        },
                        "← Back to studies",
                    ),
                ),
                _study_inner_ui(study),
            )

        # Register the single-study server logic for each study the first
        # time it's opened. Shiny's reactive scoping means stale outputs
        # for a previously-opened study still update if their inputs
        # change, but they're never visible (the renderUI swaps them
        # out), so the cost is just background reactivity.
        @reactive.effect
        def _wire_active_study():
            sid = current_sid()
            if sid is None or sid in wired_studies:
                return
            study = load(sid)
            if study is None:
                return
            _study_server(study)(input, output, session)
            wired_studies.add(sid)

    return App(app_ui, server)


def _study_inner_ui(study: Study):
    """The single-study viewer UI body without page_fluid/title wrapping.

    Mirrors scminer_viewer.app._ui_factory but skips ui.page_fluid so
    we can embed inside the browser shell.
    """
    from .app import _ui_factory as _full_factory
    # Reuse the same UI factory — the page_fluid wrapper is harmless
    # when nested (Bootstrap container-fluid inside container-fluid is
    # fine). This keeps the two surfaces in sync without duplication.
    inner = _full_factory(study)
    return inner


def run_browser(
    root_dir: str | Path,
    shard_dir: str | Path | None = None,
    host: str = "127.0.0.1",
    port: int = 8000,
    launch_browser: bool = True,
) -> None:
    """Build and run the multi-study browser, blocking until interrupted."""
    import uvicorn

    app = build_browser(root_dir, shard_dir=shard_dir)
    if launch_browser:
        import webbrowser
        webbrowser.open(f"http://{host}:{port}")
    from .app import _pick_ws_backend
    uvicorn.run(app, host=host, port=port, log_level="info",
                ws=_pick_ws_backend())
