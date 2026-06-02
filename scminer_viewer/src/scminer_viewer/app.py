"""Shiny-for-Python app mirroring the R Shiny / Vue layout."""

from __future__ import annotations

import re
from pathlib import Path

from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shinywidgets import output_widget, render_widget

from .data import Study, load_study
from .plots import (
    bubble_plot,
    cluster_plot,
    feature_plot,
    heatmap_plot,
    network_plot,
    violin_plot,
)
from .plots._common import ASPECT_CORRECT_JS


_SAFE_ID = re.compile(r"[^A-Za-z0-9_]")

# Constants shared between the static UI placeholders and the reactive
# panel-sync effect that swaps the placeholder for real gene panels.
_PLACEHOLDER_VALUE = "__no_gene__"
_REL_LABELS_VALUE = ["Expression", "TF", "SIG"]
_REL_LABELS_NETWORK = ["TF", "SIG"]
_REL_KEY = {"Expression": "Express_normalized",
            "TF":         "Activity_tf",
            "SIG":        "Activity_sig"}


def _safe_id(s: str) -> str:
    return _SAFE_ID.sub("_", s)


def _out_id(kind: str, gene: str, ct: str, rel: str) -> str:
    return (f"{kind}_{_safe_id(gene)}__{_safe_id(ct)}"
            f"__{_safe_id(rel)}")


def _empty_gene_panel() -> ui.NavPanel:
    """Placeholder nav_panel shown when no genes are selected."""
    return ui.nav_panel(
        ui.span({"class": "text-muted"}, "(no gene)"),
        ui.div({"class": "no-data-msg"}, "Add gene(s) to view this plot."),
        value=_PLACEHOLDER_VALUE,
    )


def _build_gene_panel(gene: str, kind: str,
                      cell_types: list[str]) -> ui.NavPanel:
    """Build the nested cell-type x relationship navset for one gene."""
    rels = _REL_LABELS_NETWORK if kind == "network" else _REL_LABELS_VALUE
    # Feature + Violin show an "All" tab that respects the global
    # cluster-visibility selection plus one tab per cell type. Network
    # omits "All" — a per-cluster TF/SIG sub-network is the meaningful
    # unit there, and an "All-clusters" superposition would just blur
    # the network structure.
    if kind == "network":
        cell_options = list(cell_types)
    else:
        cell_options = ["All"] + list(cell_types)
    ct_panels = []
    for ct in cell_options:
        rel_panels = [
            ui.nav_panel(rel, output_widget(_out_id(kind, gene, ct, rel)))
            for rel in rels
        ]
        ct_panels.append(
            ui.nav_panel(
                ct,
                ui.navset_card_pill(
                    *rel_panels,
                    id=f"{kind}_rel_{_safe_id(gene)}_{_safe_id(ct)}",
                ),
            )
        )
    return ui.nav_panel(
        gene,
        ui.navset_pill(
            *ct_panels,
            id=f"{kind}_ct_{_safe_id(gene)}",
        ),
        value=gene,
    )


def _ui_factory(study: Study, with_chrome: bool = True):
    cell_count_total = f"{study.n_cells:,}"

    # Build the study UI body once; `with_chrome=True` (run_app) wraps it
    # in a page_fluid with a <title>, while the multi-study browser embeds
    # the bare content inside its own page shell (with_chrome=False). The
    # injected script/style travel with the content either way so the
    # embedded plots still get aspect-lock.
    content = (
        ui.tags.script(ASPECT_CORRECT_JS),
        ui.tags.style("""
          .panel-card { border: 1px solid #dee2e6; border-radius: 6px;
                        background: #fff; margin-bottom: 10px; }
          .panel-card-header { background: #f5f7fa; padding: 8px 14px;
                               font-weight: 600;
                               border-bottom: 1px solid #dee2e6; }
          .panel-card-body { padding: 12px 14px; }
          .info-grid { display: grid; grid-template-columns: 1fr 1fr;
                       gap: 8px 14px; }
          .info-row { display: flex; align-items: center; gap: 8px;
                      font-size: 13px; }
          .info-label { font-weight: 600; color: #6c757d; min-width: 90px; }
          .info-value { font-family: monospace; font-size: 12px; }
          .no-data-msg { padding: 30px; text-align: center; color: #6c757d; }
          .cluster-table { width: 100%; font-size: 13px;
                            border-collapse: collapse; }
          .cluster-table thead th { font-weight: 600; color: #6c757d;
                                     padding: 6px 8px; text-align: center;
                                     border-bottom: 1px solid #dee2e6;
                                     background: #fbfcfd; }
          .cluster-table tbody td { padding: 6px 8px; vertical-align: middle;
                                     border-bottom: 1px solid #f1f3f5; }
          .cluster-table .ct-show, .cluster-table .ct-count,
          .cluster-table .ct-color { text-align: center; }
          .cluster-table .ct-count { font-family: monospace;
                                      font-size: 12px; }
          /* Hide the Plotly logo in every plot's modebar. */
          .modebar-btn--logo { display: none !important; }
        """),
        ui.div(
            {"class": "study-title"},
            ui.h4(study.meta.longTitle),
            ui.tags.small(
                {"class": "text-muted"},
                f"{study.meta.studyAbbr} - {study.meta.species} - "
                f"{cell_count_total} cells - {study.n_genes:,} genes",
            ),
        ),
        ui.br(),
        ui.layout_columns(
            ui.div(
                ui.div(
                    {"class": "panel-card"},
                    ui.div(
                        {"class": "panel-card-header"}, "Study Info & Controls"
                    ),
                    ui.div(
                        {"class": "panel-card-body"},
                        ui.div(
                            {"class": "info-grid"},
                            ui.div(
                                {"class": "info-row"},
                                ui.span({"class": "info-label"}, "Cells"),
                                ui.span(
                                    {"class": "info-value"},
                                    ui.output_text("cell_count", inline=True),
                                ),
                            ),
                            ui.div(
                                {"class": "info-row"},
                                ui.span({"class": "info-label"}, "Coordinate"),
                                ui.span(
                                    {"class": "info-value"},
                                    study.meta.coordinate,
                                ),
                            ),
                            ui.div(
                                {"class": "info-row"},
                                ui.span({"class": "info-label"}, "Dot Size"),
                                ui.input_numeric(
                                    "dot_size", None, value=4,
                                    min=0.5, max=20, step=0.5,
                                ),
                            ),
                            ui.div(
                                {"class": "info-row"},
                                ui.span({"class": "info-label"}, "Show Labels"),
                                ui.input_switch("show_labels", None, value=True),
                            ),
                            ui.div(
                                {"class": "info-row"},
                                ui.span({"class": "info-label"}, "Sampling %"),
                                ui.input_numeric(
                                    "sampling_percent", None,
                                    value=100, min=1, max=100, step=1,
                                ),
                            ),
                        ),
                    ),
                ),
                ui.div(
                    {"class": "panel-card"},
                    ui.div({"class": "panel-card-header"}, "Gene Selection"),
                    ui.div(
                        {"class": "panel-card-body"},
                        ui.input_selectize(
                            "gene_select", label=None,
                            choices=list(study.genes),
                            selected=(list(study.default_genes)
                                      if study.default_genes is not None
                                      else []),
                            multiple=True,
                            options={
                                "placeholder": "Type to add gene(s)...",
                                "maxOptions": 200,
                            },
                        ),
                    ),
                ),
            ),
            ui.div(
                {"class": "panel-card"},
                ui.div({"class": "panel-card-header"}, "Clusters"),
                ui.div(
                    {"class": "panel-card-body"},
                    ui.output_ui("clusters_table_ui"),
                ),
            ),
            col_widths=[8, 4],
        ),
        ui.navset_tab(
            ui.nav_panel("Cluster Plot", output_widget("cluster_plot")),
            ui.nav_panel("Heatmap", output_widget("heatmap_plot")),
            ui.nav_panel("Bubble Plot", output_widget("bubble_plot")),
            ui.nav_panel(
                "Feature Plot",
                ui.navset_tab(_empty_gene_panel(), id="feature_tabs"),
            ),
            ui.nav_panel(
                "Violin Plot",
                ui.navset_tab(_empty_gene_panel(), id="violin_tabs"),
            ),
            ui.nav_panel(
                "Network",
                ui.navset_tab(_empty_gene_panel(), id="network_tabs"),
            ),
            id="main_tabs",
        ),
    )
    if with_chrome:
        return ui.page_fluid(
            *content,
            title=f"scMINER Viewer - {study.meta.shortTitle}",
        )
    return ui.TagList(*content)


def _server_factory(study: Study):
    def server(input: Inputs, output: Outputs, session: Session) -> None:
        clusters_df = study.clusters.reset_index().assign(
            color_swatch=lambda d: d["color"].apply(
                lambda c: f'<span style="display:inline-block;width:24px;'
                          f'height:14px;background:{c};border:1px solid #ccc;"></span>'
            )
        )[["cellType", "count", "color_swatch"]].rename(
            columns={"cellType": "Cluster",
                     "count": "Cells",
                     "color_swatch": "Color"}
        )

        # All clusters selected initially.
        active_clusters = reactive.value(list(study.clusters.index))

        # --- Downsampling --------------------------------------------------
        # Shuffle each cluster's cells once per session; the reactive
        # then takes the first `pct%` from each shuffled vector. The
        # result is deterministic for a given (cluster, pct), so plots
        # don't flicker when the user navigates between tabs.
        import numpy as _np  # local import to keep top-level imports lean
        _rng = _np.random.default_rng(42)
        _cell_types = study.cells["cellType"].to_numpy()
        _sampling_order: dict[str, _np.ndarray] = {}
        for _ct in study.clusters.index:
            _idx = _np.flatnonzero(_cell_types == _ct)
            _rng.shuffle(_idx)
            _sampling_order[_ct] = _idx

        _MAX_CELLS = 65_000
        if study.n_cells > _MAX_CELLS:
            initial_pct = max(1, int(_MAX_CELLS / study.n_cells * 100))
        else:
            initial_pct = 100
        if initial_pct < 100:
            ui.update_numeric("sampling_percent", value=initial_pct,
                              session=session)

        @reactive.calc
        def sampling_mask() -> _np.ndarray:
            pct = int(input.sampling_percent() or 100)
            pct = max(1, min(100, pct))
            if pct >= 100:
                return _np.ones(study.n_cells, dtype=bool)
            out = _np.zeros(study.n_cells, dtype=bool)
            for _ct, idx in _sampling_order.items():
                n_keep = max(1, round(len(idx) * pct / 100))
                out[idx[:n_keep]] = True
            return out

        # Per-cluster cell-index vectors. Precomputed once so the
        # reactive count below only does a cheap mask-subset per sampling
        # change.
        _cluster_indices: dict[str, _np.ndarray] = {
            str(ct): _np.flatnonzero(_cell_types == ct)
            for ct in study.clusters.index
        }

        @reactive.calc
        def cluster_visible_counts() -> dict[str, int]:
            mask = sampling_mask()
            return {ct: int(mask[idx].sum())
                    for ct, idx in _cluster_indices.items()}

        @render.ui
        def clusters_table_ui():
            counts = cluster_visible_counts()
            sel = set(active_clusters())
            rows = []
            for ct, color, total in zip(
                study.clusters.index,
                study.clusters["color"],
                study.clusters["count"],
            ):
                ct_str = str(ct)
                safe = _safe_id(ct_str)
                # Inline checkbox + onclick that fires a Shiny event.
                # The cellType is read from a data-* attribute so we
                # don't have to escape it into the JS string literal.
                checkbox = ui.tags.input(
                    type="checkbox",
                    id=f"cluster_check_{safe}",
                    **{"data-celltype": ct_str},
                    checked=("checked" if ct_str in sel else None),
                    onclick=(
                        "Shiny.setInputValue('cluster_checkbox', "
                        "{ct: this.dataset.celltype, "
                        "checked: this.checked, t: Date.now()}, "
                        "{priority: 'event'});"
                    ),
                )
                swatch = ui.tags.span(
                    style=(f"display:inline-block;width:24px;height:14px;"
                           f"background:{color};border:1px solid #ccc;"),
                )
                cells_text = f"{counts.get(ct_str, 0):,} / {int(total):,}"
                rows.append(
                    ui.tags.tr(
                        ui.tags.td(checkbox,    class_="ct-show"),
                        ui.tags.td(ct_str,      class_="ct-name"),
                        ui.tags.td(cells_text,  class_="ct-count"),
                        ui.tags.td(swatch,      class_="ct-color"),
                    )
                )
            return ui.tags.table(
                ui.tags.thead(ui.tags.tr(
                    ui.tags.th("Show"),
                    ui.tags.th("Cluster"),
                    ui.tags.th("Cells"),
                    ui.tags.th("Color"),
                )),
                ui.tags.tbody(*rows),
                class_="cluster-table",
            )

        @reactive.effect
        @reactive.event(input.cluster_checkbox, ignore_init=True)
        def _track_cluster_checkbox():
            msg = input.cluster_checkbox()
            if not isinstance(msg, dict):
                return
            ct = msg.get("ct")
            if ct is None:
                return
            cur = list(active_clusters())
            if msg.get("checked"):
                if ct not in cur:
                    cur.append(ct)
            else:
                cur = [c for c in cur if c != ct]
            active_clusters.set(cur)

        @render.text
        def cell_count():
            ac = set(active_clusters())
            visible = (study.cells["cellType"].isin(ac).to_numpy()
                       & sampling_mask())
            return f"{int(visible.sum()):,} / {study.n_cells:,}"

        @render_widget
        def cluster_plot_widget():
            return cluster_plot(
                study,
                active_clusters=active_clusters(),
                dot_size=input.dot_size() or 4,
                show_labels=bool(input.show_labels()),
                cell_mask=sampling_mask(),
            )

        # Bind by output id used in the UI
        output(id="cluster_plot")(cluster_plot_widget)

        @render_widget
        def heatmap_widget():
            return heatmap_plot(
                study,
                genes=list(input.gene_select() or []),
                active_clusters=active_clusters(),
                cell_mask=sampling_mask(),
            )

        output(id="heatmap_plot")(heatmap_widget)

        @render_widget
        def bubble_widget():
            return bubble_plot(
                study,
                genes=list(input.gene_select() or []),
                active_clusters=active_clusters(),
                cell_mask=sampling_mask(),
            )

        output(id="bubble_plot")(bubble_widget)

        # --- 3-level nested panels: Gene → CellType → Relationship --------
        # Each of the three gene-driven tabs (Feature / Violin / Network)
        # is a *static* `ui.navset_tab` rendered once at page-build time.
        # A reactive effect mirrors `input.gene_select()` into those
        # navsets via `ui.insert_nav_panel` / `ui.remove_nav_panel`, so
        # removing a gene removes exactly that one nav_panel rather than
        # rebuilding the entire 3-level navset. This avoids the mass
        # widget teardown that triggered `OutputManager` close-model
        # errors in shinywidgets / @jupyter-widgets/html-manager.
        _cell_types_all = [str(c) for c in study.clusters.index]

        def _effective_clusters(global_sel, ct):
            if ct == "All":
                return list(global_sel)
            return [ct] if ct in global_sel else []

        def _register_outputs_for_gene(g: str) -> None:
            """Lazily register all (ct, rel) outputs for one gene.

            Feature + Violin include the "All" pseudo-cluster tab;
            Network is per-cluster only — matches the cell_options
            split in `_build_gene_panel`.
            """
            value_cts   = ["All"] + _cell_types_all
            network_cts = _cell_types_all

            for ct in value_cts:
                for rel in _REL_LABELS_VALUE:
                    for kind in ("feature", "violin"):
                        oid = _out_id(kind, g, ct, rel)
                        if oid in _registered:
                            continue
                        _registered.add(oid)
                        _register_value(
                            output, study,
                            kind=kind, gene=g, ct=ct, rel=rel,
                            rel_key=_REL_KEY[rel],
                            input=input,
                            active_clusters=active_clusters,
                            sampling_mask=sampling_mask,
                            effective_fn=_effective_clusters,
                        )
            for ct in network_cts:
                for rel in _REL_LABELS_NETWORK:
                    oid = _out_id("network", g, ct, rel)
                    if oid in _registered:
                        continue
                    _registered.add(oid)
                    _register_network(
                        output, study,
                        gene=g, ct=ct, rel=rel,
                        active_clusters=active_clusters,
                        effective_fn=_effective_clusters,
                    )

        # Plot-output registry (oid → already wired). Registrations
        # persist for the session, so re-adding a previously-removed
        # gene reuses its renderers without re-binding.
        _registered: set[str] = set()
        # Per-kind list of genes currently mounted in that navset.
        # Each kind is synced independently the first time its tab is
        # shown, so DOM insertion always happens in a visible-parent
        # context (otherwise Shiny stamps the new output as hidden and
        # never resumes its render effect).
        _mounted: dict[str, list[str]] = {
            "feature": [], "violin": [], "network": [],
        }
        _kind_for_tab = {
            "Feature Plot": "feature",
            "Violin Plot":  "violin",
            "Network":      "network",
        }

        def _sync_kind(kind: str, current: list[str]) -> None:
            navset_id = f"{kind}_tabs"
            prev = list(_mounted[kind])
            if current == prev:
                return
            prev_set = set(prev)
            cur_set = set(current)
            to_add = [g for g in current if g not in prev_set]
            to_remove = [g for g in prev if g not in cur_set]

            for g in to_add:
                _register_outputs_for_gene(g)

            if not prev and current:
                for g in current:
                    ui.insert_nav_panel(
                        navset_id,
                        _build_gene_panel(g, kind, _cell_types_all),
                        target=_PLACEHOLDER_VALUE, position="before",
                        session=session,
                    )
                ui.remove_nav_panel(
                    navset_id, target=_PLACEHOLDER_VALUE,
                    session=session,
                )
            elif prev and not current:
                ui.insert_nav_panel(
                    navset_id, _empty_gene_panel(),
                    target=prev[0], position="before",
                    session=session,
                )
                for g in prev:
                    ui.remove_nav_panel(
                        navset_id, target=g, session=session,
                    )
            else:
                for g in to_remove:
                    ui.remove_nav_panel(
                        navset_id, target=g, session=session,
                    )
                for g in to_add:
                    ui.insert_nav_panel(
                        navset_id,
                        _build_gene_panel(g, kind, _cell_types_all),
                        session=session,
                    )
            _mounted[kind] = list(current)

        @reactive.effect
        def _sync_visible_navset():
            # Re-runs on either gene change OR outer-tab change. Only
            # touches the navset whose outer tab is currently active —
            # so insert_nav_panel always lands in visible DOM.
            tab = input.main_tabs()
            kind = _kind_for_tab.get(tab)
            if kind is None:
                return
            _sync_kind(kind, list(input.gene_select() or []))

    return server


def _register_value(output, study, *, kind, gene, ct, rel, rel_key,
                    input, active_clusters, sampling_mask, effective_fn):
    """Bind a single (gene, ct, rel) feature/violin output."""
    gene_l, ct_l = gene, ct
    oid = (f"{kind}_{_safe_id(gene_l)}__{_safe_id(ct_l)}"
           f"__{_safe_id(rel)}")

    if kind == "feature":
        @render_widget
        def _w():
            eff = effective_fn(active_clusters(), ct_l)
            return feature_plot(
                study, gene=gene_l, relationship=rel_key,
                active_clusters=eff,
                dot_size=input.dot_size() or 4,
                cell_mask=sampling_mask(),
            )
    else:
        @render_widget
        def _w():
            eff = effective_fn(active_clusters(), ct_l)
            return violin_plot(
                study, gene=gene_l, relationship=rel_key,
                active_clusters=eff,
                cell_mask=sampling_mask(),
            )
    output(id=oid)(_w)


def _register_network(output, study, *, gene, ct, rel,
                       active_clusters, effective_fn):
    gene_l, ct_l, rel_l = gene, ct, rel
    oid = (f"network_{_safe_id(gene_l)}__{_safe_id(ct_l)}"
           f"__{_safe_id(rel_l)}")

    @render_widget
    def _w():
        eff = effective_fn(active_clusters(), ct_l)
        return network_plot(
            study, gene=gene_l, network_type=rel_l,
            active_clusters=eff,
        )
    output(id=oid)(_w)


def build_app(
    bundle_path: str | Path, shard_dir: str | Path | None = None
) -> App:
    """Build a Shiny App from a bundle file.

    Args:
        bundle_path: Path to a `.scminer.h5` bundle.
        shard_dir: Directory containing the per-gene shard tree
            (`expression_files/<sid>/...`, `activity_files/<sid>/...`).
            ``None`` (default) uses `Path(bundle_path).parent`.
    """
    study = load_study(bundle_path, shard_dir=shard_dir)
    return App(_ui_factory(study), _server_factory(study))


def run_app(
    bundle_path: str | Path,
    shard_dir: str | Path | None = None,
    host: str = "127.0.0.1",
    port: int = 8000,
    launch_browser: bool = True,
) -> None:
    """Build and run the Shiny app, blocking until interrupted."""
    import uvicorn

    app = build_app(bundle_path, shard_dir=shard_dir)
    if launch_browser:
        import webbrowser

        url = f"http://{host}:{port}"
        webbrowser.open(url)
    uvicorn.run(app, host=host, port=port, log_level="info",
                ws=_pick_ws_backend())


def _pick_ws_backend() -> str:
    """Choose a WebSocket backend that avoids the legacy `websockets`
    keepalive-ping AssertionError seen over LAN / HPC links."""
    try:
        from uvicorn.config import WS_PROTOCOLS
    except ImportError:
        return "auto"
    # Prefer the modern sans-IO implementation; fall back to wsproto if
    # only that is available; otherwise let uvicorn pick.
    for name in ("websockets-sansio", "wsproto"):
        if name in WS_PROTOCOLS:
            return name
    return "auto"
