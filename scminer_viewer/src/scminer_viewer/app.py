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


_SAFE_ID = re.compile(r"[^A-Za-z0-9_]")


def _safe_id(s: str) -> str:
    return _SAFE_ID.sub("_", s)


def _ui_factory(study: Study):
    cell_count_total = f"{study.n_cells:,}"

    return ui.page_fluid(
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
                    ui.output_data_frame("clusters_table"),
                ),
            ),
            col_widths=[8, 4],
        ),
        ui.navset_tab(
            ui.nav_panel("Cluster Plot", output_widget("cluster_plot")),
            ui.nav_panel("Heatmap", output_widget("heatmap_plot")),
            ui.nav_panel("Bubble Plot", output_widget("bubble_plot")),
            ui.nav_panel("Feature Plot", ui.output_ui("feature_panel")),
            ui.nav_panel("Violin Plot", ui.output_ui("violin_panel")),
            ui.nav_panel("Network", ui.output_ui("network_panel")),
            id="main_tabs",
        ),
        title=f"scMINER Viewer - {study.meta.shortTitle}",
    )


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

        @render.data_frame
        def clusters_table():
            return render.DataGrid(
                clusters_df,
                selection_mode="rows",
                height="320px",
                width="100%",
            )

        @reactive.effect
        def _track_cluster_selection():
            sel = clusters_table.cell_selection()
            rows = sel.get("rows", ()) if sel else ()
            if not rows:
                active_clusters.set(list(study.clusters.index))
            else:
                active_clusters.set(
                    [study.clusters.index[i] for i in rows]
                )

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
        output("cluster_plot")(cluster_plot_widget)

        @render_widget
        def heatmap_widget():
            return heatmap_plot(
                study,
                genes=list(input.gene_select() or []),
                active_clusters=active_clusters(),
                cell_mask=sampling_mask(),
            )

        output("heatmap_plot")(heatmap_widget)

        @render_widget
        def bubble_widget():
            return bubble_plot(
                study,
                genes=list(input.gene_select() or []),
                active_clusters=active_clusters(),
                cell_mask=sampling_mask(),
            )

        output("bubble_plot")(bubble_widget)

        # --- 3-level nested panels: Gene → CellType → Relationship --------
        _cell_types_all = [str(c) for c in study.clusters.index]
        _REL_LABELS_VALUE = ["Expression", "TF", "SIG"]
        _REL_LABELS_NETWORK = ["TF", "SIG"]
        _REL_KEY = {"Expression": "Express_normalized",
                    "TF":         "Activity_tf",
                    "SIG":        "Activity_sig"}

        def _out_id(kind: str, gene: str, ct: str, rel: str) -> str:
            return (f"{kind}_{_safe_id(gene)}__{_safe_id(ct)}"
                    f"__{_safe_id(rel)}")

        def _effective_clusters(global_sel, ct):
            if ct == "All":
                return list(global_sel)
            return [ct] if ct in global_sel else []

        def _three_level_panels(genes: list[str], kind: str) -> ui.TagChild:
            if not genes:
                return ui.div(
                    {"class": "no-data-msg"},
                    "Add gene(s) to view this plot.",
                )
            rels = (_REL_LABELS_NETWORK if kind == "network"
                    else _REL_LABELS_VALUE)
            cell_options = ["All"] + _cell_types_all
            gene_panels = []
            for g in genes:
                ct_panels = []
                for ct in cell_options:
                    rel_panels = []
                    for rel in rels:
                        oid = _out_id(kind, g, ct, rel)
                        rel_panels.append(
                            ui.nav_panel(rel, output_widget(oid))
                        )
                    ct_panels.append(
                        ui.nav_panel(
                            ct,
                            ui.navset_card_pill(
                                *rel_panels,
                                id=f"{kind}_rel_{_safe_id(g)}_{_safe_id(ct)}",
                            ),
                        )
                    )
                gene_panels.append(
                    ui.nav_panel(
                        g,
                        ui.navset_pill(
                            *ct_panels,
                            id=f"{kind}_ct_{_safe_id(g)}",
                        ),
                    )
                )
            return ui.navset_tab(*gene_panels, id=f"{kind}_tabs")

        @render.ui
        def feature_panel():
            return _three_level_panels(
                list(input.gene_select() or []), "feature"
            )

        @render.ui
        def violin_panel():
            return _three_level_panels(
                list(input.gene_select() or []), "violin"
            )

        @render.ui
        def network_panel():
            return _three_level_panels(
                list(input.gene_select() or []), "network"
            )

        # Per-(gene, ct, rel) plot renderers — registered dynamically.
        # Shiny only evaluates outputs that are actually visible, so the
        # combinatorial output count doesn't translate to upfront work.
        _registered: set[str] = set()

        @reactive.effect
        def _wire_dynamic_outputs():
            genes = list(input.gene_select() or [])
            cell_options = ["All"] + _cell_types_all

            for g in genes:
                for ct in cell_options:
                    for rel in _REL_LABELS_VALUE:
                        oid = _out_id("feature", g, ct, rel)
                        if oid not in _registered:
                            _registered.add(oid)
                            _register_value(
                                output, study,
                                kind="feature", gene=g, ct=ct, rel=rel,
                                rel_key=_REL_KEY[rel],
                                input=input,
                                active_clusters=active_clusters,
                                sampling_mask=sampling_mask,
                                effective_fn=_effective_clusters,
                            )
                        oid = _out_id("violin", g, ct, rel)
                        if oid not in _registered:
                            _registered.add(oid)
                            _register_value(
                                output, study,
                                kind="violin", gene=g, ct=ct, rel=rel,
                                rel_key=_REL_KEY[rel],
                                input=input,
                                active_clusters=active_clusters,
                                sampling_mask=sampling_mask,
                                effective_fn=_effective_clusters,
                            )
                    for rel in _REL_LABELS_NETWORK:
                        oid = _out_id("network", g, ct, rel)
                        if oid not in _registered:
                            _registered.add(oid)
                            _register_network(
                                output, study,
                                gene=g, ct=ct, rel=rel,
                                active_clusters=active_clusters,
                                effective_fn=_effective_clusters,
                            )

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
    output(oid)(_w)


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
    output(oid)(_w)


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
    uvicorn.run(app, host=host, port=port, log_level="info")
