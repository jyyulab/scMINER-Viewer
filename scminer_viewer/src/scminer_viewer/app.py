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
                            choices=list(study.genes), selected=[],
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
            n_visible = int(
                study.cells["cellType"].isin(ac).sum()
            )
            return f"{n_visible:,} / {study.n_cells:,}"

        @render_widget
        def cluster_plot_widget():
            return cluster_plot(
                study,
                active_clusters=active_clusters(),
                dot_size=input.dot_size() or 4,
                show_labels=bool(input.show_labels()),
            )

        # Bind by output id used in the UI
        output("cluster_plot")(cluster_plot_widget)

        @render_widget
        def heatmap_widget():
            return heatmap_plot(
                study,
                genes=list(input.gene_select() or []),
                active_clusters=active_clusters(),
            )

        output("heatmap_plot")(heatmap_widget)

        @render_widget
        def bubble_widget():
            return bubble_plot(
                study,
                genes=list(input.gene_select() or []),
                active_clusters=active_clusters(),
            )

        output("bubble_plot")(bubble_widget)

        # --- Dynamic per-gene nested panels --------------------------------
        def _nested_panels(genes: list[str], kind: str) -> ui.TagChild:
            if not genes:
                return ui.div(
                    {"class": "no-data-msg"},
                    "Add gene(s) to view this plot.",
                )
            panels = []
            for g in genes:
                gid = _safe_id(g)
                if kind == "feature":
                    panel_content = output_widget(f"feature_{gid}")
                elif kind == "violin":
                    panel_content = output_widget(f"violin_{gid}")
                else:  # network
                    panel_content = ui.div(
                        ui.input_radio_buttons(
                            f"netkind_{gid}", label=None,
                            choices=["TF", "SIG"], selected="TF",
                            inline=True,
                        ),
                        output_widget(f"network_{gid}"),
                    )
                panels.append(ui.nav_panel(g, panel_content))
            return ui.navset_tab(*panels, id=f"{kind}_subtabs")

        @render.ui
        def feature_panel():
            return _nested_panels(list(input.gene_select() or []), "feature")

        @render.ui
        def violin_panel():
            return _nested_panels(list(input.gene_select() or []), "violin")

        @render.ui
        def network_panel():
            return _nested_panels(list(input.gene_select() or []), "network")

        # Per-gene plot renderers — registered dynamically as the gene list
        # changes.
        _registered: set[str] = set()

        @reactive.effect
        def _wire_dynamic_outputs():
            genes = list(input.gene_select() or [])
            for g in genes:
                gid = _safe_id(g)
                if gid in _registered:
                    continue
                _registered.add(gid)

                def _make_feature(gene_local: str, gid_local: str):
                    @render_widget
                    def _fp():
                        return feature_plot(
                            study,
                            gene=gene_local,
                            active_clusters=active_clusters(),
                            dot_size=input.dot_size() or 4,
                        )
                    output(f"feature_{gid_local}")(_fp)
                _make_feature(g, gid)

                def _make_violin(gene_local: str, gid_local: str):
                    @render_widget
                    def _vp():
                        return violin_plot(
                            study,
                            gene=gene_local,
                            active_clusters=active_clusters(),
                        )
                    output(f"violin_{gid_local}")(_vp)
                _make_violin(g, gid)

                def _make_network(gene_local: str, gid_local: str):
                    @render_widget
                    def _np():
                        kind = input[f"netkind_{gid_local}"]() if (
                            f"netkind_{gid_local}" in dir(input)
                        ) else "TF"
                        return network_plot(
                            study,
                            gene=gene_local,
                            network_type=kind or "TF",
                            active_clusters=active_clusters(),
                        )
                    output(f"network_{gid_local}")(_np)
                _make_network(g, gid)

    return server


def build_app(bundle_path: str | Path) -> App:
    """Build a Shiny App from a bundle file."""
    study = load_study(bundle_path)
    return App(_ui_factory(study), _server_factory(study))


def run_app(
    bundle_path: str | Path,
    host: str = "127.0.0.1",
    port: int = 8000,
    launch_browser: bool = True,
) -> None:
    """Build and run the Shiny app, blocking until interrupted."""
    import uvicorn

    app = build_app(bundle_path)
    if launch_browser:
        import webbrowser

        url = f"http://{host}:{port}"
        webbrowser.open(url)
    uvicorn.run(app, host=host, port=port, log_level="info")
