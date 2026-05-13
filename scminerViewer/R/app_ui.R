.app_ui <- function(study) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny is required for run_app()")
  }
  if (!requireNamespace("bslib", quietly = TRUE)) {
    stop("bslib is required for run_app()")
  }

  short_title <- study$meta$shortTitle
  long_title  <- study$meta$longTitle

  bslib::page_fluid(
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
    title = sprintf("scMINER Viewer - %s", short_title),

    shiny::tags$style(shiny::HTML("
      .panel-card { border: 1px solid #dee2e6; border-radius: 6px;
                    background: #fff; margin-bottom: 10px; }
      .panel-card-header { background: #f5f7fa; padding: 8px 14px;
                           font-weight: 600; border-bottom: 1px solid #dee2e6; }
      .panel-card-body { padding: 12px 14px; }
      .info-grid { display: grid; grid-template-columns: 1fr 1fr;
                   gap: 8px 14px; }
      .info-row { display: flex; align-items: center; gap: 8px;
                  font-size: 13px; }
      .info-label { font-weight: 600; color: #6c757d; min-width: 90px; }
      .info-value { font-family: monospace; font-size: 12px; }
      .study-title { padding: 12px 4px 0 4px; }
      .gene-row { display: flex; gap: 8px; align-items: center; }
      .gene-row .form-group, .gene-row .selectize-control { flex: 1; }
      .no-data-msg { padding: 30px; text-align: center; color: #6c757d; }
    ")),

    shiny::div(class = "study-title",
      shiny::h4(long_title),
      shiny::tags$small(class = "text-muted",
        sprintf("%s - %s - %d cells - %d genes",
                study$meta$studyAbbr,
                study$meta$species,
                nrow(study$cells),
                length(study$genes)))
    ),

    shiny::br(),

    bslib::layout_columns(
      col_widths = c(8, 4),

      # ---- Left: Study Info & Controls + Gene Selection -------------------
      shiny::div(
        shiny::div(class = "panel-card",
          shiny::div(class = "panel-card-header", "Study Info & Controls"),
          shiny::div(class = "panel-card-body",
            shiny::div(class = "info-grid",
              shiny::div(class = "info-row",
                shiny::span(class = "info-label", "Cells"),
                shiny::span(class = "info-value",
                            shiny::textOutput("cell_count", inline = TRUE))
              ),
              shiny::div(class = "info-row",
                shiny::span(class = "info-label", "Coordinate"),
                shiny::span(class = "info-value", study$meta$coordinate)
              ),
              shiny::div(class = "info-row",
                shiny::span(class = "info-label", "Dot Size"),
                shiny::numericInput("dot_size", NULL, value = 4,
                                    min = 0.5, max = 20, step = 0.5,
                                    width = "120px")
              ),
              shiny::div(class = "info-row",
                shiny::span(class = "info-label", "Show Labels"),
                shiny::checkboxInput("show_labels", NULL, value = TRUE)
              )
            )
          )
        ),
        shiny::div(class = "panel-card",
          shiny::div(class = "panel-card-header", "Gene Selection"),
          shiny::div(class = "panel-card-body",
            shiny::div(class = "gene-row",
              shiny::selectizeInput(
                "gene_select", label = NULL, choices = NULL,
                multiple = TRUE,
                options = list(
                  placeholder = "Type to add gene(s)...",
                  plugins = list("remove_button"),
                  maxOptions = 200
                )
              )
            )
          )
        )
      ),

      # ---- Right: Clusters table ------------------------------------------
      shiny::div(class = "panel-card",
        shiny::div(class = "panel-card-header", "Clusters"),
        shiny::div(class = "panel-card-body",
          DT::DTOutput("clusters_table", width = "100%")
        )
      )
    ),

    # ---- Tabset --------------------------------------------------------------
    bslib::navset_tab(
      id = "main_tabs",
      bslib::nav_panel("Cluster Plot",
        plotly::plotlyOutput("cluster_plot", height = "640px")
      ),
      bslib::nav_panel("Heatmap",
        plotly::plotlyOutput("heatmap_plot", height = "640px")
      ),
      bslib::nav_panel("Bubble Plot",
        plotly::plotlyOutput("bubble_plot", height = "640px")
      ),
      bslib::nav_panel("Feature Plot",
        shiny::uiOutput("feature_tabs_ui")
      ),
      bslib::nav_panel("Violin Plot",
        shiny::uiOutput("violin_tabs_ui")
      ),
      bslib::nav_panel("Network",
        shiny::uiOutput("network_tabs_ui")
      )
    )
  )
}
