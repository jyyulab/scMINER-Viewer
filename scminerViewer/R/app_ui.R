#' @noRd
.app_css <- function() {
  shiny::HTML("
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
      .back-link { padding: 8px 4px; font-size: 13px; }
      .back-link a { color: #6c757d; text-decoration: none; }
      .back-link a:hover { color: #2c3e50; }
      .browser-header { padding: 16px 4px 8px 4px; }
      .browser-header h2 { margin: 0; }
      .study-card { cursor: pointer; transition: box-shadow 0.15s ease;
                    height: 100%; }
      .study-card:hover { box-shadow: 0 4px 12px rgba(0,0,0,0.08); }
      .study-card a { color: inherit; text-decoration: none;
                       display: block; height: 100%; }
      .study-card-meta { color: #6c757d; font-size: 12px; }
      .gene-row { display: flex; gap: 8px; align-items: center; }
      .gene-row .form-group, .gene-row .selectize-control { flex: 1; }
      .no-data-msg { padding: 30px; text-align: center; color: #6c757d; }
    ")
}

#' Return the single-study viewer content (no page wrapper or theme).
#'
#' Used both by [run_app()] (wrapped in `page_fluid`) and by
#' [run_browser()] (embedded inside the browser's UI shell).
#' @noRd
.app_ui_content <- function(study, with_back_link = FALSE) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny is required")
  }
  if (!requireNamespace("bslib", quietly = TRUE)) {
    stop("bslib is required")
  }

  long_title <- study$meta$longTitle

  shiny::tagList(
    if (isTRUE(with_back_link)) {
      shiny::div(class = "back-link",
        shiny::a(shiny::HTML("&larr; Back to studies"), href = "?")
      )
    } else NULL,

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

      # ---- Left: Study Info & Controls + Gene Selection -----------------
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
              ),
              shiny::div(class = "info-row",
                shiny::span(class = "info-label", "Sampling %"),
                shiny::numericInput("sampling_percent", NULL,
                                    value = 100, min = 1, max = 100,
                                    step = 1, width = "120px")
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

      # ---- Right: Clusters table ---------------------------------------
      shiny::div(class = "panel-card",
        shiny::div(class = "panel-card-header", "Clusters"),
        shiny::div(class = "panel-card-body",
          DT::DTOutput("clusters_table", width = "100%")
        )
      )
    ),

    # ---- Tabset --------------------------------------------------------
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

#' @noRd
.app_ui <- function(study) {
  bslib::page_fluid(
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
    title = sprintf("scMINER Viewer - %s", study$meta$shortTitle),
    shiny::tags$style(.app_css()),
    .app_ui_content(study, with_back_link = FALSE)
  )
}
