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
      .cluster-table { width: 100%; font-size: 13px; border-collapse: collapse; }
      .cluster-table thead th { font-weight: 600; color: #6c757d;
                                 padding: 6px 8px; text-align: center;
                                 border-bottom: 1px solid #dee2e6;
                                 background: #fbfcfd; }
      .cluster-table tbody td { padding: 6px 8px; vertical-align: middle;
                                 border-bottom: 1px solid #f1f3f5; }
      .cluster-table .ct-show, .cluster-table .ct-count,
      .cluster-table .ct-color { text-align: center; }
      .cluster-table .ct-count { font-family: monospace; font-size: 12px; }
      .study-title { padding: 12px 4px 0 4px; }
      .back-link { padding: 8px 4px; font-size: 13px; }
      .back-link a { color: #6c757d; text-decoration: none; }
      .back-link a:hover { color: #2c3e50; }
      .browser-header { padding: 16px 4px 8px 4px; }
      .browser-header h2 { margin: 0; }
      .study-card { cursor: pointer; transition: box-shadow 0.15s ease;
                    height: 100%; }
      .study-card.is-loading { opacity: 0.55; pointer-events: none; }
      .study-card:hover { box-shadow: 0 4px 12px rgba(0,0,0,0.08); }
      .study-card a { color: inherit; text-decoration: none;
                       display: block; height: 100%; }
      .study-card-meta { color: #6c757d; font-size: 12px; }
      .gene-row { display: flex; gap: 8px; align-items: center; }
      .gene-row .form-group, .gene-row .selectize-control { flex: 1; }
      .no-data-msg { padding: 30px; text-align: center; color: #6c757d; }
      /* Full-page overlay used by run_browser():
         - shown immediately on card click (so the user gets feedback
           during the URL navigation to ?study=<sid>);
         - shown on study-page load until shiny:idle fires (covers the
           initial load_study() bundle read on the new page). */
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
          font-family: system-ui, -apple-system, 'Segoe UI', sans-serif;
      }
      #study-loading-overlay .sublabel {
          font-size: 12px; color: #6c757d;
      }
      @keyframes scm-spin { to { transform: rotate(360deg); } }
      /* Hide the Plotly logo in every plot's modebar. */
      .modebar-btn--logo { display: none !important; }
    ")
}

#' Loading overlay element used by run_browser().
#'
#' Mounted in the page shell of both the index and the study page; its
#' visibility is controlled by `.loading_overlay_js()`.
#' @noRd
.loading_overlay <- function() {
  shiny::tags$div(
    id = "study-loading-overlay", role = "status",
    `aria-live` = "polite",
    shiny::tags$div(class = "spinner"),
    shiny::tags$div(class = "label", id = "study-loading-label",
                     "Loading study…"),
    shiny::tags$div(class = "sublabel",
                     "Reading the bundle and warming caches")
  )
}

#' Inline JS that wires up the loading overlay for the browser shell.
#'
#' On the index page: clicking a card shows the overlay (the underlying
#' `?study=<sid>` navigation then completes the transition). On the
#' study page (`?study=` present in URL): overlay is shown immediately at
#' page load and cleared once Shiny finishes rendering the study UI
#' (`shiny:idle`). `shiny:error` clears the overlay so a failed load
#' doesn't lock the page.
#' @noRd
.loading_overlay_js <- function() {
  shiny::HTML("
    (function() {
      function setLabel(text) {
        var lab = document.getElementById('study-loading-label');
        if (lab) lab.textContent = text;
      }
      function show()  {
        var ov = document.getElementById('study-loading-overlay');
        if (ov) ov.classList.add('active');
      }
      function hide() {
        var ov = document.getElementById('study-loading-overlay');
        if (ov) { ov.classList.remove('active'); ov.style.display = 'none'; }
        document.querySelectorAll('.study-card.is-loading')
          .forEach(function(c) { c.classList.remove('is-loading'); });
      }
      // If we landed directly on a study URL, show the overlay during
      // the initial load_study() round-trip.
      var sid = new URLSearchParams(window.location.search).get('study');
      if (sid) {
        setLabel('Loading study…');
        show();
      }
      // Shiny dispatches `shiny:idle` / `shiny:error` via
      // `jQuery(document).trigger(...)`. Those custom events are NOT
      // delivered to native `addEventListener` listeners -- they only
      // reach `jQuery(document).on(...)` handlers -- so register through
      // jQuery (bundled with Shiny). Retry once if jQuery loads after
      // this IIFE.
      function bindShinyEvents() {
        if (typeof window.jQuery === 'undefined') {
          setTimeout(bindShinyEvents, 50);
          return;
        }
        window.jQuery(document).on('shiny:idle shiny:error', hide);
        // shiny:value fires per output completion -- clear as soon as
        // page_content is delivered, even if other outputs lag.
        window.jQuery(document).on('shiny:value', function(ev) {
          if (ev.name === 'page_content') hide();
        });
      }
      bindShinyEvents();
      // Backup signal that doesn't depend on jQuery event dispatch at
      // all: shiny.js toggles a `shiny-busy` class on <html> while a
      // flush is in flight. Hide whenever that class transitions off.
      if (window.MutationObserver) {
        new MutationObserver(function() {
          if (!document.documentElement.classList.contains('shiny-busy')) {
            hide();
          }
        }).observe(document.documentElement, {
          attributes: true, attributeFilter: ['class']
        });
      }
      // Index-page card click: dim the card and show overlay, then let
      // the <a href=\"?study=...\"> navigation continue.
      document.addEventListener('click', function(e) {
        var a = e.target.closest('.study-card a[href*=\"?study=\"]');
        if (!a) return;
        var card = a.closest('.study-card');
        if (card) {
          card.classList.add('is-loading');
          var h = card.querySelector('h5');
          if (h) setLabel(\"Loading '\" + h.textContent.trim() + \"'…\");
        }
        show();
      });
    })();
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
                  #plugins = list("remove_button"),
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
          shiny::uiOutput("clusters_table_ui")
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
    .page_chrome_css(),
    shiny::tags$style(.app_css()),
    .page_header(),
    shiny::div(class = "scv-content",
               .app_ui_content(study, with_back_link = FALSE)),
    .page_footer()
  )
}
