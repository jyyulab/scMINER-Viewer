#' Discover scMINER study bundles under a root directory.
#'
#' Scans `root_dir` for subdirectories named after a studyID and looks
#' for `<studyID>/<studyID>.scminer.h5` inside each (the layout produced
#' by [prepare_study_data()]). Returns one row per discovered bundle.
#'
#' @param root_dir Directory containing one subfolder per study.
#'
#' @return A data.frame with columns `studyID`, `studyAbbr`, `shortTitle`,
#'   `longTitle`, `species`, `n_cells`, `n_genes`, `n_clusters`,
#'   `bundle_path`, `study_dir`. Empty data.frame if no bundles are found.
#' @export
discover_studies <- function(root_dir) {
  if (!dir.exists(root_dir)) {
    return(data.frame(
      studyID = character(0), studyAbbr = character(0),
      shortTitle = character(0), longTitle = character(0),
      species = character(0),
      n_cells = integer(0), n_genes = integer(0), n_clusters = integer(0),
      bundle_path = character(0), study_dir = character(0),
      stringsAsFactors = FALSE
    ))
  }
  subdirs <- list.dirs(root_dir, recursive = FALSE)
  entries <- list()
  for (d in subdirs) {
    sid    <- basename(d)
    bundle <- file.path(d, paste0(sid, ".scminer.h5"))
    if (!file.exists(bundle)) next
    s <- tryCatch(load_study(bundle), error = function(e) NULL)
    if (is.null(s)) next
    entries[[length(entries) + 1L]] <- data.frame(
      studyID     = as.character(s$meta$studyID),
      studyAbbr   = as.character(s$meta$studyAbbr),
      shortTitle  = as.character(s$meta$shortTitle),
      longTitle   = as.character(s$meta$longTitle),
      species     = as.character(s$meta$species),
      n_cells     = nrow(s$cells),
      n_genes     = length(s$genes),
      n_clusters  = nrow(s$clusters),
      bundle_path = normalizePath(bundle, mustWork = FALSE),
      study_dir   = normalizePath(d,      mustWork = FALSE),
      stringsAsFactors = FALSE
    )
  }
  if (length(entries) == 0L) {
    return(data.frame(
      studyID = character(0), studyAbbr = character(0),
      shortTitle = character(0), longTitle = character(0),
      species = character(0),
      n_cells = integer(0), n_genes = integer(0), n_clusters = integer(0),
      bundle_path = character(0), study_dir = character(0),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, entries)
  out[order(out$studyID), , drop = FALSE]
}

#' Render the card-grid index of studies (no page wrapper).
#' @noRd
.browser_index_content <- function(studies) {
  if (!requireNamespace("bslib", quietly = TRUE)) {
    stop("bslib is required for run_browser()")
  }
  if (nrow(studies) == 0) {
    return(shiny::tagList(
      shiny::div(class = "browser-header",
                 shiny::h2("Studies")),
      shiny::div(class = "no-data-msg",
                 "No studies found. Run prepare_study() to create one.")
    ))
  }
  cards <- lapply(seq_len(nrow(studies)), function(i) {
    e <- studies[i, , drop = FALSE]
    href <- paste0("?study=", utils::URLencode(e$studyID, reserved = TRUE))
    bslib::card(
      class = "study-card",
      bslib::card_body(
        shiny::a(href = href,
          shiny::h5(e$shortTitle, class = "mb-1"),
          shiny::div(class = "study-card-meta",
            paste0(e$studyAbbr, " - ", e$species,
                   " - ", format(e$n_cells, big.mark = ","),
                   " cells x ", format(e$n_genes, big.mark = ","),
                   " genes x ", e$n_clusters, " clusters")
          ),
          shiny::p(class = "mt-2 mb-0 small text-muted",
                   shiny::tagAppendAttributes(
                     shiny::span(e$longTitle),
                     title = e$longTitle
                   ))
        )
      )
    )
  })
  grid <- do.call(bslib::layout_column_wrap,
                  c(list(width = "320px", gap = "12px"), cards))
  shiny::tagList(
    shiny::div(class = "browser-header",
      shiny::h2("Studies"),
      shiny::tags$small(class = "text-muted",
        sprintf("%d %s available", nrow(studies),
                if (nrow(studies) == 1) "study" else "studies"))
    ),
    grid
  )
}

#' Launch the multi-study scMINER Viewer browser.
#'
#' Serves a card-grid landing page of all studies discovered under
#' `root_dir` (subfolder-per-study layout produced by
#' [prepare_study_data()]). Clicking a card navigates to
#' `?study=<studyID>`, which loads that study's bundle and renders the
#' standard single-study viewer with a "← Back to studies" link.
#'
#' @param root_dir Directory containing one `<studyID>/<studyID>.scminer.h5`
#'   per study.
#' @param shard_dir Directory containing the per-gene shard trees
#'   (`expression_files/<studyID>/...`, `activity_files/<studyID>/...`).
#'   If `NULL` (default), each bundle defaults to `dirname(bundle_path)`
#'   — appropriate when bundles and shards are co-located inside each
#'   study folder. Pass an explicit path when the shards live elsewhere
#'   (e.g. an existing `data/example/` tree alongside fresh per-study
#'   bundle directories).
#' @param host Host interface to bind to.
#' @param port Port to bind to (NULL = pick a free one).
#' @param launch_browser Logical; open the system browser.
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return Called for side effect. Returns `invisible(NULL)`.
#' @export
run_browser <- function(root_dir,
                        shard_dir = NULL,
                        host = "127.0.0.1",
                        port = NULL,
                        launch_browser = interactive(),
                        ...) {
  for (pkg in c("shiny", "bslib", "plotly", "DT", "visNetwork")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package %s is required for run_browser()", pkg))
    }
  }
  root_dir <- normalizePath(root_dir, mustWork = TRUE)
  studies  <- discover_studies(root_dir)
  if (nrow(studies) == 0L) {
    warning("No studies found at: ", root_dir,
            " (expected <studyID>/<studyID>.scminer.h5 subdirs)",
            call. = FALSE)
  }

  ui <- bslib::page_fluid(
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
    title = "scMINER Viewer",
    .page_chrome_css(),
    shiny::tags$style(.app_css()),
    .page_header(),
    .loading_overlay(),
    shiny::tags$script(.loading_overlay_js()),
    shiny::div(class = "scv-content",
               shiny::uiOutput("page_content")),
    .page_footer()
  )

  server <- function(input, output, session) {
    # Read the URL once at session start. Each ?study=<id> click navigates
    # to a fresh URL, which Shiny treats as a new session — so we don't
    # need to handle in-session study switching.
    chosen <- shiny::reactive({
      qs <- shiny::parseQueryString(session$clientData$url_search %||% "")
      sid <- qs$study
      if (is.null(sid) || !nzchar(sid)) return(NULL)
      hit <- studies[studies$studyID == sid, , drop = FALSE]
      if (nrow(hit) == 0L) return(NULL)
      tryCatch(load_study(hit$bundle_path[1], shard_dir = shard_dir),
               error = function(e) NULL)
    })

    output$page_content <- shiny::renderUI({
      s <- chosen()
      if (is.null(s)) {
        .browser_index_content(studies)
      } else {
        .app_ui_content(s, with_back_link = TRUE)
      }
    })

    # Wire up the study-specific server logic. We use observeEvent gated
    # on chosen() so the logic runs exactly once for the chosen study.
    wired <- shiny::reactiveVal(FALSE)
    shiny::observe({
      if (isTRUE(wired())) return()
      s <- chosen()
      if (is.null(s)) return()
      .app_server_logic(s, input, output, session)
      wired(TRUE)
    })
  }

  app <- shiny::shinyApp(ui = ui, server = server)
  shiny::runApp(app, host = host, port = port,
                launch.browser = launch_browser, ...)
}

#' Build a multi-study browser app without launching it.
#'
#' Useful for tests and embedding. See [run_browser()] for details.
#'
#' @inheritParams run_browser
#' @return A `shiny.appobj`.
#' @export
build_browser <- function(root_dir, shard_dir = NULL) {
  for (pkg in c("shiny", "bslib", "plotly", "DT", "visNetwork")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package %s is required for build_browser()", pkg))
    }
  }
  root_dir <- normalizePath(root_dir, mustWork = TRUE)
  studies  <- discover_studies(root_dir)
  ui <- bslib::page_fluid(
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
    title = "scMINER Viewer",
    .page_chrome_css(),
    shiny::tags$style(.app_css()),
    .page_header(),
    .loading_overlay(),
    shiny::tags$script(.loading_overlay_js()),
    shiny::div(class = "scv-content",
               shiny::uiOutput("page_content")),
    .page_footer()
  )
  server <- function(input, output, session) {
    chosen <- shiny::reactive({
      qs <- shiny::parseQueryString(session$clientData$url_search %||% "")
      sid <- qs$study
      if (is.null(sid) || !nzchar(sid)) return(NULL)
      hit <- studies[studies$studyID == sid, , drop = FALSE]
      if (nrow(hit) == 0L) return(NULL)
      tryCatch(load_study(hit$bundle_path[1], shard_dir = shard_dir),
               error = function(e) NULL)
    })
    output$page_content <- shiny::renderUI({
      s <- chosen()
      if (is.null(s)) .browser_index_content(studies)
      else .app_ui_content(s, with_back_link = TRUE)
    })
    wired <- shiny::reactiveVal(FALSE)
    shiny::observe({
      if (isTRUE(wired())) return()
      s <- chosen()
      if (is.null(s)) return()
      .app_server_logic(s, input, output, session)
      wired(TRUE)
    })
  }
  shiny::shinyApp(ui = ui, server = server)
}
