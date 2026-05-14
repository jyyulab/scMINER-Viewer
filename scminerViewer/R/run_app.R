#' Launch the scMINER Viewer Shiny app.
#'
#' Reads the bundle once at startup and serves the UI matching the layout
#' at <https://scminer.stjude.org/study/Tregs> (core plots only — see the
#' implementation plan for what's intentionally out of scope in v1).
#'
#' @param bundle_path Path to a `.scminer.h5` bundle (produced by
#'   [write_bundle()] or [prepare_study_data()]).
#' @param shard_dir Directory containing the per-gene shard tree
#'   (`expression_files/<studyID>/`, `activity_files/<studyID>/`). If
#'   `NULL` (default), uses `dirname(bundle_path)`. Pass an explicit
#'   path when the bundle is not co-located with its shards.
#' @param host Host interface to bind to. Defaults to `"127.0.0.1"`.
#' @param port Port to bind to. Defaults to a random free port (NULL).
#' @param launch_browser Logical; open the system browser. Default is
#'   `interactive()`.
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return Called for side effect (launches the Shiny app). Returns
#'   `invisible(NULL)`.
#' @export
run_app <- function(bundle_path,
                    shard_dir = NULL,
                    host = "127.0.0.1",
                    port = NULL,
                    launch_browser = interactive(),
                    ...) {
  for (pkg in c("shiny", "bslib", "plotly", "DT", "visNetwork")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf(
        "Package %s is required for run_app(). Install it with install.packages('%s').",
        pkg, pkg
      ))
    }
  }
  study <- load_study(bundle_path, shard_dir = shard_dir)
  app <- shiny::shinyApp(
    ui     = .app_ui(study),
    server = .app_server(study)
  )
  shiny::runApp(app, host = host, port = port,
                launch.browser = launch_browser, ...)
}

#' Build a Shiny app object without launching it.
#'
#' Useful for tests, embedding in larger apps, or deploying via shinyServer.
#'
#' @inheritParams run_app
#' @return A `shiny.appobj`.
#' @export
build_app <- function(bundle_path, shard_dir = NULL) {
  study <- load_study(bundle_path, shard_dir = shard_dir)
  shiny::shinyApp(
    ui     = .app_ui(study),
    server = .app_server(study)
  )
}
