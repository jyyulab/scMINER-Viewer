#' Path to the demo scMINER bundle shipped with the package.
#'
#' Returns the on-disk path to `2327.scminer.h5` (the Tex CD8+ T-cell
#' study) under `inst/extdata/`, together with its co-located shard
#' tree. Pass the returned path straight to [load_study()] or
#' [run_app()].
#'
#' Only ~200 gene shards are bundled (a curated mix of canonical
#' T-cell / exhaustion markers plus a fill from each modality's index):
#' ~67 expression, ~67 SIG activity, ~66 TF activity. See
#' `inst/extdata/demo_genes.csv` for the manifest. Selecting a gene
#' outside that list in the app yields an empty plot.
#'
#' @return Absolute path to the demo `.scminer.h5` bundle.
#' @export
#' @examples
#' \dontrun{
#'   study <- load_study(demo_bundle_path())
#'   print(study)
#'   run_demo()
#' }
demo_bundle_path <- function() {
  path <- system.file("extdata", "2327.scminer.h5",
                      package = "scminerViewer")
  if (!nzchar(path)) {
    stop("Demo bundle not found. Reinstall scminerViewer, or run ",
         "`Rscript inst/scripts/build_demo_data.R` from the repo root.")
  }
  path
}

#' Manifest of genes that have shards in the bundled demo.
#'
#' Returns a `data.frame` with columns `modality` (one of
#' `"expression"`, `"activity_sig"`, `"activity_tf"`) and `gene`. Use
#' this to drive demo screenshots or to filter inputs to the curated
#' set when scripting against the demo bundle.
#'
#' @return A `data.frame` with `modality` and `gene` columns.
#' @export
demo_genes <- function() {
  path <- system.file("extdata", "demo_genes.csv",
                      package = "scminerViewer")
  if (!nzchar(path)) {
    stop("demo_genes.csv not found in package; reinstall scminerViewer.")
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

#' Launch the Shiny app on the bundled demo study.
#'
#' Convenience wrapper around `run_app(demo_bundle_path())` — the
#' fastest way to see what the viewer looks like without preparing
#' your own study data.
#'
#' @inheritParams run_app
#' @return Called for side effect; invisibly returns `NULL`.
#' @export
run_demo <- function(host = "127.0.0.1",
                     port = NULL,
                     launch_browser = interactive(),
                     ...) {
  run_app(bundle_path    = demo_bundle_path(),
          host           = host,
          port           = port,
          launch_browser = launch_browser,
          ...)
}
