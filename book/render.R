#!/usr/bin/env Rscript
# Render the scMINER Viewer Guide.
#
# Usage (from the project root OR the book/ dir):
#   Rscript book/render.R              # gitbook only
#   Rscript book/render.R pdf          # pdf only (needs xelatex)
#   Rscript book/render.R all          # both
#
# Output goes to ../docs/ by default (per output_dir in _bookdown.yml),
# so the rendered site lives at <project_root>/docs/ and is ready to
# serve from GitHub Pages (Settings → Pages → main → /docs).

if (basename(getwd()) != "book") setwd("book")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[1] else "gitbook"

# Resolve the output directory from _bookdown.yml so we can print where
# the rendered files actually landed (the bookdown::render_book function
# returns only the index path).
.resolve_output_dir <- function() {
  out <- "_book"
  if (requireNamespace("yaml", quietly = TRUE) &&
      file.exists("_bookdown.yml")) {
    cfg <- tryCatch(yaml::read_yaml("_bookdown.yml"),
                    error = function(e) NULL)
    if (!is.null(cfg) && !is.null(cfg$output_dir)) out <- cfg$output_dir
  }
  normalizePath(out, mustWork = FALSE)
}

if (mode %in% c("gitbook", "all")) {
  bookdown::render_book("index.Rmd",
                        output_format = "bookdown::gitbook")
}
if (mode %in% c("pdf", "all")) {
  bookdown::render_book("index.Rmd",
                        output_format = "bookdown::pdf_book")
}

out_dir <- .resolve_output_dir()
cat("\nOutput in: ", out_dir, "\n", sep = "")
if (dir.exists(out_dir)) {
  idx <- file.path(out_dir, "index.html")
  if (file.exists(idx)) {
    cat("Open: file://", normalizePath(idx, mustWork = TRUE), "\n",
        sep = "")
  }
}
