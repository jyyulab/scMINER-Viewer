#!/usr/bin/env Rscript
# Render the scMINER Viewer Guide to book/_book/ (gitbook + pdf).
#
# Usage (from the project root or the book/ dir):
#   Rscript book/render.R              # gitbook only
#   Rscript book/render.R pdf          # pdf only (needs xelatex)
#   Rscript book/render.R all          # both
setwd(file.path(getwd(), if (basename(getwd()) == "book") "" else "book"))
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[1] else "gitbook"

if (mode %in% c("gitbook", "all")) {
  bookdown::render_book("index.Rmd",
                        output_format = "bookdown::gitbook")
}
if (mode %in% c("pdf", "all")) {
  bookdown::render_book("index.Rmd",
                        output_format = "bookdown::pdf_book")
}
cat("Output in: ", normalizePath("_book"), "\n", sep = "")
