#' Load a scMINER study bundle from HDF5.
#'
#' Returns a list with class `scminer_study` containing all groups present in
#' the bundle. Optional groups (expression / activity / network) are returned
#' as `NULL` if absent.
#'
#' @param bundle_path Path to a `.scminer.h5` file written by [write_bundle()].
#'
#' @return An object of class `scminer_study` (a list).
#' @export
load_study <- function(bundle_path) {
  if (!file.exists(bundle_path)) {
    stop("Bundle not found: ", bundle_path)
  }
  file <- hdf5r::H5File$new(bundle_path, mode = "r")
  on.exit({ try(file$close_all(), silent = TRUE) }, add = TRUE)

  meta     <- .read_meta(file)
  cells    <- .read_cells(file)
  clusters <- .read_clusters(file)
  genes    <- file[["genes/symbol"]]$read()

  expression   <- if (file$exists("expression"))   .read_sparse(file, "expression")   else NULL
  activity_tf  <- if (file$exists("activity_tf"))  .read_sparse(file, "activity_tf")  else NULL
  activity_sig <- if (file$exists("activity_sig")) .read_sparse(file, "activity_sig") else NULL
  network_tf   <- if (file$exists("network_tf"))   .read_network(file, "network_tf")  else NULL
  network_sig  <- if (file$exists("network_sig"))  .read_network(file, "network_sig") else NULL

  structure(
    list(
      meta = meta,
      cells = cells,
      clusters = clusters,
      genes = genes,
      expression = expression,
      activity_tf = activity_tf,
      activity_sig = activity_sig,
      network_tf = network_tf,
      network_sig = network_sig,
      bundle_path = normalizePath(bundle_path, mustWork = FALSE)
    ),
    class = "scminer_study"
  )
}

.read_meta <- function(file) {
  grp <- file[["meta"]]
  list(
    studyID       = grp[["studyID"]]$read(),
    studyAbbr     = grp[["studyAbbr"]]$read(),
    longTitle     = grp[["longTitle"]]$read(),
    shortTitle    = grp[["shortTitle"]]$read(),
    species       = grp[["species"]]$read(),
    coordinate    = grp[["coordinate"]]$read(),
    bundleVersion = grp[["bundleVersion"]]$read()
  )
}

.read_cells <- function(file) {
  grp <- file[["cells"]]
  data.frame(
    cellID    = grp[["cellID"]]$read(),
    cellType  = grp[["cellType"]]$read(),
    cellGroup = grp[["cellGroup"]]$read(),
    coord1    = grp[["coord1"]]$read(),
    coord2    = grp[["coord2"]]$read(),
    stringsAsFactors = FALSE
  )
}

.read_clusters <- function(file) {
  grp <- file[["clusters"]]
  out <- data.frame(
    cellType = grp[["cellType"]]$read(),
    count    = grp[["count"]]$read(),
    color    = grp[["color"]]$read(),
    stringsAsFactors = FALSE
  )
  if (grp$exists("label_1")) out$label_1 <- grp[["label_1"]]$read()
  if (grp$exists("label_2")) out$label_2 <- grp[["label_2"]]$read()
  out
}

.read_sparse <- function(file, name) {
  grp <- file[[name]]
  data    <- grp[["data"]]$read()
  indices <- grp[["indices"]]$read()
  indptr  <- grp[["indptr"]]$read()
  shape   <- grp[["shape"]]$read()
  G <- as.integer(shape[1])
  N <- as.integer(shape[2])
  if (length(data) == 0L) {
    return(Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(G, N), repr = "C"
    ))
  }
  row_lengths <- diff(indptr)
  row_idx <- rep.int(seq_len(G), as.integer(row_lengths))
  Matrix::sparseMatrix(
    i = row_idx,
    j = as.integer(indices) + 1L,
    x = as.numeric(data),
    dims = c(G, N),
    repr = "C"
  )
}

.read_network <- function(file, name) {
  grp <- file[[name]]
  data.frame(
    source   = grp[["source"]]$read(),
    target   = grp[["target"]]$read(),
    cellType = grp[["cellType"]]$read(),
    mi       = grp[["mi"]]$read(),
    pearson  = grp[["pearson"]]$read(),
    spearman = grp[["spearman"]]$read(),
    rho      = grp[["rho"]]$read(),
    pvalue   = grp[["pvalue"]]$read(),
    stringsAsFactors = FALSE
  )
}

#' @export
print.scminer_study <- function(x, ...) {
  cat("<scminer_study>\n")
  cat("  Study:      ", x$meta$studyAbbr, " - ", x$meta$shortTitle, "\n", sep = "")
  cat("  ID:         ", x$meta$studyID,   "\n", sep = "")
  cat("  Species:    ", x$meta$species,    "\n", sep = "")
  cat("  Coordinate: ", x$meta$coordinate, "\n", sep = "")
  cat("  Cells:      ", nrow(x$cells),     "\n", sep = "")
  cat("  Genes:      ", length(x$genes),   "\n", sep = "")
  cat("  Clusters:   ", nrow(x$clusters),  "\n", sep = "")
  mark <- function(v) if (!is.null(v)) "yes" else "-"
  cat("  Expression: ", mark(x$expression),   "\n", sep = "")
  cat("  ActivityTF: ", mark(x$activity_tf),  "\n", sep = "")
  cat("  ActivitySIG:", mark(x$activity_sig), "\n", sep = "")
  cat("  NetworkTF:  ", if (!is.null(x$network_tf))  nrow(x$network_tf)  else "-", "\n", sep = "")
  cat("  NetworkSIG: ", if (!is.null(x$network_sig)) nrow(x$network_sig) else "-", "\n", sep = "")
  invisible(x)
}
