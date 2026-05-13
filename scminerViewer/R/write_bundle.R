#' Write a self-contained scMINER study bundle (HDF5).
#'
#' The bundle is a single `.scminer.h5` file consumable by both the R Shiny app
#' and the companion Python webui. Sparse matrices are stored in CSR layout
#' (rows = genes, columns = cells) with zero-based indices, matching the
#' scipy.sparse convention so Python can construct `csr_matrix` directly.
#'
#' @param bundle_path Output path; should end in `.scminer.h5`.
#' @param meta Named list with `studyID`, `studyAbbr`, `longTitle`,
#'   `shortTitle`, `species`, `coordinate`.
#' @param cells data.frame with columns `cellID`, `cellType`, `cellGroup`
#'   (optional, defaults to `cellType`), `coord1`, `coord2`.
#' @param clusters data.frame with columns `cellType`, `count`, `color`
#'   (optional), `label_1` (optional), `label_2` (optional).
#' @param genes Character vector of gene symbols. Length must match
#'   `nrow(expression)`.
#' @param expression Optional `Matrix` (or coercible) sized `length(genes)` by
#'   `nrow(cells)`.
#' @param activity_tf Optional matrix in the same orientation.
#' @param activity_sig Optional matrix in the same orientation.
#' @param network_tf Optional data.frame with `source`, `target`, `cellType`,
#'   `mi`, `pearson`, `spearman`, `rho`, `pvalue`.
#' @param network_sig Optional data.frame with the same columns.
#' @param overwrite Logical; allow overwriting `bundle_path`.
#'
#' @return `bundle_path` (invisibly).
#' @export
write_bundle <- function(bundle_path,
                         meta,
                         cells,
                         clusters,
                         genes,
                         expression = NULL,
                         activity_tf = NULL,
                         activity_sig = NULL,
                         network_tf = NULL,
                         network_sig = NULL,
                         overwrite = FALSE) {
  if (file.exists(bundle_path)) {
    if (!isTRUE(overwrite)) {
      stop("Bundle file already exists: ", bundle_path,
           " (pass overwrite = TRUE to replace)")
    }
    file.remove(bundle_path)
  }
  dir.create(dirname(bundle_path), showWarnings = FALSE, recursive = TRUE)

  .validate_meta(meta)
  .validate_cells(cells)
  .validate_clusters(clusters)
  if (!is.character(genes) || length(genes) == 0) {
    stop("`genes` must be a non-empty character vector")
  }

  n_cells <- nrow(cells)
  n_genes <- length(genes)
  if (!is.null(expression))   .validate_matrix(expression,   n_genes, n_cells, "expression")
  if (!is.null(activity_tf))  .validate_matrix(activity_tf,  n_genes, n_cells, "activity_tf")
  if (!is.null(activity_sig)) .validate_matrix(activity_sig, n_genes, n_cells, "activity_sig")

  file <- hdf5r::H5File$new(bundle_path, mode = "w")
  on.exit({ try(file$close_all(), silent = TRUE) }, add = TRUE)

  .write_meta(file, meta)
  .write_cells(file, cells)
  .write_clusters(file, clusters)
  .write_genes(file, genes)
  if (!is.null(expression))   .write_sparse(file, "expression",   expression)
  if (!is.null(activity_tf))  .write_sparse(file, "activity_tf",  activity_tf)
  if (!is.null(activity_sig)) .write_sparse(file, "activity_sig", activity_sig)
  if (!is.null(network_tf))   .write_network(file, "network_tf",  network_tf)
  if (!is.null(network_sig))  .write_network(file, "network_sig", network_sig)

  invisible(bundle_path)
}

.write_meta <- function(file, meta) {
  grp <- file$create_group("meta")
  for (k in c("studyID", "studyAbbr", "longTitle", "shortTitle",
              "species", "coordinate")) {
    grp[[k]] <- as.character(meta[[k]] %||% "")
  }
  grp[["bundleVersion"]] <- .bundle_version
}

.write_cells <- function(file, cells) {
  grp <- file$create_group("cells")
  grp[["cellID"]]    <- as.character(cells$cellID)
  grp[["cellType"]]  <- as.character(cells$cellType)
  grp[["cellGroup"]] <- as.character(cells$cellGroup %||% cells$cellType)
  grp[["coord1"]]    <- as.numeric(cells$coord1)
  grp[["coord2"]]    <- as.numeric(cells$coord2)
}

.write_clusters <- function(file, clusters) {
  grp <- file$create_group("clusters")
  grp[["cellType"]] <- as.character(clusters$cellType)
  grp[["count"]]    <- as.integer(clusters$count)
  grp[["color"]]    <- as.character(clusters$color %||%
                                      rep("#888888", nrow(clusters)))
  if (!is.null(clusters$label_1)) grp[["label_1"]] <- as.numeric(clusters$label_1)
  if (!is.null(clusters$label_2)) grp[["label_2"]] <- as.numeric(clusters$label_2)
}

.write_genes <- function(file, genes) {
  grp <- file$create_group("genes")
  grp[["symbol"]] <- as.character(genes)
}

.write_sparse <- function(file, name, mat) {
  if (!inherits(mat, "Matrix")) {
    mat <- methods::as(as.matrix(mat), "CsparseMatrix")
  }
  rmat <- methods::as(mat, "RsparseMatrix")
  grp <- file$create_group(name)
  grp[["data"]]    <- as.numeric(rmat@x)
  grp[["indices"]] <- as.integer(rmat@j)
  grp[["indptr"]]  <- as.integer(rmat@p)
  grp[["shape"]]   <- as.integer(dim(rmat))
}

.write_network <- function(file, name, df) {
  required <- c("source", "target", "cellType",
                "mi", "pearson", "spearman", "rho", "pvalue")
  miss <- setdiff(required, colnames(df))
  if (length(miss) > 0) {
    stop(sprintf("network `%s` missing columns: %s",
                 name, paste(miss, collapse = ", ")))
  }
  grp <- file$create_group(name)
  grp[["source"]]   <- as.character(df$source)
  grp[["target"]]   <- as.character(df$target)
  grp[["cellType"]] <- as.character(df$cellType)
  grp[["mi"]]       <- as.numeric(df$mi)
  grp[["pearson"]]  <- as.numeric(df$pearson)
  grp[["spearman"]] <- as.numeric(df$spearman)
  grp[["rho"]]      <- as.numeric(df$rho)
  grp[["pvalue"]]   <- as.numeric(df$pvalue)
}

.validate_meta <- function(meta) {
  if (!is.list(meta)) stop("`meta` must be a list")
  required <- c("studyID", "studyAbbr", "longTitle", "shortTitle",
                "species", "coordinate")
  miss <- setdiff(required, names(meta))
  if (length(miss) > 0) {
    stop("`meta` missing fields: ", paste(miss, collapse = ", "))
  }
}

.validate_cells <- function(cells) {
  required <- c("cellID", "cellType", "coord1", "coord2")
  miss <- setdiff(required, colnames(cells))
  if (length(miss) > 0) {
    stop("`cells` missing columns: ", paste(miss, collapse = ", "))
  }
}

.validate_clusters <- function(clusters) {
  required <- c("cellType", "count")
  miss <- setdiff(required, colnames(clusters))
  if (length(miss) > 0) {
    stop("`clusters` missing columns: ", paste(miss, collapse = ", "))
  }
}

.validate_matrix <- function(mat, n_genes, n_cells, name) {
  d <- dim(mat)
  if (is.null(d) || d[1] != n_genes || d[2] != n_cells) {
    stop(sprintf("`%s` must be %d (genes) x %d (cells), got %s",
                 name, n_genes, n_cells,
                 if (is.null(d)) "<no dim>" else paste(d, collapse = "x")))
  }
}
