#' Write a self-contained scMINER study metadata bundle (HDF5).
#'
#' The bundle is a small `.scminer.h5` file containing study metadata,
#' cell + cluster info, the master gene list, networks, and per-matrix
#' **gene indexes** (which genes are available in `expression` /
#' `activity_tf` / `activity_sig`). Expression and activity matrix
#' values are *not* stored in the bundle — they are read on demand from
#' the on-disk shard tree by [gene_values()].
#'
#' @param bundle_path Output path; should end in `.scminer.h5`.
#' @param meta Named list with `studyID`, `studyAbbr`, `longTitle`,
#'   `shortTitle`, `species`, `coordinate`.
#' @param cells data.frame with columns `cellID`, `cellType`, `cellGroup`
#'   (optional, defaults to `cellType`), `coord1`, `coord2`.
#' @param clusters data.frame with columns `cellType`, `count`, `color`
#'   (optional), `label_1` (optional), `label_2` (optional).
#' @param genes Master character vector of gene symbols (all genes in
#'   the study, regardless of which matrices they appear in).
#' @param expression_genes,activity_tf_genes,activity_sig_genes Optional
#'   character vectors listing the genes that have shards in each matrix.
#'   Typically the `GeneSymbol` column of the corresponding manifest CSV.
#' @param default_genes Optional character vector; if present, the app
#'   pre-selects these genes on startup.
#' @param network_tf,network_sig Optional data.frames.
#' @param overwrite Logical; allow overwriting `bundle_path`.
#'
#' @return `bundle_path` (invisibly).
#' @export
write_bundle <- function(bundle_path,
                         meta,
                         cells,
                         clusters,
                         genes,
                         expression_genes  = NULL,
                         activity_tf_genes = NULL,
                         activity_sig_genes = NULL,
                         default_genes     = NULL,
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

  file <- hdf5r::H5File$new(bundle_path, mode = "w")
  on.exit({ try(file$close_all(), silent = TRUE) }, add = TRUE)

  .write_meta(file, meta)
  .write_cells(file, cells)
  .write_clusters(file, clusters)
  .write_genes(file, genes)
  .write_index(file, expression_genes, activity_tf_genes, activity_sig_genes)
  .write_defaults(file, default_genes)
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

.write_index <- function(file, exp_g, tf_g, sig_g) {
  if (is.null(exp_g) && is.null(tf_g) && is.null(sig_g)) return(invisible(NULL))
  grp <- file$create_group("index")
  if (!is.null(exp_g))  grp[["expression"]]   <- as.character(exp_g)
  if (!is.null(tf_g))   grp[["activity_tf"]]  <- as.character(tf_g)
  if (!is.null(sig_g))  grp[["activity_sig"]] <- as.character(sig_g)
}

.write_defaults <- function(file, default_genes) {
  if (is.null(default_genes) || length(default_genes) == 0) {
    return(invisible(NULL))
  }
  grp <- file$create_group("defaults")
  grp[["genes"]] <- as.character(default_genes)
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
