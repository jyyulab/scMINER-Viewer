#' Load a scMINER study bundle from HDF5.
#'
#' Returns an S3 `scminer_study` object with study metadata, cell +
#' cluster info, the master gene list, per-matrix gene indexes, networks,
#' and a path to the shard tree from which [gene_values()] reads
#' individual gene shards on demand.
#'
#' @param bundle_path Path to a `.scminer.h5` file written by
#'   [write_bundle()].
#' @param shard_dir Directory containing the `expression_files/<studyID>/`
#'   and `activity_files/<studyID>/` shard trees. Defaults to
#'   `dirname(bundle_path)` (i.e. the bundle is co-located with its
#'   shards).
#'
#' @return An S3 object of class `scminer_study`.
#' @export
load_study <- function(bundle_path, shard_dir = NULL) {
  if (!file.exists(bundle_path)) {
    stop("Bundle not found: ", bundle_path)
  }
  if (is.null(shard_dir)) shard_dir <- dirname(bundle_path)
  shard_dir <- normalizePath(shard_dir, mustWork = FALSE)

  file <- hdf5r::H5File$new(bundle_path, mode = "r")
  on.exit({ try(file$close_all(), silent = TRUE) }, add = TRUE)

  meta     <- .read_meta(file)
  cells    <- .read_cells(file)
  clusters <- .read_clusters(file)
  genes    <- file[["genes/symbol"]]$read()

  expression_index   <- .read_optional(file, "index/expression")
  activity_tf_index  <- .read_optional(file, "index/activity_tf")
  activity_sig_index <- .read_optional(file, "index/activity_sig")
  default_genes      <- .read_optional(file, "defaults/genes")

  network_tf  <- if (file$exists("network_tf"))  .read_network(file, "network_tf")  else NULL
  network_sig <- if (file$exists("network_sig")) .read_network(file, "network_sig") else NULL

  structure(
    list(
      meta               = meta,
      cells              = cells,
      clusters           = clusters,
      genes              = genes,
      expression_index   = expression_index,
      activity_tf_index  = activity_tf_index,
      activity_sig_index = activity_sig_index,
      default_genes      = default_genes,
      network_tf         = network_tf,
      network_sig        = network_sig,
      bundle_path        = normalizePath(bundle_path, mustWork = FALSE),
      shard_dir          = shard_dir,
      .cache             = new.env(parent = emptyenv()),
      .meta_cache        = new.env(parent = emptyenv())
    ),
    class = "scminer_study"
  )
}

#' Read one gene's values from the shard tree.
#'
#' Lazily reads `<gene>.csv.gz` from the on-disk shard tree, parses it,
#' and aligns the values to `study$cells$cellID` (cells missing from the
#' shard header are returned as `NA`). Results are cached on the study
#' object so repeated calls for the same gene are free.
#'
#' @param study An `scminer_study` returned by [load_study()].
#' @param gene Gene symbol.
#' @param relationship One of `"Express_normalized"` (expression),
#'   `"Activity_tf"`, `"Activity_sig"`.
#'
#' @return Numeric vector of length `nrow(study$cells)`, or `NULL` if
#'   the gene isn't available for the requested relationship.
#' @export
gene_values <- function(study, gene, relationship = "Express_normalized") {
  stopifnot(inherits(study, "scminer_study"))
  if (length(gene) != 1L || !nzchar(gene)) return(NULL)
  spec <- .shard_spec(study, relationship)
  if (is.null(spec)) return(NULL)
  if (!is.null(spec$index) && !(gene %in% spec$index)) return(NULL)

  key <- paste(spec$tag, gene, sep = ":")
  cached <- study$.cache[[key]]
  if (!is.null(cached)) return(cached)

  perm_info <- .shard_perm(study, spec)
  if (is.null(perm_info)) return(NULL)

  letter <- .shard_letter(gene)
  shard_name <- gsub("/", "_", gene, fixed = TRUE)
  shard_path <- file.path(study$shard_dir, spec$shard_subdir,
                          letter, paste0(shard_name, ".csv.gz"))
  if (!file.exists(shard_path)) return(NULL)

  vals <- tryCatch(
    data.table::fread(shard_path, header = FALSE, sep = ",",
                      data.table = FALSE),
    error = function(e) NULL
  )
  if (is.null(vals) || nrow(vals) == 0 || ncol(vals) == 0) return(NULL)
  shard_vals <- as.numeric(vals[1, ])
  if (length(shard_vals) != perm_info$n_shard_cells) {
    warning(sprintf(
      "[%s] gene %s: shard has %d values but meta has %d cells",
      relationship, gene, length(shard_vals), perm_info$n_shard_cells),
      call. = FALSE)
    return(NULL)
  }
  aligned <- shard_vals[perm_info$perm]
  study$.cache[[key]] <- aligned
  aligned
}

.shard_spec <- function(study, relationship) {
  study_id <- as.character(study$meta$studyID)
  switch(relationship,
    "Express_normalized" = list(
      tag           = "exp",
      shard_subdir  = file.path("expression_files", study_id),
      meta_path     = file.path(study$shard_dir, "expression_files",
                                study_id, "meta.csv"),
      index         = study$expression_index
    ),
    "Activity_tf" = list(
      tag           = "tf",
      shard_subdir  = file.path("activity_files", study_id, "TF"),
      meta_path     = file.path(study$shard_dir, "activity_files",
                                study_id, "meta.csv"),
      index         = study$activity_tf_index
    ),
    "Activity_sig" = list(
      tag           = "sig",
      shard_subdir  = file.path("activity_files", study_id, "SIG"),
      meta_path     = file.path(study$shard_dir, "activity_files",
                                study_id, "meta.csv"),
      index         = study$activity_sig_index
    ),
    NULL
  )
}

.shard_perm <- function(study, spec) {
  cached <- study$.meta_cache[[spec$tag]]
  if (!is.null(cached)) return(cached)
  if (!file.exists(spec$meta_path)) {
    warning(sprintf("meta.csv not found: %s", spec$meta_path),
            call. = FALSE)
    study$.meta_cache[[spec$tag]] <- NULL
    return(NULL)
  }
  header_line <- tryCatch(readLines(spec$meta_path, n = 1, warn = FALSE),
                          error = function(e) character(0))
  if (length(header_line) == 0 || !nzchar(header_line)) {
    warning(sprintf("meta.csv is empty: %s", spec$meta_path),
            call. = FALSE)
    return(NULL)
  }
  shard_cells <- trimws(strsplit(header_line, ",", fixed = TRUE)[[1]])
  out <- list(
    perm          = match(study$cells$cellID, shard_cells),
    n_shard_cells = length(shard_cells)
  )
  study$.meta_cache[[spec$tag]] <- out
  out
}

.read_optional <- function(file, path) {
  # hdf5r's exists() raises when the parent group is missing; walk the
  # path manually so we can return NULL silently.
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  node <- file
  for (p in parts) {
    if (!node$exists(p)) return(NULL)
    node <- node[[p]]
  }
  node$read()
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
  count <- function(v) if (is.null(v)) "-" else as.character(length(v))
  cat("<scminer_study>\n")
  cat("  Study:      ", x$meta$studyAbbr, " - ", x$meta$shortTitle, "\n", sep = "")
  cat("  ID:         ", x$meta$studyID,   "\n", sep = "")
  cat("  Species:    ", x$meta$species,    "\n", sep = "")
  cat("  Coordinate: ", x$meta$coordinate, "\n", sep = "")
  cat("  Cells:      ", nrow(x$cells),     "\n", sep = "")
  cat("  Genes:      ", length(x$genes),   "\n", sep = "")
  cat("  Clusters:   ", nrow(x$clusters),  "\n", sep = "")
  cat("  Expression: ", count(x$expression_index),   " genes\n", sep = "")
  cat("  ActivityTF: ", count(x$activity_tf_index),  " genes\n", sep = "")
  cat("  ActivitySIG:", count(x$activity_sig_index), " genes\n", sep = "")
  cat("  Defaults:   ", count(x$default_genes),      " genes\n", sep = "")
  cat("  NetworkTF:  ", if (!is.null(x$network_tf))  nrow(x$network_tf)  else "-", " edges\n", sep = "")
  cat("  NetworkSIG: ", if (!is.null(x$network_sig)) nrow(x$network_sig) else "-", " edges\n", sep = "")
  cat("  Shard dir:  ", x$shard_dir, "\n", sep = "")
  invisible(x)
}
