#' Prepare an scMINER study: emit graph-import layout and/or HDF5 bundle.
#'
#' Three entry points:
#'
#' \itemize{
#'   \item [prepare_study_data()] — lowest-level form; takes already-extracted
#'     `meta`, `cells`, `clusters`, `genes`, sparse matrices, and network
#'     data frames. No Biobase / yaml dependency.
#'   \item [prepare_study_from_eset()] — accepts a Biobase ExpressionSet
#'     (one for expression, optional one for activity) and pulls out
#'     pData / fData / exprs internally. Requires the `Biobase` package
#'     (suggested).
#'   \item [prepare_study()] — drop-in for the original `main.R`: reads a
#'     YAML config, loads the referenced RDS / TSV files, then calls
#'     `prepare_study_from_eset`. Requires the `yaml` and `Biobase`
#'     packages.
#' }
#'
#' All three accept `emit = c("graph", "bundle")`; either or both. The
#' bundle is always written as `<out_dir>/<studyID>.scminer.h5`.
#'
#' @name prepare_study
NULL

#' @rdname prepare_study
#' @param out_dir Output root directory.
#' @param meta Named list with `studyID`, `studyAbbr`, `longTitle`,
#'   `shortTitle`, `species`, `coordinate`.
#' @param cells data.frame: `cellID`, `cellType`, `cellGroup` (optional),
#'   `coord1`, `coord2`.
#' @param clusters data.frame: `cellType`, `count` (optional, computed
#'   from `cells` if absent), `color`, `label_1`, `label_2`.
#' @param genes Character vector of gene symbols.
#' @param expression Optional `Matrix` (G x N).
#' @param activity_tf,activity_sig Optional matrices, same orientation.
#' @param network_tf,network_sig Optional data.frames.
#' @param emit Character vector — any subset of `c("graph", "bundle")`.
#' @param verbose Logical; emit progress messages while writing shards.
#' @export
prepare_study_data <- function(out_dir,
                                meta,
                                cells,
                                clusters = NULL,
                                genes,
                                expression = NULL,
                                activity_tf = NULL,
                                activity_sig = NULL,
                                network_tf = NULL,
                                network_sig = NULL,
                                emit = c("graph", "bundle"),
                                verbose = FALSE) {
  emit <- match.arg(emit, c("graph", "bundle"), several.ok = TRUE)
  stopifnot(is.list(meta), is.data.frame(cells), is.character(genes))

  if (is.null(clusters) || nrow(clusters) == 0) {
    cnt <- as.data.frame(table(cells$cellType), stringsAsFactors = FALSE)
    clusters <- data.frame(
      cellType = cnt$Var1,
      count    = as.integer(cnt$Freq),
      color    = rep("#888888", nrow(cnt)),
      stringsAsFactors = FALSE
    )
  } else if (is.null(clusters$count) || all(is.na(clusters$count))) {
    cnt <- as.data.frame(table(cells$cellType), stringsAsFactors = FALSE)
    clusters$count <- cnt$Freq[match(clusters$cellType, cnt$Var1)]
    clusters$count[is.na(clusters$count)] <- 0L
  }

  # Align matrix row names to the master gene list so shard writes can
  # use them directly. Callers can pass matrices without rownames as long
  # as the row order matches `genes`.
  attach_rownames <- function(mat, name) {
    if (is.null(mat)) return(NULL)
    if (nrow(mat) != length(genes)) {
      stop(sprintf("%s row count (%d) does not match length(genes) (%d)",
                   name, nrow(mat), length(genes)))
    }
    rownames(mat) <- genes
    mat
  }
  expression   <- attach_rownames(expression,   "expression")
  activity_tf  <- attach_rownames(activity_tf,  "activity_tf")
  activity_sig <- attach_rownames(activity_sig, "activity_sig")

  if ("graph" %in% emit) {
    if (isTRUE(verbose)) message("Writing graph-import layout to ", out_dir)
    .ensure_graph_tree(out_dir)
    .write_graph_study(out_dir, meta)
    .write_graph_clusters(out_dir, meta, clusters)
    .write_graph_genes(out_dir, meta, genes)
    .write_graph_cells(out_dir, meta, cells)
    .write_graph_networks(out_dir, meta, network_tf, network_sig)

    study_id <- as.character(meta$studyID)
    exp_root <- file.path("expression_files", study_id)
    act_root <- file.path("activity_files",   study_id)

    if (!is.null(expression)) {
      .write_graph_shards(out_dir, meta, expression,
                          kind          = exp_root,
                          meta_kind     = exp_root,
                          manifest_dir  = "study_gene_expression",
                          manifest_name = "expression",
                          cell_ids      = cells$cellID,
                          type_label    = "Expression",
                          verbose       = verbose)
    }
    if (!is.null(activity_tf)) {
      .write_graph_shards(out_dir, meta, activity_tf,
                          kind          = file.path(act_root, "TF"),
                          meta_kind     = act_root,
                          manifest_dir  = "study_gene_tf",
                          manifest_name = "activity_tf",
                          cell_ids      = cells$cellID,
                          type_label    = "TF",
                          verbose       = verbose)
    }
    if (!is.null(activity_sig)) {
      .write_graph_shards(out_dir, meta, activity_sig,
                          kind          = file.path(act_root, "SIG"),
                          meta_kind     = act_root,
                          manifest_dir  = "study_gene_sig",
                          manifest_name = "activity_sig",
                          cell_ids      = cells$cellID,
                          type_label    = "SIG",
                          verbose       = verbose)
    }
  }

  bundle_path <- NULL
  if ("bundle" %in% emit) {
    bundle_path <- file.path(out_dir,
                             paste0(meta$studyID, ".scminer.h5"))
    if (isTRUE(verbose)) message("Writing bundle to ", bundle_path)
    write_bundle(
      bundle_path  = bundle_path,
      meta         = meta,
      cells        = cells,
      clusters     = clusters,
      genes        = genes,
      expression   = expression,
      activity_tf  = activity_tf,
      activity_sig = activity_sig,
      network_tf   = network_tf,
      network_sig  = network_sig,
      overwrite    = TRUE
    )
  }

  invisible(list(out_dir = out_dir, bundle_path = bundle_path))
}

#' @rdname prepare_study
#' @param expression_eset A Biobase ExpressionSet for the expression matrix.
#'   `pData(eset)` provides the cell metadata; `fData(eset)` provides gene
#'   symbols; `exprs(eset)` provides the value matrix.
#' @param activity_eset Optional Biobase ExpressionSet for the activity
#'   matrix. Row names must end in `_TF` / `.TF` or `_SIG` / `.SIG`; rows
#'   are split into `activity_tf` and `activity_sig` accordingly.
#' @param networks_path Optional path to a tab-separated networks file with
#'   columns `source, target, NetworkType, [studyID,] CellGroup, mi,
#'   pearson, spearman, rho, pvalue` (the format consumed by the original
#'   `h_networks.R`).
#' @param cell_id_col,cell_type_col,cell_group_col,coordinate_col
#'   Column names within `pData(expression_eset)`.
#' @param gene_symbol_col Column name within `fData(expression_eset)`.
#' @export
prepare_study_from_eset <- function(out_dir,
                                     expression_eset,
                                     activity_eset = NULL,
                                     networks_path = NULL,
                                     meta,
                                     cell_id_col      = "cellID",
                                     cell_type_col    = "cellGroup",
                                     cell_group_col   = "cellGroup",
                                     coordinate_col   = "UMAP",
                                     gene_symbol_col  = "geneSymbol",
                                     clusters         = NULL,
                                     emit             = c("graph", "bundle"),
                                     verbose          = FALSE) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase is required for prepare_study_from_eset(); ",
         "install it from Bioconductor.")
  }

  p_data <- as.data.frame(Biobase::pData(expression_eset),
                          stringsAsFactors = FALSE)
  p_data[] <- lapply(p_data, function(x) {
    if (is.factor(x)) as.character(x) else x
  })
  p_data[[cell_id_col]] <- rownames(p_data)

  coord1_col <- paste0(coordinate_col, "_1")
  coord2_col <- paste0(coordinate_col, "_2")
  required <- c(cell_id_col, cell_type_col, coord1_col, coord2_col)
  miss <- setdiff(required, colnames(p_data))
  if (length(miss) > 0) {
    stop("pData missing columns: ", paste(miss, collapse = ", "))
  }
  cells <- data.frame(
    cellID    = as.character(p_data[[cell_id_col]]),
    cellType  = as.character(p_data[[cell_type_col]]),
    cellGroup = as.character(p_data[[cell_group_col]] %||%
                              p_data[[cell_type_col]]),
    coord1    = as.numeric(p_data[[coord1_col]]),
    coord2    = as.numeric(p_data[[coord2_col]]),
    stringsAsFactors = FALSE
  )
  meta$coordinate <- meta$coordinate %||% coordinate_col

  f_data <- as.data.frame(Biobase::fData(expression_eset),
                          stringsAsFactors = FALSE)
  if (!gene_symbol_col %in% colnames(f_data)) {
    stop("fData missing gene-symbol column: ", gene_symbol_col)
  }
  genes <- as.character(f_data[[gene_symbol_col]])

  expr <- Biobase::exprs(expression_eset)
  if (!inherits(expr, "Matrix")) {
    expr <- methods::as(as.matrix(expr), "CsparseMatrix")
  }
  rownames(expr) <- genes

  activity_tf  <- NULL
  activity_sig <- NULL
  if (!is.null(activity_eset)) {
    act <- Biobase::exprs(activity_eset)
    act_rows <- rownames(act)
    is_tf  <- grepl("[._]TF$",  act_rows)
    is_sig <- grepl("[._]SIG$", act_rows)
    if (!any(is_tf) && !any(is_sig)) {
      warning("activity_eset has no rows ending in _TF/.TF or _SIG/.SIG; ",
              "skipping activity matrices.")
    }
    strip_suffix <- function(rows) {
      sub("[._](TF|SIG)$", "", rows)
    }
    if (any(is_tf)) {
      sub_act <- act[is_tf, , drop = FALSE]
      rownames(sub_act) <- strip_suffix(rownames(sub_act))
      if (!inherits(sub_act, "Matrix")) {
        sub_act <- methods::as(as.matrix(sub_act), "CsparseMatrix")
      }
      # Reindex rows to the master gene list (zero rows for genes absent
      # from the activity matrix)
      activity_tf <- .reindex_rows(sub_act, genes)
    }
    if (any(is_sig)) {
      sub_act <- act[is_sig, , drop = FALSE]
      rownames(sub_act) <- strip_suffix(rownames(sub_act))
      if (!inherits(sub_act, "Matrix")) {
        sub_act <- methods::as(as.matrix(sub_act), "CsparseMatrix")
      }
      activity_sig <- .reindex_rows(sub_act, genes)
    }
  }

  networks <- NULL
  if (!is.null(networks_path)) {
    networks <- .read_networks_file(networks_path)
  }
  network_tf  <- networks$tf
  network_sig <- networks$sig

  prepare_study_data(
    out_dir      = out_dir,
    meta         = meta,
    cells        = cells,
    clusters     = clusters,
    genes        = genes,
    expression   = expr,
    activity_tf  = activity_tf,
    activity_sig = activity_sig,
    network_tf   = network_tf,
    network_sig  = network_sig,
    emit         = emit,
    verbose      = verbose
  )
}

#' @rdname prepare_study
#' @param config_path Path to a YAML config compatible with the original
#'   `scMINER-portal-datapre-R/main.R`.
#' @export
prepare_study <- function(config_path,
                           emit = c("graph", "bundle"),
                           verbose = FALSE) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("yaml is required for prepare_study(); install.packages('yaml').")
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase is required for prepare_study(); ",
         "install it from Bioconductor.")
  }
  cfg <- yaml::yaml.load_file(config_path)
  study <- cfg$study

  out_dir <- cfg$output
  meta <- list(
    studyID    = as.character(study$ID),
    studyAbbr  = as.character(study$studyAbbr),
    longTitle  = as.character(study$longTitle),
    shortTitle = as.character(study$shortTitle),
    species    = as.character(cfg$species %||% ""),
    coordinate = as.character(cfg$coordinate %||% "UMAP")
  )

  expression_eset <- readRDS(cfg$input$expression)
  activity_eset <- if (!is.null(cfg$input$activity)) {
    readRDS(cfg$input$activity)
  } else NULL

  prepare_study_from_eset(
    out_dir          = out_dir,
    expression_eset  = expression_eset,
    activity_eset    = activity_eset,
    networks_path    = cfg$input$networks,
    meta             = meta,
    cell_id_col      = cfg$cellID    %||% "cellID",
    cell_type_col    = cfg$cellType  %||% "cellGroup",
    cell_group_col   = cfg$cellGroup %||% (cfg$cellType %||% "cellGroup"),
    coordinate_col   = cfg$coordinate %||% "UMAP",
    gene_symbol_col  = cfg$geneSymbol %||% "geneSymbol",
    emit             = emit,
    verbose          = verbose
  )
}

# Reindex a matrix's rows onto a master gene list, dropping rows whose
# names aren't in the master and inserting zero rows for missing genes.
.reindex_rows <- function(mat, master_genes) {
  src <- rownames(mat)
  idx <- match(master_genes, src)
  # Build an empty CsparseMatrix and copy in
  G <- length(master_genes)
  N <- ncol(mat)
  out <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                              x = numeric(0), dims = c(G, N),
                              repr = "C")
  rownames(out) <- master_genes
  colnames(out) <- colnames(mat)
  present <- which(!is.na(idx))
  if (length(present) > 0) {
    out[present, ] <- mat[idx[present], , drop = FALSE]
  }
  out
}

# Read the networks TSV used by the original h_networks.R: header line +
# rows with columns source, target, NetworkType, [studyID,] CellGroup,
# mi, pearson, spearman, rho, pvalue.
.read_networks_file <- function(path) {
  df <- utils::read.table(path, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE,
                          comment.char = "", quote = "")
  cols <- colnames(df)
  src_col <- cols[1]
  tgt_col <- cols[2]
  if (!"NetworkType" %in% cols) {
    stop("Networks file ", path, " missing NetworkType column")
  }
  if (!"CellGroup" %in% cols) {
    stop("Networks file ", path, " missing CellGroup column")
  }
  # Score columns: take the four numeric cols after CellGroup if they
  # match mi/pearson/spearman/rho/pvalue; otherwise positional.
  score_cols <- c("mi", "pearson", "spearman", "rho", "pvalue")
  missing <- setdiff(score_cols, cols)
  if (length(missing) > 0) {
    # Fall back to positional: 5..9 of the original h_networks layout
    df_pos <- utils::read.table(path, header = TRUE, sep = "\t",
                                stringsAsFactors = FALSE,
                                comment.char = "", quote = "")
    if (ncol(df_pos) >= 9) {
      df$mi       <- df_pos[[5]]
      df$pearson  <- df_pos[[6]]
      df$spearman <- df_pos[[7]]
      df$rho      <- df_pos[[8]]
      df$pvalue   <- df_pos[[9]]
    } else {
      stop("Networks file ", path,
           " is missing required score columns: ",
           paste(missing, collapse = ", "))
    }
  }
  split_df <- split(df, df$NetworkType)
  to_canon <- function(sub) {
    if (is.null(sub) || nrow(sub) == 0) return(NULL)
    data.frame(
      source   = as.character(sub[[src_col]]),
      target   = as.character(sub[[tgt_col]]),
      cellType = as.character(sub$CellGroup),
      mi       = as.numeric(sub$mi),
      pearson  = as.numeric(sub$pearson),
      spearman = as.numeric(sub$spearman),
      rho      = as.numeric(sub$rho),
      pvalue   = as.numeric(sub$pvalue),
      stringsAsFactors = FALSE
    )
  }
  list(tf  = to_canon(split_df[["TF"]]),
       sig = to_canon(split_df[["SIG"]]))
}
