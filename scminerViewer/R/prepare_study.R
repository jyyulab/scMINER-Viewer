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
                                default_genes = NULL,
                                cluster_palette = "npg",
                                emit = c("graph", "bundle"),
                                verbose = FALSE) {
  emit <- match.arg(emit, c("graph", "bundle"), several.ok = TRUE)
  stopifnot(is.list(meta), is.data.frame(cells), is.character(genes))

  # Fill in counts, colours (ggsci palette per `cluster_palette`), and
  # label centroids (mean coord1/coord2 per cellType) for any clusters
  # column that's missing — see fill_clusters().
  clusters <- fill_clusters(cells, clusters, palette = cluster_palette)

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

  # Every study gets its own directory inside out_dir, so multiple
  # studies under one root never collide (and run_browser() can find
  # them via the <studyID>/<studyID>.scminer.h5 convention).
  study_id  <- as.character(meta$studyID)
  study_out <- file.path(out_dir, study_id)
  dir.create(study_out, showWarnings = FALSE, recursive = TRUE)

  if ("graph" %in% emit) {
    if (isTRUE(verbose)) message("Writing graph-import layout to ", study_out)
    .ensure_graph_tree(study_out)
    .write_graph_study(study_out, meta)
    .write_graph_clusters(study_out, meta, clusters)
    .write_graph_genes(study_out, meta, genes)
    .write_graph_cells(study_out, meta, cells)
    .write_graph_networks(study_out, meta, network_tf, network_sig)

    exp_root <- file.path("expression_files", study_id)
    act_root <- file.path("activity_files",   study_id)

    if (!is.null(expression)) {
      .write_graph_shards(study_out, meta, expression,
                          kind          = exp_root,
                          meta_kind     = exp_root,
                          manifest_dir  = "study_gene_expression",
                          manifest_name = "expression",
                          cell_ids      = cells$cellID,
                          type_label    = "Expression",
                          verbose       = verbose)
    }
    if (!is.null(activity_tf)) {
      .write_graph_shards(study_out, meta, activity_tf,
                          kind          = file.path(act_root, "TF"),
                          meta_kind     = act_root,
                          manifest_dir  = "study_gene_tf",
                          manifest_name = "activity_tf",
                          cell_ids      = cells$cellID,
                          type_label    = "TF",
                          verbose       = verbose)
    }
    if (!is.null(activity_sig)) {
      .write_graph_shards(study_out, meta, activity_sig,
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
    bundle_path <- file.path(study_out,
                             paste0(study_id, ".scminer.h5"))
    if (isTRUE(verbose)) message("Writing bundle to ", bundle_path)
    write_bundle(
      bundle_path        = bundle_path,
      meta               = meta,
      cells              = cells,
      clusters           = clusters,
      genes              = genes,
      expression_genes   = if (!is.null(expression))   rownames(expression)   else NULL,
      activity_tf_genes  = if (!is.null(activity_tf))  rownames(activity_tf)  else NULL,
      activity_sig_genes = if (!is.null(activity_sig)) rownames(activity_sig) else NULL,
      default_genes      = default_genes,
      network_tf         = network_tf,
      network_sig        = network_sig,
      overwrite          = TRUE
    )
  }

  invisible(list(
    out_dir     = study_out,    # the per-study subdir actually written to
    root_dir    = out_dir,      # the root the user passed
    bundle_path = bundle_path
  ))
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
#'
#' @details Internally, `prepare_study_from_eset()` is a thin wrapper
#' over the staged helpers [extract_cells()], [extract_genes()],
#' [extract_expression()], [extract_activity()] and [read_networks()].
#' Call those directly when you want to inspect or debug a single
#' stage in isolation.
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
                                     cluster_palette  = "npg",
                                     default_genes    = NULL,
                                     emit             = c("graph", "bundle"),
                                     verbose          = FALSE) {
  cells <- extract_cells(
    expression_eset,
    cell_id_col    = cell_id_col,
    cell_type_col  = cell_type_col,
    cell_group_col = cell_group_col,
    coordinate_col = coordinate_col
  )
  genes <- extract_genes(expression_eset, gene_symbol_col = gene_symbol_col)
  expr  <- extract_expression(expression_eset, genes = genes)
  act   <- if (!is.null(activity_eset)) {
    extract_activity(activity_eset, master_genes = genes)
  } else list(tf = NULL, sig = NULL)
  nets  <- if (!is.null(networks_path)) {
    read_networks(networks_path)
  } else list(tf = NULL, sig = NULL)
  meta$coordinate <- meta$coordinate %||% coordinate_col
  default_genes <- validate_default_genes(default_genes, genes)

  prepare_study_data(
    out_dir         = out_dir,
    meta            = meta,
    cells           = cells,
    clusters        = clusters,
    cluster_palette = cluster_palette,
    genes           = genes,
    expression      = expr,
    activity_tf     = act$tf,
    activity_sig    = act$sig,
    network_tf      = nets$tf,
    network_sig     = nets$sig,
    default_genes   = default_genes,
    emit            = emit,
    verbose         = verbose
  )
}

#' @rdname prepare_study
#' @param config_path Path to a YAML config.
#'
#' @details `prepare_study()` is a thin wrapper that calls
#' [load_study_config()] to parse + validate the YAML, then `readRDS()`
#' on the input files, then [prepare_study_from_eset()]. Call those
#' two pieces directly when you want to inspect the parsed config or
#' the loaded ExpressionSets before writing anything to disk.
#' @export
prepare_study <- function(config_path,
                           emit = c("graph", "bundle"),
                           verbose = FALSE) {
  cfg <- load_study_config(config_path)

  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase is required for prepare_study(); ",
         "install it from Bioconductor.")
  }
  expression_eset <- readRDS(cfg$input$expression)
  activity_eset <- if (!is.null(cfg$input$activity)) {
    readRDS(cfg$input$activity)
  } else NULL

  meta <- list(
    studyID    = as.character(cfg$study$ID),
    studyAbbr  = as.character(cfg$study$studyAbbr),
    longTitle  = as.character(cfg$study$longTitle),
    shortTitle = as.character(cfg$study$shortTitle),
    species    = cfg$species,
    coordinate = cfg$coordinate
  )

  prepare_study_from_eset(
    out_dir          = cfg$output,
    expression_eset  = expression_eset,
    activity_eset    = activity_eset,
    networks_path    = cfg$input$networks,
    meta             = meta,
    cell_id_col      = cfg$cellID,
    cell_type_col    = cfg$cellType,
    cell_group_col   = cfg$cellGroup,
    coordinate_col   = cfg$coordinate,
    gene_symbol_col  = cfg$geneSymbol,
    cluster_palette  = cfg$cluster_palette,
    default_genes    = cfg$default_genes,
    emit             = emit,
    verbose          = verbose
  )
}

# ============================================================================
# Staged helpers — call these directly for granular control / debug.
# Each helper does one step of the pipeline and returns an inspectable R
# structure (no side effects).
# ============================================================================

#' Parse + validate a `prepare_study` YAML config.
#'
#' Reads the YAML file at `config_path`, validates required keys
#' (`output`, `study.{ID,studyAbbr,longTitle,shortTitle}`,
#' `input.expression`), and fills in all optional keys with their
#' defaults. Pure parsing — does not touch the RDS / TSV files
#' referenced by `input.*`.
#'
#' @param config_path Path to a YAML config.
#' @return A named list with the parsed config + defaults applied.
#' @seealso [prepare_study()] for the full orchestration.
#' @export
load_study_config <- function(config_path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("yaml is required for load_study_config(); install.packages('yaml').")
  }
  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path)
  }
  cfg <- yaml::yaml.load_file(config_path)

  for (k in c("output", "study", "input")) {
    if (is.null(cfg[[k]])) {
      stop(sprintf("Config %s is missing required top-level key: '%s'",
                   config_path, k))
    }
  }
  for (k in c("ID", "studyAbbr", "longTitle", "shortTitle")) {
    if (is.null(cfg$study[[k]])) {
      stop(sprintf("Config %s is missing required key: study.%s",
                   config_path, k))
    }
  }
  if (is.null(cfg$input$expression)) {
    stop(sprintf("Config %s is missing required key: input.expression",
                 config_path))
  }

  # Fill defaults (preserving the existing values where set)
  cfg$species         <- as.character(cfg$species         %||% "")
  cfg$coordinate      <- as.character(cfg$coordinate      %||% "UMAP")
  cfg$cellID          <- as.character(cfg$cellID          %||% "cellID")
  cfg$cellType        <- as.character(cfg$cellType        %||% "cellGroup")
  cfg$cellGroup       <- as.character(cfg$cellGroup       %||% cfg$cellType)
  cfg$geneSymbol      <- as.character(cfg$geneSymbol      %||% "geneSymbol")
  cfg$cluster_palette <- as.character(cfg$cluster_palette %||% "npg")
  cfg$output          <- as.character(cfg$output)

  # Default genes: the app loads these on launch (mirrors the scMINER
  # Portal's `preGenes` field). Accepts any of three YAML shapes:
  #   default_genes: [GeneA, GeneB, GeneC]      # list
  #   default_genes: "GeneA, GeneB, GeneC"      # comma- or whitespace-
  #                                             # separated string
  #   defaults: { genes: [...] }                # nested form
  # The legacy `preGenes:` key is accepted as an alias.
  cfg$default_genes <- parse_default_genes(
    cfg$default_genes      %||%
    cfg$defaults$genes     %||%
    cfg$preGenes           %||%
    NULL
  )
  cfg
}

#' Normalise a YAML `default_genes:` value into a character vector.
#'
#' Accepts either:
#'   * an R list / character vector (returned with whitespace trimmed
#'     and empty entries dropped),
#'   * a single string of comma-, semicolon-, whitespace-, or newline-
#'     separated gene symbols, or
#'   * `NULL` (returns `NULL`).
#'
#' Duplicates are dropped (first occurrence wins). Used by
#' [load_study_config()]; exported so callers can apply the same
#' parsing rules to ad-hoc input.
#'
#' @param x A YAML-derived value (list, character vector, or string).
#' @return `NULL` (if `x` is empty) or a deduplicated character vector
#'   of gene symbols.
#' @export
parse_default_genes <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  if (!is.character(x)) x <- as.character(x)
  # Split anything that still has internal separators.
  parts <- unlist(strsplit(x, "[,;[:space:]\n]+"), use.names = FALSE)
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) return(NULL)
  unique(parts)
}

#' Filter a default-genes list to those actually present in a master
#' gene set, warning about any drops.
#'
#' Both [prepare_study()] and [prepare_study_from_eset()] call this
#' before passing `default_genes` to [prepare_study_data()] / the
#' bundle writer so the on-disk `defaults/genes` only references
#' genes the app can actually load.
#'
#' @param default_genes Character vector or `NULL`.
#' @param master_genes The bundle's master gene list.
#' @param warn If `TRUE`, emit a `warning()` listing missing genes.
#' @return `NULL` if `default_genes` is empty, otherwise the subset
#'   of `default_genes` present in `master_genes` (in user-specified
#'   order).
#' @export
validate_default_genes <- function(default_genes, master_genes,
                                     warn = TRUE) {
  if (is.null(default_genes) || length(default_genes) == 0L) return(NULL)
  default_genes <- as.character(default_genes)
  master_genes  <- as.character(master_genes)
  hit  <- default_genes %in% master_genes
  miss <- default_genes[!hit]
  if (length(miss) > 0L && isTRUE(warn)) {
    warning(sprintf(
      "default_genes: %d/%d not in master gene list (dropping): %s",
      length(miss), length(default_genes),
      paste(utils::head(miss, 6L), collapse = ", ")),
      call. = FALSE)
  }
  kept <- default_genes[hit]
  if (length(kept) == 0L) NULL else kept
}

#' Extract the cells data.frame from an `ExpressionSet`.
#'
#' Pulls `pData()`, coerces factor columns to character, and selects /
#' renames the cellID, cellType, cellGroup, and coordinate columns into
#' the canonical `cells` data.frame consumed by [prepare_study_data()].
#'
#' Errors loudly if any required pData column is missing — useful when
#' debugging mis-named source columns.
#'
#' @param expression_eset A Biobase `ExpressionSet`.
#' @param cell_id_col,cell_type_col,cell_group_col Column names in
#'   `pData(expression_eset)`. `cell_id_col` defaults to the eset's
#'   rownames if not already a column.
#' @param coordinate_col Stem of the coord columns; the eset must have
#'   `<coordinate_col>_1` and `<coordinate_col>_2`.
#' @return data.frame with columns `cellID`, `cellType`, `cellGroup`,
#'   `coord1`, `coord2`.
#' @export
extract_cells <- function(expression_eset,
                           cell_id_col    = "cellID",
                           cell_type_col  = "cellGroup",
                           cell_group_col = cell_type_col,
                           coordinate_col = "UMAP") {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase is required for extract_cells(); ",
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
  group_vals <- p_data[[cell_group_col]] %||% p_data[[cell_type_col]]
  data.frame(
    cellID    = as.character(p_data[[cell_id_col]]),
    cellType  = as.character(p_data[[cell_type_col]]),
    cellGroup = as.character(group_vals),
    coord1    = as.numeric(p_data[[coord1_col]]),
    coord2    = as.numeric(p_data[[coord2_col]]),
    stringsAsFactors = FALSE
  )
}

#' Extract the master gene symbol vector from an `ExpressionSet`.
#'
#' @param expression_eset A Biobase `ExpressionSet`.
#' @param gene_symbol_col Column name in `fData(expression_eset)`.
#' @return Character vector of gene symbols (length = number of rows).
#' @export
extract_genes <- function(expression_eset, gene_symbol_col = "geneSymbol") {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase is required for extract_genes(); ",
         "install it from Bioconductor.")
  }
  f_data <- as.data.frame(Biobase::fData(expression_eset),
                          stringsAsFactors = FALSE)
  if (!gene_symbol_col %in% colnames(f_data)) {
    stop("fData missing gene-symbol column: ", gene_symbol_col)
  }
  as.character(f_data[[gene_symbol_col]])
}

#' Extract the expression matrix from an `ExpressionSet`.
#'
#' Returns a CSC sparse `Matrix` (genes × cells). If `genes` is given,
#' the matrix's rownames are set from it — useful when the eset's
#' fData has more useful symbols than the matrix's existing rownames.
#'
#' @param expression_eset A Biobase `ExpressionSet`.
#' @param genes Optional character vector; replaces the matrix's
#'   rownames if `nrow(exprs(eset)) == length(genes)`.
#' @return A `dgCMatrix` of shape `(nrow(eset), ncol(eset))`.
#' @export
extract_expression <- function(expression_eset, genes = NULL) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase is required for extract_expression(); ",
         "install it from Bioconductor.")
  }
  expr <- Biobase::exprs(expression_eset)
  # NOTE: we deliberately do NOT convert to CsparseMatrix here.
  # .write_graph_shards iterates row-by-row and handles dgeMatrix /
  # dgCMatrix / base matrix uniformly via mat[i, , drop = TRUE], so
  # forcing a dense->sparse cast would (a) double-allocate, and (b)
  # blow up on studies whose nnz approaches the 2^31-1 dgCMatrix cap
  # (e.g. > 100k cells x 20k genes scRNA-seq). Keep whatever class
  # exprs(eset) returned.
  if (!is.null(genes)) {
    if (nrow(expr) != length(genes)) {
      stop(sprintf("Expression rows (%d) != length(genes) (%d)",
                   nrow(expr), length(genes)))
    }
    rownames(expr) <- genes
  }
  expr
}

#' Extract activity matrices (TF + SIG) from an `ExpressionSet`.
#'
#' Splits rows by suffix (`_TF` / `.TF` → TF; `_SIG` / `.SIG` → SIG),
#' strips the suffix from row names, and reindexes both sub-matrices
#' onto the master `genes` list (rows not present in the activity
#' source become zero rows).
#'
#' @param activity_eset A Biobase `ExpressionSet` whose row names carry
#'   the `_TF` / `_SIG` suffixes.
#' @param master_genes Character vector — the bundle's master gene
#'   list. Both output matrices have this many rows in this order.
#' @return A list with elements `tf` and `sig`, each either a sparse
#'   `Matrix` (rows in `master_genes` order) or `NULL` if no rows of
#'   that kind were present.
#' @export
extract_activity <- function(activity_eset, master_genes) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase is required for extract_activity(); ",
         "install it from Bioconductor.")
  }
  act <- Biobase::exprs(activity_eset)
  act_rows <- rownames(act)
  is_tf  <- grepl("[._]TF$",  act_rows)
  is_sig <- grepl("[._]SIG$", act_rows)
  if (!any(is_tf) && !any(is_sig)) {
    warning("activity_eset has no rows ending in _TF/.TF or _SIG/.SIG; ",
            "returning NULL for both kinds.", call. = FALSE)
    return(list(tf = NULL, sig = NULL))
  }
  strip_suffix <- function(rows) sub("[._](TF|SIG)$", "", rows)

  build_one <- function(mask) {
    if (!any(mask)) return(NULL)
    sub_act <- act[mask, , drop = FALSE]
    rownames(sub_act) <- strip_suffix(rownames(sub_act))
    # Same rationale as extract_expression: keep the source class
    # rather than forcing a dense->sparse cast that allocates twice
    # and trips the 2^31-1 dgCMatrix cap on large activity matrices.
    # .reindex_rows builds its own sparse output of the master shape.
    .reindex_rows(sub_act, master_genes)
  }
  list(tf = build_one(is_tf), sig = build_one(is_sig))
}

#' Read a scMINER networks TSV.
#'
#' Parses a tab-separated networks file with columns
#' `source, target, NetworkType, [studyID,] CellGroup, mi, pearson,
#' spearman, rho, pvalue` (the format consumed by the original
#' `h_networks.R`) and splits it into the canonical TF / SIG data.frames
#' [prepare_study_data()] expects.
#'
#' @param path Path to the networks TSV.
#' @return A list with elements `tf` and `sig`, each a data.frame with
#'   columns `source, target, cellType, mi, pearson, spearman, rho,
#'   pvalue`, or `NULL` if that NetworkType wasn't present in the file.
#' @export
read_networks <- function(path) {
  .read_networks_file(path)
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
