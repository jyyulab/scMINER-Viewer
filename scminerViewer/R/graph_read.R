#' Read a study from the existing graph-import directory layout.
#'
#' Reconstructs the inputs that [write_bundle()] consumes from the on-disk
#' layout under `data/`.
#'
#' Expression and activity matrices are sharded by the first lowercase
#' letter of each gene (or `nm` for non-alphabetic first chars), with a
#' single cell-ID header per matrix type:
#'
#' \itemize{
#'   \item `<shard_dir>/expression_files/meta.csv` (cell-ID column order
#'     for every expression shard)
#'   \item `<shard_dir>/expression_files/<letter>/<gene>.csv.gz`
#'   \item `<shard_dir>/activity_files/TF/meta.csv` and
#'     `<shard_dir>/activity_files/TF/<letter>/<gene>.csv.gz`
#'   \item `<shard_dir>/activity_files/SIG/meta.csv` and
#'     `<shard_dir>/activity_files/SIG/<letter>/<gene>.csv.gz`
#' }
#'
#' Set the `load_expression`, `load_activity_tf`, `load_activity_sig`
#' flags to include those matrices. The manifest CSVs under
#' `study_gene_expression/`, `study_gene_tf/`, `study_gene_sig/` are used
#' only to enumerate which genes have shards; their `File` / `FileHeader`
#' columns are ignored (paths are derived from the gene name). Missing
#' shards or per-letter meta files are tolerated and reported as warnings.
#'
#' @param data_dir Root directory containing `Study/`, `Cell/`, `Gene/`,
#'   `Network_TF_Activity/`, `Network_SIG_Activity/`, `study_meta/`,
#'   `study_gene_expression/`, `study_gene_tf/`, `study_gene_sig/`.
#' @param study_id Study identifier (e.g., `"2327"`).
#' @param shard_dir Root directory under which `expression_files/` and
#'   `activity_files/` live. Defaults to `data_dir`.
#' @param load_expression,load_activity_tf,load_activity_sig Logical;
#'   load the corresponding matrix from per-gene shards.
#' @param verbose Logical; emit progress messages while reading shards.
#'
#' @return A list with `meta`, `cells`, `clusters`, `genes`, `expression`,
#'   `activity_tf`, `activity_sig`, `network_tf`, `network_sig`. Matrix
#'   slots are `NULL` when the corresponding `load_*` flag is `FALSE` or
#'   the manifest is not found.
#' @export
read_graph_study <- function(data_dir, study_id,
                             shard_dir = data_dir,
                             load_expression = FALSE,
                             load_activity_tf = FALSE,
                             load_activity_sig = FALSE,
                             verbose = FALSE) {
  study_id <- as.character(study_id)

  study_file <- file.path(data_dir, "Study", paste0(study_id, "_study.tsv"))
  if (!file.exists(study_file)) {
    stop("Study file not found: ", study_file)
  }
  study_row <- strsplit(readLines(study_file, n = 1, warn = FALSE),
                        "\t", fixed = TRUE)[[1]]
  if (length(study_row) < 4) {
    stop("Study file ", study_file, " has fewer than 4 fields")
  }

  # Study meta lives in a per-study file: study_meta/<studyID>_study_meta.csv.
  # The older shared-file form (study_meta/study_meta.csv) is still accepted
  # as a fallback so old exports keep working.
  meta_csv_path <- file.path(data_dir, "study_meta",
                             paste0(study_id, "_study_meta.csv"))
  if (!file.exists(meta_csv_path)) {
    fallback <- file.path(data_dir, "study_meta", "study_meta.csv")
    if (file.exists(fallback)) meta_csv_path <- fallback
  }
  meta_csv <- if (file.exists(meta_csv_path)) {
    read.csv(meta_csv_path, stringsAsFactors = FALSE)
  } else {
    data.frame()
  }
  if (nrow(meta_csv) > 0 && !is.null(meta_csv$StudyID)) {
    meta_csv <- meta_csv[as.character(meta_csv$StudyID) == study_id, ,
                         drop = FALSE]
  }

  cell_file <- file.path(data_dir, "Cell", paste0(study_id, "_n_cell.tsv"))
  if (!file.exists(cell_file)) {
    stop("Cell file not found: ", cell_file)
  }
  cell_df <- data.table::fread(cell_file, header = FALSE, sep = "\t",
                               data.table = FALSE)
  # The graph export writes the cellID column twice (artifact of the
  # original sapply over a data.frame row); 7-column form is what 2327
  # produces. We tolerate the older 6-column form too.
  if (ncol(cell_df) >= 7) {
    cells <- data.frame(
      cellID    = as.character(cell_df[[2]]),
      cellType  = as.character(cell_df[[3]]),
      cellGroup = as.character(cell_df[[4]]),
      coord1    = as.numeric(cell_df[[5]]),
      coord2    = as.numeric(cell_df[[6]]),
      stringsAsFactors = FALSE
    )
    coordinate_name <- unique(as.character(cell_df[[7]]))[1]
  } else if (ncol(cell_df) == 6) {
    cells <- data.frame(
      cellID    = as.character(cell_df[[1]]),
      cellType  = as.character(cell_df[[2]]),
      cellGroup = as.character(cell_df[[3]]),
      coord1    = as.numeric(cell_df[[4]]),
      coord2    = as.numeric(cell_df[[5]]),
      stringsAsFactors = FALSE
    )
    coordinate_name <- unique(as.character(cell_df[[6]]))[1]
  } else {
    stop("Unexpected Cell tsv column count: ", ncol(cell_df))
  }

  species_from_meta <- if (nrow(meta_csv) > 0 && !is.null(meta_csv$Species)) {
    as.character(meta_csv$Species[1])
  } else NA_character_

  meta <- list(
    studyID    = study_row[1],
    studyAbbr  = study_row[2],
    longTitle  = study_row[3],
    shortTitle = study_row[4],
    species    = if (!is.na(species_from_meta)) species_from_meta else "",
    coordinate = if (!is.na(coordinate_name)) coordinate_name else "UMAP"
  )

  cnt <- as.data.frame(table(cells$cellType), stringsAsFactors = FALSE)
  if (nrow(meta_csv) > 0) {
    clusters <- data.frame(
      cellType = meta_csv$CellType,
      count    = cnt$Freq[match(meta_csv$CellType, cnt$Var1)],
      color    = meta_csv$Color,
      stringsAsFactors = FALSE
    )
    if (!is.null(meta_csv$Label_1)) {
      clusters$label_1 <- as.numeric(meta_csv$Label_1)
    }
    if (!is.null(meta_csv$Label_2)) {
      clusters$label_2 <- as.numeric(meta_csv$Label_2)
    }
  } else {
    clusters <- data.frame(
      cellType = cnt$Var1,
      count    = cnt$Freq,
      color    = rep("#888888", nrow(cnt)),
      stringsAsFactors = FALSE
    )
  }
  clusters$count[is.na(clusters$count)] <- 0L

  gene_file <- file.path(data_dir, "Gene", paste0(study_id, "_n_gene.tsv"))
  if (!file.exists(gene_file)) {
    stop("Gene file not found: ", gene_file)
  }
  genes <- readLines(gene_file, warn = FALSE)

  network_tf  <- .read_graph_network(
    file.path(data_dir, "Network_TF_Activity",
              paste0(study_id, "_TF.tsv"))
  )
  network_sig <- .read_graph_network(
    file.path(data_dir, "Network_SIG_Activity",
              paste0(study_id, "_SIG.tsv"))
  )

  # Try to enrich species from one of the matrix manifests when missing
  if (!nzchar(meta$species)) {
    for (m in c(file.path(data_dir, "study_gene_expression",
                          paste0(study_id, "_expression.csv")),
                file.path(data_dir, "study_gene_tf",
                          paste0(study_id, "_activity_tf.csv")),
                file.path(data_dir, "study_gene_sig",
                          paste0(study_id, "_activity_sig.csv")))) {
      if (file.exists(m)) {
        first <- tryCatch(
          utils::read.csv(m, stringsAsFactors = FALSE, nrows = 1),
          error = function(e) NULL
        )
        if (!is.null(first) && !is.null(first$Species) && nrow(first) > 0) {
          meta$species <- as.character(first$Species[1])
          break
        }
      }
    }
  }

  # Shard tree: expression_files/<studyID>/{meta.csv, <letter>/<gene>.csv.gz}
  # and activity_files/<studyID>/{meta.csv, TF/<letter>/..., SIG/<letter>/...}.
  # Activity TF and SIG share the same meta.csv (one level above the kind dir).
  exp_root  <- file.path("expression_files", study_id)
  act_root  <- file.path("activity_files",   study_id)
  act_meta  <- file.path(act_root, "meta.csv")
  exp_meta  <- file.path(exp_root, "meta.csv")

  expression <- if (isTRUE(load_expression)) {
    .read_shard_matrix(
      shard_dir     = shard_dir,
      shard_root    = exp_root,
      meta_rel_path = exp_meta,
      manifest_path = file.path(data_dir, "study_gene_expression",
                                paste0(study_id, "_expression.csv")),
      genes         = genes,
      cells         = cells$cellID,
      label         = "expression",
      verbose       = verbose
    )
  } else NULL

  activity_tf <- if (isTRUE(load_activity_tf)) {
    .read_shard_matrix(
      shard_dir     = shard_dir,
      shard_root    = file.path(act_root, "TF"),
      meta_rel_path = act_meta,
      manifest_path = file.path(data_dir, "study_gene_tf",
                                paste0(study_id, "_activity_tf.csv")),
      genes         = genes,
      cells         = cells$cellID,
      label         = "activity_tf",
      verbose       = verbose
    )
  } else NULL

  activity_sig <- if (isTRUE(load_activity_sig)) {
    .read_shard_matrix(
      shard_dir     = shard_dir,
      shard_root    = file.path(act_root, "SIG"),
      meta_rel_path = act_meta,
      manifest_path = file.path(data_dir, "study_gene_sig",
                                paste0(study_id, "_activity_sig.csv")),
      genes         = genes,
      cells         = cells$cellID,
      label         = "activity_sig",
      verbose       = verbose
    )
  } else NULL

  list(
    meta         = meta,
    cells        = cells,
    clusters     = clusters,
    genes        = genes,
    expression   = expression,
    activity_tf  = activity_tf,
    activity_sig = activity_sig,
    network_tf   = network_tf,
    network_sig  = network_sig
  )
}

.read_graph_network <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- data.table::fread(path, header = FALSE, sep = "\t",
                          data.table = FALSE)
  if (ncol(df) < 10) return(NULL)
  data.frame(
    source   = as.character(df[[1]]),
    target   = as.character(df[[2]]),
    cellType = as.character(df[[5]]),
    mi       = as.numeric(df[[6]]),
    pearson  = as.numeric(df[[7]]),
    spearman = as.numeric(df[[8]]),
    rho      = as.numeric(df[[9]]),
    pvalue   = as.numeric(df[[10]]),
    stringsAsFactors = FALSE
  )
}

.read_shard_matrix <- function(shard_dir, shard_root, meta_rel_path,
                               manifest_path, genes, cells, label, verbose) {
  if (!file.exists(manifest_path)) {
    warning(sprintf("[%s] manifest not found: %s -- skipping",
                    label, manifest_path), call. = FALSE)
    return(NULL)
  }
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  if (!"GeneSymbol" %in% colnames(manifest)) {
    warning(sprintf("[%s] manifest missing GeneSymbol column -- skipping",
                    label), call. = FALSE)
    return(NULL)
  }
  if (nrow(manifest) == 0) {
    warning(sprintf("[%s] manifest is empty -- skipping", label),
            call. = FALSE)
    return(NULL)
  }

  shard_root_path <- file.path(shard_dir, shard_root)
  if (!dir.exists(shard_root_path)) {
    warning(sprintf("[%s] shard root not found: %s -- skipping",
                    label, shard_root_path), call. = FALSE)
    return(NULL)
  }

  # Header path is independent of shard_root: activity TF and SIG share a
  # single meta.csv one level above their kind directory.
  meta_path <- file.path(shard_dir, meta_rel_path)
  if (!file.exists(meta_path)) {
    warning(sprintf("[%s] cell-header meta.csv not found at %s -- skipping",
                    label, meta_path), call. = FALSE)
    return(NULL)
  }
  header_line <- tryCatch(readLines(meta_path, n = 1, warn = FALSE),
                          error = function(e) character(0))
  if (length(header_line) == 0 || !nzchar(header_line)) {
    warning(sprintf("[%s] %s is empty -- skipping", label, meta_path),
            call. = FALSE)
    return(NULL)
  }
  shard_cells <- trimws(strsplit(header_line, ",", fixed = TRUE)[[1]])
  n_shard_cells <- length(shard_cells)
  perm <- match(cells, shard_cells)
  missing_in_shard <- which(is.na(perm))
  if (length(missing_in_shard) > 0) {
    warning(sprintf(
      "[%s] %d cell ID(s) from Cell/ not present in %s (e.g. %s); those columns will be zero",
      label, length(missing_in_shard), meta_path,
      paste(utils::head(cells[missing_in_shard], 3), collapse = ", ")
    ), call. = FALSE)
  }

  gene_to_row <- stats::setNames(seq_along(genes), genes)
  n_genes <- length(genes)
  n_cells <- length(cells)

  triplets_i <- vector("list", nrow(manifest))
  triplets_j <- vector("list", nrow(manifest))
  triplets_x <- vector("list", nrow(manifest))
  n_loaded <- 0L
  n_missing_shard <- 0L

  for (k in seq_len(nrow(manifest))) {
    gene <- manifest$GeneSymbol[k]
    row_idx <- gene_to_row[[gene]]
    if (is.null(row_idx) || is.na(row_idx)) {
      n_missing_shard <- n_missing_shard + 1L
      next
    }
    letter <- .shard_letter(gene)
    shard_path <- file.path(shard_root_path, letter,
                            paste0(gene, ".csv.gz"))
    if (!file.exists(shard_path)) {
      n_missing_shard <- n_missing_shard + 1L
      next
    }
    vals_df <- tryCatch(
      data.table::fread(shard_path, header = FALSE, sep = ",",
                        data.table = FALSE),
      error = function(e) NULL
    )
    if (is.null(vals_df) || nrow(vals_df) == 0 || ncol(vals_df) == 0) {
      n_missing_shard <- n_missing_shard + 1L
      next
    }
    shard_vals <- as.numeric(vals_df[1, ])
    if (length(shard_vals) != n_shard_cells) {
      warning(sprintf(
        "[%s] gene %s: shard has %d values but meta.csv has %d -- skipping",
        label, gene, length(shard_vals), n_shard_cells),
        call. = FALSE)
      n_missing_shard <- n_missing_shard + 1L
      next
    }
    reordered <- shard_vals[perm]
    nz <- which(!is.na(reordered) & reordered != 0)
    if (length(nz) > 0) {
      triplets_i[[k]] <- rep.int(as.integer(row_idx), length(nz))
      triplets_j[[k]] <- as.integer(nz)
      triplets_x[[k]] <- reordered[nz]
    }
    n_loaded <- n_loaded + 1L
    if (isTRUE(verbose) && (k %% 1000L == 0L)) {
      message(sprintf("  [%s] %d/%d shards loaded",
                      label, k, nrow(manifest)))
    }
  }

  if (isTRUE(verbose)) {
    message(sprintf("[%s] %d shards loaded, %d skipped",
                    label, n_loaded, n_missing_shard))
  }

  all_i <- unlist(triplets_i, use.names = FALSE)
  all_j <- unlist(triplets_j, use.names = FALSE)
  all_x <- unlist(triplets_x, use.names = FALSE)
  if (is.null(all_i)) {
    all_i <- integer(0); all_j <- integer(0); all_x <- numeric(0)
  }
  Matrix::sparseMatrix(
    i = all_i, j = all_j, x = all_x,
    dims = c(n_genes, n_cells), repr = "C"
  )
}
