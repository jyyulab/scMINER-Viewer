#' Read a study from the existing graph-import directory layout.
#'
#' Reconstructs the inputs that [write_bundle()] / [prepare_study_data()]
#' consume from the on-disk layout under `data/`. Matrix values are
#' *not* loaded (those are read on demand by [gene_values()] from the
#' sharded `expression_files/<studyID>/` and `activity_files/<studyID>/`
#' trees). Instead, per-matrix **gene indexes** are read from the three
#' manifest CSVs and returned for the bundle to embed.
#'
#' @param data_dir Root directory containing `Study/`, `Cell/`, `Gene/`,
#'   `Network_TF_Activity/`, `Network_SIG_Activity/`, `study_meta/`,
#'   `study_gene_expression/`, `study_gene_tf/`, `study_gene_sig/`.
#' @param study_id Study identifier (e.g., `"2327"`).
#'
#' @return A list with `meta`, `cells`, `clusters`, `genes`,
#'   `expression_genes`, `activity_tf_genes`, `activity_sig_genes`
#'   (character vectors from the matching manifest), and `network_tf`,
#'   `network_sig` data.frames. Any slot is `NULL` when its source file
#'   isn't present.
#' @export
read_graph_study <- function(data_dir, study_id) {
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

  # Per-matrix gene indexes are read from the manifest CSVs (we ignore
  # the manifest's File / FileHeader columns since the reader derives
  # those from gene + kind).
  expression_genes <- .read_manifest_genes(
    file.path(data_dir, "study_gene_expression",
              paste0(study_id, "_expression.csv"))
  )
  activity_tf_genes <- .read_manifest_genes(
    file.path(data_dir, "study_gene_tf",
              paste0(study_id, "_activity_tf.csv"))
  )
  activity_sig_genes <- .read_manifest_genes(
    file.path(data_dir, "study_gene_sig",
              paste0(study_id, "_activity_sig.csv"))
  )

  list(
    meta               = meta,
    cells              = cells,
    clusters           = clusters,
    genes              = genes,
    expression_genes   = expression_genes,
    activity_tf_genes  = activity_tf_genes,
    activity_sig_genes = activity_sig_genes,
    network_tf         = network_tf,
    network_sig        = network_sig
  )
}

.read_manifest_genes <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || !"GeneSymbol" %in% colnames(df) || nrow(df) == 0) {
    return(NULL)
  }
  unique(as.character(df$GeneSymbol))
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

