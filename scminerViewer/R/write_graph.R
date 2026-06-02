# Internal helpers that write the graph-import layout produced by
# scMINER-portal-datapre-R. Ported from h_metadata.R, h_activity.R,
# h_networks.R; per-row write loops vectorised; hard-coded `.libPaths`
# and YAML I/O removed (those live in prepare_study.R).
#
# All helpers take plain R structures (data.frame, character vector,
# Matrix) — no Biobase dependency at this layer.

.graph_subdirs <- c("Header", "Study", "Gene", "Cell",
                    "Study_Contains_Gene", "Study_Contains_Cell",
                    "Network_TF_Activity", "Network_SIG_Activity",
                    "study_gene_expression", "study_gene_tf",
                    "study_gene_sig", "study_meta")

.shard_letters <- c(letters, "nm")

.ensure_graph_tree <- function(out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (sub in .graph_subdirs) {
    dir.create(file.path(out_dir, sub),
               showWarnings = FALSE, recursive = TRUE)
  }
}

.ensure_shard_tree <- function(out_dir, kind) {
  # kind: "expression_files" | "activity_files/TF" | "activity_files/SIG"
  base <- file.path(out_dir, kind)
  dir.create(base, showWarnings = FALSE, recursive = TRUE)
  for (letter in .shard_letters) {
    dir.create(file.path(base, letter),
               showWarnings = FALSE, recursive = TRUE)
  }
  base
}

# --- Study / clusters --------------------------------------------------------

.write_graph_study <- function(out_dir, meta) {
  study_id <- as.character(meta$studyID)
  head_title <- c("studyID:ID(Study)", "studyAbbr",
                  "longTitle", "shortTitle")
  writeLines(paste(head_title, collapse = "\t"),
             file.path(out_dir, "Header",
                       paste0(study_id, "_study.header.tsv")))
  writeLines(paste(c(meta$studyID, meta$studyAbbr,
                     meta$longTitle, meta$shortTitle),
                   collapse = "\t"),
             file.path(out_dir, "Study",
                       paste0(study_id, "_study.tsv")))
}

.write_graph_clusters <- function(out_dir, meta, clusters) {
  if (is.null(clusters) || nrow(clusters) == 0) return(invisible(NULL))
  study_id  <- as.character(meta$studyID)
  study_abbr <- as.character(meta$studyAbbr)
  if (!"NetworkCellType" %in% colnames(clusters)) {
    clusters$NetworkCellType <- clusters$cellType
  }
  if (!"cellGroup" %in% colnames(clusters)) {
    clusters$cellGroup <- clusters$cellType
  }
  if (!"color" %in% colnames(clusters)) {
    clusters$color <- "#888888"
  }
  if (!"label_1" %in% colnames(clusters)) clusters$label_1 <- 0
  if (!"label_2" %in% colnames(clusters)) clusters$label_2 <- 0
  out <- data.frame(
    StudyID         = study_id,
    StudyAbbr       = study_abbr,
    CellType        = as.character(clusters$cellType),
    CellGroup       = as.character(clusters$cellGroup),
    Color           = as.character(clusters$color),
    Label_1         = as.numeric(clusters$label_1),
    Label_2         = as.numeric(clusters$label_2),
    NetworkCellType = as.character(clusters$NetworkCellType),
    stringsAsFactors = FALSE
  )
  # Per-study cluster meta file. The older shared layout used a single
  # `study_meta.csv` that all studies appended to; the new layout writes
  # one file per study, so we always emit the header.
  path <- file.path(out_dir, "study_meta",
                    paste0(study_id, "_study_meta.csv"))
  utils::write.table(out, file = path,
                     sep = ",", quote = FALSE, row.names = FALSE,
                     col.names = TRUE, append = FALSE)
}

# --- Genes -------------------------------------------------------------------

.write_graph_genes <- function(out_dir, meta, genes) {
  study_id <- as.character(meta$studyID)
  writeLines("geneSymbol:ID(Gene)",
             file.path(out_dir, "Header",
                       paste0(study_id, "_n_gene.header.tsv")))
  writeLines(as.character(genes),
             file.path(out_dir, "Gene",
                       paste0(study_id, "_n_gene.tsv")))
  writeLines(":START_ID(Study)\t:END_ID(Gene)",
             file.path(out_dir, "Header",
                       paste0(study_id, "_r_study_gene.header.tsv")))
  writeLines(paste(study_id, as.character(genes), sep = "\t"),
             file.path(out_dir, "Study_Contains_Gene",
                       paste0(study_id, "_r_study_gene.tsv")))
}

# --- Cells -------------------------------------------------------------------

.write_graph_cells <- function(out_dir, meta, cells) {
  study_id   <- as.character(meta$studyID)
  coordinate <- as.character(meta$coordinate %||% "UMAP")

  head_title <- c(
    "cellID:ID(Cell)", "cellType", "cellGroup",
    paste0(coordinate, "_1:float"),
    paste0(coordinate, "_2:float"),
    "coordinateName"
  )
  writeLines(paste(head_title, collapse = "\t"),
             file.path(out_dir, "Header",
                       paste0(study_id, "_n_cell.header.tsv")))

  cell_id    <- as.character(cells$cellID)
  cell_type  <- as.character(cells$cellType)
  cell_group <- as.character(cells$cellGroup %||% cells$cellType)
  coord1     <- as.character(cells$coord1)
  coord2     <- as.character(cells$coord2)

  # Preserve the original 7-column form: cellID is written twice (the
  # original sapply over a data.frame row produced this layout, and the
  # downstream Java backend expects it). Vectorised here.
  lines <- paste(cell_id, cell_id, cell_type, cell_group,
                 coord1, coord2, coordinate, sep = "\t")
  writeLines(lines,
             file.path(out_dir, "Cell",
                       paste0(study_id, "_n_cell.tsv")))

  writeLines(":START_ID(Study)\t:END_ID(Cell)",
             file.path(out_dir, "Header",
                       paste0(study_id, "_r_study_cell.header.tsv")))
  writeLines(paste(study_id, cell_id, sep = "\t"),
             file.path(out_dir, "Study_Contains_Cell",
                       paste0(study_id, "_r_study_cell.tsv")))
}

# --- Networks ----------------------------------------------------------------

.write_graph_networks <- function(out_dir, meta, network_tf, network_sig) {
  study_id <- as.character(meta$studyID)
  for (kind in c("TF", "SIG")) {
    df <- if (kind == "TF") network_tf else network_sig
    if (is.null(df) || nrow(df) == 0) next
    sub <- if (kind == "TF") "Network_TF_Activity" else "Network_SIG_Activity"
    writeLines(
      paste(":START_ID(Gene)", ":END_ID(Gene)",
            "relationshipType:string", "studyID:string", "cellType:string",
            "mi:float", "pearson:float", "spearman:float",
            "rho:float", "pvalue:float", sep = "\t"),
      file.path(out_dir, "Header",
                paste0(study_id, "_", kind, ".header.tsv"))
    )
    lines <- paste(
      as.character(df$source), as.character(df$target),
      kind, study_id, as.character(df$cellType),
      df$mi, df$pearson, df$spearman, df$rho, df$pvalue,
      sep = "\t"
    )
    writeLines(lines,
               file.path(out_dir, sub,
                         paste0(study_id, "_", kind, ".tsv")))
  }
}

# --- Per-gene shards (expression + activity) ---------------------------------

.write_graph_shards <- function(out_dir, meta, mat, kind, meta_kind,
                                manifest_dir, manifest_name, cell_ids,
                                type_label, verbose = FALSE,
                                progress = FALSE) {
  # mat: Matrix (G x N), rows = genes, columns = cells
  # kind:      where per-letter shard dirs live, relative to out_dir.
  #            e.g. "expression_files/<studyID>" or "activity_files/<studyID>/TF"
  # meta_kind: where the cell-header meta.csv lives, relative to out_dir.
  #            For expression this is the same as `kind`. For activity the
  #            meta is one level above (TF and SIG share one meta.csv).
  # manifest_dir: e.g. "study_gene_expression"
  # manifest_name: e.g. "expression"
  # type_label: written into manifest's Type column
  if (is.null(mat) || nrow(mat) == 0) return(invisible(NULL))

  study_id <- as.character(meta$studyID)
  abbr     <- as.character(meta$studyAbbr)
  species  <- as.character(meta$species %||% "")

  shard_base <- .ensure_shard_tree(out_dir, kind)

  # Cell-header may live above the shard root (activity TF/SIG share one).
  meta_dir <- file.path(out_dir, meta_kind)
  dir.create(meta_dir, showWarnings = FALSE, recursive = TRUE)
  meta_csv_rel <- file.path(meta_kind, "meta.csv")
  writeLines(paste(as.character(cell_ids), collapse = ","),
             file.path(meta_dir, "meta.csv"))

  manifest_path <- file.path(out_dir, manifest_dir,
                             paste0(study_id, "_", manifest_name, ".csv"))
  manifest_header <- c("GeneSymbol", "Species", "StudyID", "StudyAbbr",
                       "Type", "FileHeader", "File")
  con <- file(manifest_path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste(manifest_header, collapse = ","), con)

  genes <- rownames(mat)
  if (is.null(genes)) {
    stop("Matrix passed to .write_graph_shards must have rownames (gene symbols)")
  }
  ncells <- ncol(mat)

  # A live progress bar over the per-gene shard loop — this is the slow
  # part of prepare_study (one gzipped file per gene). Falls back to the
  # every-1000 verbose message when `progress` is off.
  pb <- NULL
  if (isTRUE(progress) && nrow(mat) > 0) {
    pb <- utils::txtProgressBar(min = 0, max = nrow(mat), style = 3)
    on.exit(close(pb), add = TRUE)
  }

  is_sparse <- inherits(mat, "Matrix")
  for (i in seq_len(nrow(mat))) {
    gene <- genes[i]
    name <- gsub("/", "_", gene, fixed = TRUE)
    letter <- .shard_letter(gene)
    rel <- file.path(kind, letter, paste0(name, ".csv.gz"))
    writeLines(
      paste(c(gene, species, study_id, abbr, type_label,
              meta_csv_rel, rel), collapse = ","),
      con
    )
    if (is_sparse) {
      row_vals <- as.numeric(mat[i, , drop = TRUE])
    } else {
      row_vals <- as.numeric(mat[i, ])
    }
    csv_path <- file.path(shard_base, letter, paste0(name, ".csv"))
    # Use fwrite to preserve full double precision (the original
    # h_activity.R also relied on fwrite for this). Pre-coerce to
    # data.frame to suppress the "coerced from matrix" notice.
    row_df <- as.data.frame(matrix(row_vals, nrow = 1L))
    data.table::fwrite(
      row_df,
      file      = csv_path,
      row.names = FALSE,
      col.names = FALSE,
      quote     = FALSE,
      sep       = ","
    )
    R.utils::gzip(csv_path, destname = paste0(csv_path, ".gz"),
                  overwrite = TRUE, remove = TRUE)
    if (!is.null(pb)) {
      utils::setTxtProgressBar(pb, i)
    } else if (isTRUE(verbose) && (i %% 1000L == 0L)) {
      message(sprintf("  [%s] %d/%d shards written",
                      manifest_name, i, nrow(mat)))
    }
  }
}
