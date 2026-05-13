# Build a tiny synthetic graph layout + shard tree in tempdir, then verify
# read_graph_study reconstructs the expression / activity matrices in the
# right cell column order and tolerates missing shards.

make_fake_layout <- function(root, study_id = "1234",
                              cells, gene_expression, gene_activity_tf,
                              gene_activity_sig,
                              cell_order_in_shards = NULL,
                              drop_genes = character(0)) {
  if (is.null(cell_order_in_shards)) {
    cell_order_in_shards <- cells$cellID
  }
  for (sub in c("Study", "Cell", "Gene", "study_meta",
                "Network_TF_Activity", "Network_SIG_Activity",
                "study_gene_expression", "study_gene_tf",
                "study_gene_sig")) {
    dir.create(file.path(root, sub), recursive = TRUE,
               showWarnings = FALSE)
  }

  writeLines(paste(c(study_id, "syn", "Synthetic study", "Syn"),
                   collapse = "\t"),
             file.path(root, "Study", paste0(study_id, "_study.tsv")))

  # Cell tsv: 7-column form with cellID duplicated in col 1 & 2
  lines <- apply(cells, 1, function(r) {
    paste(c(r["cellID"], r["cellID"], r["cellType"], r["cellGroup"],
            r["coord1"], r["coord2"], "UMAP"),
          collapse = "\t")
  })
  writeLines(lines,
             file.path(root, "Cell", paste0(study_id, "_n_cell.tsv")))

  genes <- unique(c(rownames(gene_expression),
                    rownames(gene_activity_tf),
                    rownames(gene_activity_sig)))
  writeLines(genes,
             file.path(root, "Gene", paste0(study_id, "_n_gene.tsv")))

  # study_meta.csv
  meta_rows <- unique(cells[, c("cellType", "cellGroup")])
  meta_rows$StudyID  <- study_id
  meta_rows$StudyAbbr <- "syn"
  meta_rows$Color    <- "#aabbcc"
  meta_rows$Label_1  <- 0
  meta_rows$Label_2  <- 0
  meta_rows$NetworkCellType <- meta_rows$cellType
  meta_rows$Species  <- "Homo sapiens"
  utils::write.csv(meta_rows[, c("StudyID", "StudyAbbr", "cellType",
                                 "cellGroup", "Color", "Label_1",
                                 "Label_2", "NetworkCellType", "Species")],
                   file.path(root, "study_meta",
                             paste0(study_id, "_study_meta.csv")),
                   row.names = FALSE)
  # rename to match expected column casing
  meta_path <- file.path(root, "study_meta",
                         paste0(study_id, "_study_meta.csv"))
  txt <- readLines(meta_path)
  txt[1] <- gsub("\"cellType\"", "\"CellType\"", txt[1], fixed = TRUE)
  txt[1] <- gsub("\"cellGroup\"", "\"CellGroup\"", txt[1], fixed = TRUE)
  writeLines(txt, meta_path)

  # Networks (minimal)
  net_line <- function(s, t, ty, ct) {
    paste(c(s, t, ty, study_id, ct, "0.1", "0.2", "0.2", "0.2", "0.5"),
          collapse = "\t")
  }
  writeLines(c(net_line("g1", "g2", "TF", cells$cellType[1])),
             file.path(root, "Network_TF_Activity",
                       paste0(study_id, "_TF.tsv")))
  writeLines(c(net_line("g1", "g2", "SIG", cells$cellType[1])),
             file.path(root, "Network_SIG_Activity",
                       paste0(study_id, "_SIG.tsv")))

  # Shard tree + manifests. Layout under root/:
  #   expression_files/<studyID>/{meta.csv, <letter>/<gene>.csv.gz}
  #   activity_files/<studyID>/{meta.csv (shared), TF/<letter>/..., SIG/<letter>/...}
  write_shards <- function(mat, shard_root, meta_root, manifest_dir,
                            manifest_name, kind) {
    base <- file.path(root, shard_root)
    dir.create(base, recursive = TRUE, showWarnings = FALSE)
    meta_base <- file.path(root, meta_root)
    dir.create(meta_base, recursive = TRUE, showWarnings = FALSE)
    writeLines(paste(cell_order_in_shards, collapse = ","),
               file.path(meta_base, "meta.csv"))
    manifest <- data.frame(
      GeneSymbol = character(0), Species = character(0),
      StudyID = character(0), StudyAbbr = character(0),
      Type = character(0), FileHeader = character(0), File = character(0),
      stringsAsFactors = FALSE
    )
    perm <- match(cell_order_in_shards, colnames(mat))
    mat_shard <- mat[, perm, drop = FALSE]
    for (g in rownames(mat_shard)) {
      if (g %in% drop_genes) next
      first <- tolower(substr(g, 1, 1))
      dir_id <- if (grepl("[a-z]", first)) first else "nm"
      dir.create(file.path(base, dir_id), recursive = TRUE,
                 showWarnings = FALSE)
      csv_path <- file.path(base, dir_id, paste0(g, ".csv"))
      writeLines(paste(format(mat_shard[g, ], scientific = FALSE,
                              trim = TRUE),
                       collapse = ","), csv_path)
      gz_path <- paste0(csv_path, ".gz")
      R.utils::gzip(csv_path, destname = gz_path, overwrite = TRUE,
                    remove = TRUE)
      rel_file <- file.path(shard_root, dir_id, paste0(g, ".csv.gz"))
      manifest <- rbind(manifest, data.frame(
        GeneSymbol = g, Species = "Homo sapiens",
        StudyID = study_id, StudyAbbr = "syn", Type = kind,
        FileHeader = file.path(meta_root, "meta.csv"),
        File = rel_file,
        stringsAsFactors = FALSE
      ))
    }
    utils::write.csv(manifest,
                     file.path(root, manifest_dir,
                               paste0(study_id, "_", manifest_name, ".csv")),
                     row.names = FALSE)
  }
  exp_root <- file.path("expression_files", study_id)
  act_root <- file.path("activity_files",   study_id)
  write_shards(gene_expression,   shard_root = exp_root,
               meta_root = exp_root,
               manifest_dir = "study_gene_expression",
               manifest_name = "expression", kind = "Expression")
  write_shards(gene_activity_tf,  shard_root = file.path(act_root, "TF"),
               meta_root = act_root,
               manifest_dir = "study_gene_tf",
               manifest_name = "activity_tf", kind = "TF")
  write_shards(gene_activity_sig, shard_root = file.path(act_root, "SIG"),
               meta_root = act_root,
               manifest_dir = "study_gene_sig",
               manifest_name = "activity_sig", kind = "SIG")

  invisible(root)
}

# Skip the whole file if R.utils is not available (it is in DESCRIPTION of
# the original pipeline, used here only by tests).
if (!requireNamespace("R.utils", quietly = TRUE)) {
  test_that("shard tests skipped without R.utils", { skip("R.utils not installed") })
} else {

test_that("read_graph_study reconstructs expression matrix from shards", {
  root <- tempfile("shard_root_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  cells <- data.frame(
    cellID    = c("c1", "c2", "c3", "c4"),
    cellType  = c("A", "A", "B", "B"),
    cellGroup = c("A", "A", "B", "B"),
    coord1    = c(0.1, 0.2, 0.3, 0.4),
    coord2    = c(1.1, 1.2, 1.3, 1.4),
    stringsAsFactors = FALSE
  )
  expr <- matrix(c(0, 0.5, 0, 0,
                   0,   0, 2.5, 0,
                   1.1, 0, 0, 1.2),
                 nrow = 3, byrow = TRUE,
                 dimnames = list(c("Apc", "Bmp2", "Cdk1"),
                                 c("c1", "c2", "c3", "c4")))
  act_tf  <- expr * 2
  rownames(act_tf) <- c("Apc", "Bmp2", "Cdk1")
  act_sig <- expr * 3
  rownames(act_sig) <- c("Apc", "Bmp2", "Cdk1")

  make_fake_layout(root, study_id = "1234",
                   cells = cells,
                   gene_expression = expr,
                   gene_activity_tf = act_tf,
                   gene_activity_sig = act_sig)

  s <- read_graph_study(root, "1234",
                        load_expression = TRUE,
                        load_activity_tf = TRUE,
                        load_activity_sig = TRUE)

  expect_equal(s$meta$species, "Homo sapiens")
  expect_equal(nrow(s$cells), 4)
  expect_equal(s$cells$cellID, c("c1", "c2", "c3", "c4"))
  expect_equal(s$genes, c("Apc", "Bmp2", "Cdk1"))

  expect_equal(as.matrix(s$expression),   expr,        ignore_attr = TRUE)
  expect_equal(as.matrix(s$activity_tf),  expr * 2,    ignore_attr = TRUE)
  expect_equal(as.matrix(s$activity_sig), expr * 3,    ignore_attr = TRUE)
})

test_that("shard reader reorders columns when meta cell order differs", {
  root <- tempfile("shard_root_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  cells <- data.frame(
    cellID    = c("c1", "c2", "c3", "c4"),
    cellType  = c("A", "A", "B", "B"),
    cellGroup = c("A", "A", "B", "B"),
    coord1    = c(0.1, 0.2, 0.3, 0.4),
    coord2    = c(1.1, 1.2, 1.3, 1.4),
    stringsAsFactors = FALSE
  )
  expr <- matrix(c(0, 0.5, 0, 0,
                   0,   0, 2.5, 0,
                   1.1, 0, 0, 1.2),
                 nrow = 3, byrow = TRUE,
                 dimnames = list(c("Apc", "Bmp2", "Cdk1"),
                                 c("c1", "c2", "c3", "c4")))

  # Shards written in a different cell order than Cell/ TSV
  shuffled <- c("c3", "c1", "c4", "c2")
  make_fake_layout(root, study_id = "1234",
                   cells = cells,
                   gene_expression = expr,
                   gene_activity_tf = expr,
                   gene_activity_sig = expr,
                   cell_order_in_shards = shuffled)

  s <- read_graph_study(root, "1234",
                        load_expression = TRUE)
  # The reader must align columns back to the Cell/ order
  expect_equal(as.matrix(s$expression), expr, ignore_attr = TRUE)
})

test_that("missing shards leave gene rows zero with a warning", {
  root <- tempfile("shard_root_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  cells <- data.frame(
    cellID    = c("c1", "c2"),
    cellType  = c("A", "A"),
    cellGroup = c("A", "A"),
    coord1    = c(0, 1), coord2 = c(0, 1),
    stringsAsFactors = FALSE
  )
  expr <- matrix(c(1, 2, 3, 4),
                 nrow = 2, byrow = TRUE,
                 dimnames = list(c("Apc", "Bmp2"), c("c1", "c2")))

  make_fake_layout(root, "1234",
                   cells = cells,
                   gene_expression = expr,
                   gene_activity_tf = expr,
                   gene_activity_sig = expr,
                   drop_genes = "Bmp2")

  s <- read_graph_study(root, "1234", load_expression = TRUE)
  # The Bmp2 shard was deleted; manifest still lists it but file is absent.
  # Reader should produce a 2-row matrix with row 2 all zero.
  m <- as.matrix(s$expression)
  expect_equal(m[1, ], c(c1 = 1, c2 = 2), ignore_attr = TRUE)
  expect_equal(m[2, ], c(c1 = 0, c2 = 0), ignore_attr = TRUE)
})

test_that("missing manifest yields NULL matrix slot and a warning", {
  root <- tempfile("shard_root_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  cells <- data.frame(
    cellID = c("c1"), cellType = "A", cellGroup = "A",
    coord1 = 0, coord2 = 0, stringsAsFactors = FALSE
  )
  expr <- matrix(1, nrow = 1, ncol = 1,
                 dimnames = list("Apc", "c1"))
  make_fake_layout(root, "1234",
                   cells = cells,
                   gene_expression = expr,
                   gene_activity_tf = expr,
                   gene_activity_sig = expr)

  # Remove the expression manifest
  unlink(file.path(root, "study_gene_expression",
                   "1234_expression.csv"))

  expect_warning(
    s <- read_graph_study(root, "1234", load_expression = TRUE),
    "manifest not found"
  )
  expect_null(s$expression)
})

}  # end requireNamespace gate
