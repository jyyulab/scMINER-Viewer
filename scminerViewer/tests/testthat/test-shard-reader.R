# Build a tiny synthetic data layout in tempdir, then verify the lazy
# gene_values() accessor reads the right shard, aligns columns to
# cells$cellID, caches results, and handles missing shards.

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

  meta_rows <- unique(cells[, c("cellType", "cellGroup")])
  meta_rows$StudyID <- study_id
  meta_rows$StudyAbbr <- "syn"
  meta_rows$Color <- "#aabbcc"
  meta_rows$Label_1 <- 0
  meta_rows$Label_2 <- 0
  meta_rows$NetworkCellType <- meta_rows$cellType
  meta_rows$Species <- "Homo sapiens"
  utils::write.csv(
    meta_rows[, c("StudyID", "StudyAbbr", "cellType", "cellGroup",
                  "Color", "Label_1", "Label_2", "NetworkCellType",
                  "Species")],
    file.path(root, "study_meta",
              paste0(study_id, "_study_meta.csv")),
    row.names = FALSE
  )
  meta_path <- file.path(root, "study_meta",
                         paste0(study_id, "_study_meta.csv"))
  txt <- readLines(meta_path)
  txt[1] <- gsub("\"cellType\"", "\"CellType\"", txt[1], fixed = TRUE)
  txt[1] <- gsub("\"cellGroup\"", "\"CellGroup\"", txt[1], fixed = TRUE)
  writeLines(txt, meta_path)

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

if (!requireNamespace("R.utils", quietly = TRUE)) {
  test_that("shard tests skipped without R.utils", { skip("R.utils not installed") })
} else {

# Helper: build the synthetic data + bundle and return the loaded study.
build_synthetic <- function(seed = 7, ...) {
  set.seed(seed)
  root <- tempfile("layout_")
  dir.create(root)

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
  act_tf  <- expr * 2; rownames(act_tf)  <- rownames(expr)
  act_sig <- expr * 3; rownames(act_sig) <- rownames(expr)
  make_fake_layout(root, "1234",
                   cells = cells,
                   gene_expression = expr,
                   gene_activity_tf = act_tf,
                   gene_activity_sig = act_sig, ...)
  list(root = root, cells = cells, expr = expr,
       act_tf = act_tf, act_sig = act_sig)
}

test_that("gene_values() lazily reads expression / activity shards", {
  fx <- build_synthetic()
  on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
  s <- read_graph_study(fx$root, "1234")
  bundle <- file.path(fx$root, "1234.scminer.h5")
  write_bundle(
    bundle, meta = s$meta, cells = s$cells, clusters = s$clusters,
    genes = s$genes,
    expression_genes   = s$expression_genes,
    activity_tf_genes  = s$activity_tf_genes,
    activity_sig_genes = s$activity_sig_genes,
    network_tf = s$network_tf, network_sig = s$network_sig
  )
  study <- load_study(bundle)

  expect_setequal(study$expression_index, rownames(fx$expr))

  for (g in rownames(fx$expr)) {
    vals_exp <- gene_values(study, g, "Express_normalized")
    expect_equal(vals_exp, as.numeric(fx$expr[g, ]), ignore_attr = TRUE)
    vals_tf <- gene_values(study, g, "Activity_tf")
    expect_equal(vals_tf, as.numeric(fx$act_tf[g, ]), ignore_attr = TRUE)
    vals_sig <- gene_values(study, g, "Activity_sig")
    expect_equal(vals_sig, as.numeric(fx$act_sig[g, ]), ignore_attr = TRUE)
  }

  # Cached on repeat call
  first  <- gene_values(study, "Apc")
  second <- gene_values(study, "Apc")
  expect_identical(first, second)
})

test_that("gene_values() returns NULL for unknown gene", {
  fx <- build_synthetic()
  on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
  s <- read_graph_study(fx$root, "1234")
  bundle <- file.path(fx$root, "1234.scminer.h5")
  write_bundle(
    bundle, meta = s$meta, cells = s$cells, clusters = s$clusters,
    genes = s$genes,
    expression_genes = s$expression_genes
  )
  study <- load_study(bundle)
  expect_null(gene_values(study, "NotAGene"))
})

test_that("gene_values() handles shuffled shard cell order", {
  fx <- build_synthetic(cell_order_in_shards = c("c3", "c1", "c4", "c2"))
  on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
  s <- read_graph_study(fx$root, "1234")
  bundle <- file.path(fx$root, "1234.scminer.h5")
  write_bundle(
    bundle, meta = s$meta, cells = s$cells, clusters = s$clusters,
    genes = s$genes,
    expression_genes = s$expression_genes
  )
  study <- load_study(bundle)
  vals <- gene_values(study, "Apc", "Express_normalized")
  expect_equal(vals, as.numeric(fx$expr["Apc", ]), ignore_attr = TRUE)
})

test_that("missing shard file yields NULL not an error", {
  fx <- build_synthetic(drop_genes = "Bmp2")
  on.exit(unlink(fx$root, recursive = TRUE), add = TRUE)
  s <- read_graph_study(fx$root, "1234")
  bundle <- file.path(fx$root, "1234.scminer.h5")
  write_bundle(
    bundle, meta = s$meta, cells = s$cells, clusters = s$clusters,
    genes = s$genes,
    expression_genes = s$expression_genes
  )
  study <- load_study(bundle)
  # Bmp2 is in the manifest/index, but its shard was deleted
  expect_null(gene_values(study, "Bmp2", "Express_normalized"))
  # Apc still works
  expect_equal(gene_values(study, "Apc"),
               as.numeric(fx$expr["Apc", ]),
               ignore_attr = TRUE)
})

}  # end requireNamespace gate
