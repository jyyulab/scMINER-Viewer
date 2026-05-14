# Tests for the staged helpers introduced for granular control / debug:
#   load_study_config(), extract_cells(), extract_genes(),
#   extract_expression(), extract_activity(), read_networks().

test_that("load_study_config fills in defaults", {
  skip_if_not_installed("yaml")
  tmp <- tempfile(fileext = ".yml")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(
    "output: out",
    "study:",
    "  ID: \"42\"",
    "  studyAbbr: tst",
    "  longTitle: A test study",
    "  shortTitle: Tst",
    "input:",
    "  expression: /tmp/expr.rds"
  ), tmp)

  cfg <- load_study_config(tmp)
  expect_equal(cfg$output, "out")
  expect_equal(cfg$study$ID, "42")
  # Defaults
  expect_equal(cfg$species, "")
  expect_equal(cfg$coordinate, "UMAP")
  expect_equal(cfg$cellID, "cellID")
  expect_equal(cfg$cellType, "cellGroup")
  expect_equal(cfg$cellGroup, "cellGroup")  # falls back to cellType
  expect_equal(cfg$geneSymbol, "geneSymbol")
})

test_that("load_study_config preserves explicit values", {
  skip_if_not_installed("yaml")
  tmp <- tempfile(fileext = ".yml")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(
    "output: out",
    "species: Homo sapiens",
    "coordinate: TSNE",
    "cellID: barcode",
    "cellType: ctype",
    "cellGroup: cgroup",
    "geneSymbol: gene_name",
    "study:",
    "  ID: \"7\"",
    "  studyAbbr: hs",
    "  longTitle: A second study",
    "  shortTitle: HS",
    "input:",
    "  expression: /a.rds",
    "  activity:   /b.rds",
    "  networks:   /c.tsv"
  ), tmp)
  cfg <- load_study_config(tmp)
  expect_equal(cfg$species, "Homo sapiens")
  expect_equal(cfg$coordinate, "TSNE")
  expect_equal(cfg$cellID, "barcode")
  expect_equal(cfg$cellType, "ctype")
  expect_equal(cfg$cellGroup, "cgroup")
  expect_equal(cfg$geneSymbol, "gene_name")
  expect_equal(cfg$input$activity, "/b.rds")
  expect_equal(cfg$input$networks, "/c.tsv")
})

test_that("load_study_config rejects missing required keys", {
  skip_if_not_installed("yaml")
  tmp <- tempfile(fileext = ".yml")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("output: out", tmp)
  expect_error(load_study_config(tmp), "missing required top-level key: 'study'")

  writeLines(c("output: out",
               "study:",
               "  ID: \"1\"",
               "input:",
               "  expression: /x.rds"), tmp)
  expect_error(load_study_config(tmp), "study.studyAbbr")
})

test_that("load_study_config errors on missing file", {
  expect_error(load_study_config(tempfile(fileext = ".yml")),
               "Config file not found")
})

# ----- Biobase-dependent extract_* tests -----------------------------------
# These need a real Biobase ExpressionSet. Skip cleanly if unavailable.

build_fake_eset <- function(n_cells = 6, n_genes = 4, with_activity = FALSE) {
  cell_ids <- paste0("c", seq_len(n_cells))
  gene_ids <- paste0("g", seq_len(n_genes))
  expr_mat <- matrix(seq_len(n_cells * n_genes),
                     nrow = n_genes, ncol = n_cells,
                     dimnames = list(gene_ids, cell_ids))

  p_data <- data.frame(
    cellGroup = rep(c("A", "B"), length.out = n_cells),
    UMAP_1    = seq_len(n_cells) * 1.0,
    UMAP_2    = seq_len(n_cells) * 0.5,
    row.names = cell_ids,
    stringsAsFactors = FALSE
  )
  f_data <- data.frame(
    geneSymbol = gene_ids,
    row.names  = gene_ids,
    stringsAsFactors = FALSE
  )

  eset <- Biobase::ExpressionSet(
    assayData = expr_mat,
    phenoData = methods::new("AnnotatedDataFrame", data = p_data),
    featureData = methods::new("AnnotatedDataFrame", data = f_data)
  )
  if (with_activity) {
    # Activity has rows named gN_TF and gN_SIG
    act_rows <- c(paste0(gene_ids, "_TF"), paste0(gene_ids, "_SIG"))
    act_mat <- matrix(seq_along(act_rows) * 0.1,
                      nrow = length(act_rows), ncol = n_cells,
                      dimnames = list(act_rows, cell_ids))
    act <- Biobase::ExpressionSet(
      assayData = act_mat,
      phenoData = methods::new("AnnotatedDataFrame", data = p_data),
      featureData = methods::new("AnnotatedDataFrame",
                                  data = data.frame(row.names = act_rows))
    )
    return(list(expression = eset, activity = act))
  }
  eset
}

test_that("extract_cells / genes / expression handle a typical ExpressionSet", {
  skip_if_not_installed("Biobase")
  eset <- build_fake_eset()

  cells <- extract_cells(eset)
  expect_equal(nrow(cells), 6)
  expect_equal(colnames(cells),
               c("cellID", "cellType", "cellGroup", "coord1", "coord2"))
  expect_equal(cells$cellID, paste0("c", 1:6))
  expect_equal(cells$cellType, rep(c("A", "B"), 3))

  genes <- extract_genes(eset)
  expect_equal(genes, paste0("g", 1:4))

  expr <- extract_expression(eset, genes = genes)
  expect_equal(dim(expr), c(4, 6))
  expect_equal(rownames(expr), genes)
})

test_that("extract_cells fails loudly on missing pData columns", {
  skip_if_not_installed("Biobase")
  eset <- build_fake_eset()
  expect_error(
    extract_cells(eset, coordinate_col = "PCA"),
    "PCA_1"
  )
  expect_error(
    extract_cells(eset, cell_type_col = "missing_col"),
    "missing_col"
  )
})

test_that("extract_activity splits by _TF / _SIG suffix", {
  skip_if_not_installed("Biobase")
  pair <- build_fake_eset(with_activity = TRUE)
  genes <- extract_genes(pair$expression)
  act   <- extract_activity(pair$activity, master_genes = genes)
  expect_true(!is.null(act$tf))
  expect_true(!is.null(act$sig))
  expect_equal(nrow(act$tf),  length(genes))
  expect_equal(nrow(act$sig), length(genes))
  # First TF row corresponds to "g1_TF" with value 0.1
  expect_equal(as.numeric(act$tf[1, ])[1], 0.1)
})

test_that("extract_activity warns when no _TF/_SIG rows present", {
  skip_if_not_installed("Biobase")
  eset <- build_fake_eset()      # rows named g1..g4, no suffix
  expect_warning(
    res <- extract_activity(eset, master_genes = paste0("g", 1:4)),
    "_TF/.TF or _SIG/.SIG"
  )
  expect_null(res$tf)
  expect_null(res$sig)
})

test_that("read_networks splits by NetworkType", {
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c(
    "source\ttarget\tNetworkType\tCellGroup\tmi\tpearson\tspearman\trho\tpvalue",
    "g1\tg2\tTF\tA\t0.5\t0.4\t0.4\t0.4\t0.01",
    "g1\tg3\tTF\tB\t0.3\t0.2\t0.2\t0.2\t0.10",
    "g2\tg3\tSIG\tA\t0.7\t0.6\t0.6\t0.6\t0.001"
  ), tmp)

  nets <- read_networks(tmp)
  expect_equal(nrow(nets$tf),  2L)
  expect_equal(nrow(nets$sig), 1L)
  expect_equal(nets$tf$source, c("g1", "g1"))
  expect_equal(nets$sig$cellType, "A")
  expect_equal(nets$sig$mi, 0.7)
})
