test_that("prepare_study_data emits graph layout and bundle round-trips", {
  s <- synthetic_study()
  out_dir <- tempfile("prep_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  res <- prepare_study_data(
    out_dir      = out_dir,
    meta         = s$meta,
    cells        = s$cells,
    clusters     = s$clusters,
    genes        = s$genes,
    expression   = s$expression,
    activity_tf  = s$activity_tf,
    activity_sig = s$activity_sig,
    network_tf   = s$network_tf,
    network_sig  = s$network_sig,
    emit         = c("graph", "bundle")
  )

  expect_true(file.exists(res$bundle_path))
  sid <- s$meta$studyID
  study_dir <- file.path(out_dir, sid)
  expect_equal(res$out_dir, study_dir)
  expect_equal(res$bundle_path,
               file.path(study_dir, paste0(sid, ".scminer.h5")))

  # --- graph layout (per-study subdir) -----------------------------------
  study_tsv <- file.path(study_dir, "Study", paste0(sid, "_study.tsv"))
  expect_true(file.exists(study_tsv))
  row <- strsplit(readLines(study_tsv), "\t", fixed = TRUE)[[1]]
  expect_equal(row, c(s$meta$studyID, s$meta$studyAbbr,
                      s$meta$longTitle, s$meta$shortTitle))

  cell_tsv <- file.path(study_dir, "Cell", paste0(sid, "_n_cell.tsv"))
  cell_lines <- readLines(cell_tsv)
  expect_equal(length(cell_lines), nrow(s$cells))
  first <- strsplit(cell_lines[1], "\t", fixed = TRUE)[[1]]
  expect_length(first, 7L)
  expect_equal(first[1], first[2])          # cellID duplicated
  expect_equal(first[1], s$cells$cellID[1])
  expect_equal(first[7], s$meta$coordinate)

  gene_lines <- readLines(file.path(study_dir, "Gene",
                                    paste0(sid, "_n_gene.tsv")))
  expect_equal(gene_lines, s$genes)

  cluster_csv <- file.path(study_dir, "study_meta",
                           paste0(sid, "_study_meta.csv"))
  expect_true(file.exists(cluster_csv))
  cluster_df <- utils::read.csv(cluster_csv, stringsAsFactors = FALSE)
  expect_setequal(cluster_df$CellType, s$clusters$cellType)
  expect_setequal(cluster_df$Color,    s$clusters$color)

  tf_tsv <- file.path(study_dir, "Network_TF_Activity",
                      paste0(sid, "_TF.tsv"))
  tf_lines <- readLines(tf_tsv)
  expect_equal(length(tf_lines), nrow(s$network_tf))
  fields <- strsplit(tf_lines[1], "\t", fixed = TRUE)[[1]]
  expect_length(fields, 10L)
  expect_equal(fields[3], "TF")
  expect_equal(fields[4], sid)

  # Manifests + shards
  exp_manifest <- file.path(study_dir, "study_gene_expression",
                            paste0(sid, "_expression.csv"))
  expect_true(file.exists(exp_manifest))
  manifest <- utils::read.csv(exp_manifest, stringsAsFactors = FALSE)
  expect_equal(nrow(manifest), length(s$genes))
  expect_setequal(manifest$GeneSymbol, s$genes)

  meta_csv <- file.path(study_dir, "expression_files", sid, "meta.csv")
  expect_true(file.exists(meta_csv))
  shard_cells <- strsplit(readLines(meta_csv, warn = FALSE),
                          ",", fixed = TRUE)[[1]]
  expect_equal(shard_cells, s$cells$cellID)

  # Activity TF and SIG share a single meta.csv one level above their kind dirs.
  expect_true(file.exists(file.path(study_dir, "activity_files", sid,
                                    "meta.csv")))

  # Spot-check one shard
  sample_gene <- s$genes[5]
  letter <- tolower(substr(sample_gene, 1, 1))
  if (!grepl("^[a-z]$", letter)) letter <- "nm"
  shard_path <- file.path(study_dir, "expression_files", sid, letter,
                          paste0(sample_gene, ".csv.gz"))
  expect_true(file.exists(shard_path))
  shard_vals <- as.numeric(strsplit(readLines(shard_path, warn = FALSE),
                                    ",", fixed = TRUE)[[1]])
  expect_equal(shard_vals,
               as.numeric(s$expression[sample_gene, ]))

  # --- Bundle stores indexes; values come lazily from disk ----------------
  loaded <- load_study(res$bundle_path)
  expect_equal(loaded$meta$studyID, s$meta$studyID)
  expect_equal(loaded$cells$cellID, s$cells$cellID)
  expect_equal(loaded$genes, s$genes)
  expect_setequal(loaded$expression_index,   rownames(s$expression))
  expect_setequal(loaded$activity_tf_index,  rownames(s$activity_tf))
  expect_setequal(loaded$activity_sig_index, rownames(s$activity_sig))

  # Lazy read of a single gene
  vals <- gene_values(loaded, sample_gene, "Express_normalized")
  expect_equal(vals,
               as.numeric(s$expression[sample_gene, ]),
               ignore_attr = TRUE)

  # --- Re-read the graph layout via read_graph_study (indexes only) -------
  # read_graph_study autodetects the wrapped layout under <root>/<sid>/.
  reread <- read_graph_study(out_dir, s$meta$studyID)
  expect_equal(reread$cells$cellID, s$cells$cellID)
  expect_equal(reread$genes,        s$genes)
  expect_equal(nrow(reread$network_tf),  nrow(s$network_tf))
  expect_equal(nrow(reread$network_sig), nrow(s$network_sig))
  expect_setequal(reread$expression_genes,   rownames(s$expression))
  expect_setequal(reread$activity_tf_genes,  rownames(s$activity_tf))
  expect_setequal(reread$activity_sig_genes, rownames(s$activity_sig))
})

test_that("prepare_study_data emits only bundle when requested", {
  s <- synthetic_study()
  out_dir <- tempfile("prep_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  res <- prepare_study_data(
    out_dir = out_dir, meta = s$meta, cells = s$cells,
    clusters = s$clusters, genes = s$genes,
    expression = s$expression,
    emit = "bundle"
  )
  sid <- s$meta$studyID
  expect_true(file.exists(res$bundle_path))
  expect_false(dir.exists(file.path(out_dir, sid, "Cell")))
  expect_false(dir.exists(file.path(out_dir, sid,
                                     "expression_files", sid)))
})

test_that("prepare_study_data emits only graph layout when requested", {
  s <- synthetic_study()
  out_dir <- tempfile("prep_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  res <- prepare_study_data(
    out_dir = out_dir, meta = s$meta, cells = s$cells,
    clusters = s$clusters, genes = s$genes,
    expression = s$expression,
    emit = "graph"
  )
  sid <- s$meta$studyID
  expect_null(res$bundle_path)
  expect_true(file.exists(file.path(out_dir, sid, "Cell",
                                     paste0(sid, "_n_cell.tsv"))))
  expect_true(file.exists(file.path(out_dir, sid,
                                     "expression_files", sid,
                                     "meta.csv")))
})

test_that("prepare_study_data fills in cluster counts when missing", {
  s <- synthetic_study()
  out_dir <- tempfile("prep_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  no_counts <- s$clusters
  no_counts$count <- NA_integer_
  res <- prepare_study_data(
    out_dir = out_dir, meta = s$meta, cells = s$cells,
    clusters = no_counts, genes = s$genes,
    emit = "bundle"
  )
  loaded <- load_study(res$bundle_path)
  expected <- as.data.frame(table(s$cells$cellType),
                            stringsAsFactors = FALSE)
  for (i in seq_len(nrow(loaded$clusters))) {
    ct <- loaded$clusters$cellType[i]
    expect_equal(loaded$clusters$count[i],
                 expected$Freq[expected$Var1 == ct])
  }
})
