test_that("bundle round-trips metadata + indexes + networks", {
  s <- synthetic_study()
  tmp <- tempfile(fileext = ".scminer.h5")
  on.exit(unlink(tmp), add = TRUE)

  write_bundle(
    bundle_path        = tmp,
    meta               = s$meta,
    cells              = s$cells,
    clusters           = s$clusters,
    genes              = s$genes,
    expression_genes   = rownames(s$expression),
    activity_tf_genes  = rownames(s$activity_tf),
    activity_sig_genes = rownames(s$activity_sig),
    default_genes      = s$genes[1:3],
    network_tf         = s$network_tf,
    network_sig        = s$network_sig
  )
  expect_true(file.exists(tmp))

  out <- load_study(tmp)
  expect_s3_class(out, "scminer_study")
  expect_equal(out$meta$studyID,       s$meta$studyID)
  expect_equal(out$meta$bundleVersion, 1L)

  expect_equal(out$cells$cellID,    s$cells$cellID)
  expect_equal(out$cells$cellType,  s$cells$cellType)
  expect_equal(out$cells$coord1,    s$cells$coord1)

  ord <- order(out$clusters$cellType)
  ord_s <- order(s$clusters$cellType)
  expect_equal(out$clusters$cellType[ord], s$clusters$cellType[ord_s])
  expect_equal(out$clusters$count[ord],    s$clusters$count[ord_s])

  expect_equal(out$genes, s$genes)
  expect_equal(out$expression_index,   rownames(s$expression))
  expect_equal(out$activity_tf_index,  rownames(s$activity_tf))
  expect_equal(out$activity_sig_index, rownames(s$activity_sig))
  expect_equal(out$default_genes, s$genes[1:3])

  expect_equal(out$network_tf$source, s$network_tf$source)
  expect_equal(out$network_tf$mi,     s$network_tf$mi)
})

test_that("bundle round-trips with all optional groups omitted", {
  s <- synthetic_study()
  tmp <- tempfile(fileext = ".scminer.h5")
  on.exit(unlink(tmp), add = TRUE)

  write_bundle(
    bundle_path = tmp,
    meta        = s$meta,
    cells       = s$cells,
    clusters    = s$clusters,
    genes       = s$genes
  )

  out <- load_study(tmp)
  expect_null(out$expression_index)
  expect_null(out$activity_tf_index)
  expect_null(out$activity_sig_index)
  expect_null(out$default_genes)
  expect_null(out$network_tf)
  expect_null(out$network_sig)
  expect_equal(out$cells$cellID, s$cells$cellID)
  expect_equal(out$genes, s$genes)
})

test_that("write_bundle refuses to overwrite by default", {
  s <- synthetic_study()
  tmp <- tempfile(fileext = ".scminer.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_bundle(tmp, s$meta, s$cells, s$clusters, s$genes)
  expect_error(
    write_bundle(tmp, s$meta, s$cells, s$clusters, s$genes),
    "already exists"
  )
  expect_silent(
    write_bundle(tmp, s$meta, s$cells, s$clusters, s$genes,
                 overwrite = TRUE)
  )
})

test_that("write_bundle errors when network is missing required columns", {
  s <- synthetic_study()
  tmp <- tempfile(fileext = ".scminer.h5")
  on.exit(unlink(tmp), add = TRUE)
  bad_net <- s$network_tf[, c("source", "target", "mi")]
  expect_error(
    write_bundle(tmp, s$meta, s$cells, s$clusters, s$genes,
                 network_tf = bad_net),
    "missing columns"
  )
})

test_that("load_study infers shard_dir from bundle path", {
  s <- synthetic_study()
  bundle_dir <- tempfile("bundle_")
  dir.create(bundle_dir)
  on.exit(unlink(bundle_dir, recursive = TRUE), add = TRUE)
  bundle_path <- file.path(bundle_dir,
                           paste0(s$meta$studyID, ".scminer.h5"))
  write_bundle(bundle_path, s$meta, s$cells, s$clusters, s$genes)
  out <- load_study(bundle_path)
  expect_equal(normalizePath(out$shard_dir),
               normalizePath(bundle_dir))
})
