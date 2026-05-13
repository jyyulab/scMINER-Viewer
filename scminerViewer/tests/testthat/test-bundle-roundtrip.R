test_that("bundle round-trips a synthetic study with full matrices", {
  s <- synthetic_study()
  tmp <- tempfile(fileext = ".scminer.h5")
  on.exit(unlink(tmp), add = TRUE)

  write_bundle(
    bundle_path  = tmp,
    meta         = s$meta,
    cells        = s$cells,
    clusters     = s$clusters,
    genes        = s$genes,
    expression   = s$expression,
    activity_tf  = s$activity_tf,
    activity_sig = s$activity_sig,
    network_tf   = s$network_tf,
    network_sig  = s$network_sig
  )
  expect_true(file.exists(tmp))

  out <- load_study(tmp)
  expect_s3_class(out, "scminer_study")

  expect_equal(out$meta$studyID,    s$meta$studyID)
  expect_equal(out$meta$studyAbbr,  s$meta$studyAbbr)
  expect_equal(out$meta$species,    s$meta$species)
  expect_equal(out$meta$coordinate, s$meta$coordinate)
  expect_equal(out$meta$bundleVersion, 1L)

  expect_equal(out$cells$cellID,    s$cells$cellID)
  expect_equal(out$cells$cellType,  s$cells$cellType)
  expect_equal(out$cells$cellGroup, s$cells$cellGroup)
  expect_equal(out$cells$coord1,    s$cells$coord1)
  expect_equal(out$cells$coord2,    s$cells$coord2)

  ord <- order(out$clusters$cellType)
  ord_s <- order(s$clusters$cellType)
  expect_equal(out$clusters$cellType[ord], s$clusters$cellType[ord_s])
  expect_equal(out$clusters$count[ord],    s$clusters$count[ord_s])
  expect_equal(out$clusters$color[ord],    s$clusters$color[ord_s])
  expect_equal(out$clusters$label_1[ord],  s$clusters$label_1[ord_s])

  expect_equal(out$genes, s$genes)

  expect_equal(as.matrix(out$expression),   as.matrix(s$expression),
               ignore_attr = TRUE)
  expect_equal(as.matrix(out$activity_tf),  as.matrix(s$activity_tf),
               ignore_attr = TRUE)
  expect_equal(as.matrix(out$activity_sig), as.matrix(s$activity_sig),
               ignore_attr = TRUE)

  expect_equal(out$network_tf$source,  s$network_tf$source)
  expect_equal(out$network_tf$target,  s$network_tf$target)
  expect_equal(out$network_tf$mi,      s$network_tf$mi)
  expect_equal(out$network_tf$pvalue,  s$network_tf$pvalue)
  expect_equal(nrow(out$network_sig),  nrow(s$network_sig))
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
  expect_null(out$expression)
  expect_null(out$activity_tf)
  expect_null(out$activity_sig)
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

test_that("write_bundle errors on gene/cell dimension mismatches", {
  s <- synthetic_study()
  tmp <- tempfile(fileext = ".scminer.h5")
  on.exit(unlink(tmp), add = TRUE)

  expect_error(
    write_bundle(tmp, s$meta, s$cells, s$clusters, s$genes,
                 expression = s$expression[1:5, ]),
    "expression"
  )
  expect_error(
    write_bundle(tmp, s$meta, s$cells, s$clusters, s$genes,
                 activity_tf = s$activity_tf[, 1:5]),
    "activity_tf"
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

test_that("sparse round-trip preserves an empty matrix", {
  s <- synthetic_study()
  tmp <- tempfile(fileext = ".scminer.h5")
  on.exit(unlink(tmp), add = TRUE)
  empty <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                dims = c(length(s$genes), nrow(s$cells)),
                                repr = "C")
  write_bundle(tmp, s$meta, s$cells, s$clusters, s$genes,
               expression = empty)
  out <- load_study(tmp)
  expect_equal(dim(out$expression), c(length(s$genes), nrow(s$cells)))
  expect_equal(length(out$expression@x), 0L)
})
