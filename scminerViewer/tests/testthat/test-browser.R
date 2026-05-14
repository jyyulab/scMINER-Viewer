# Tests for discover_studies() + build_browser(): multi-study root with
# the <studyID>/<studyID>.scminer.h5 layout written by prepare_study_data.

test_that("discover_studies returns empty data.frame for missing/empty dir", {
  d <- tempfile("empty_")
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  expect_equal(nrow(discover_studies(d)), 0L)
  dir.create(d)
  expect_equal(nrow(discover_studies(d)), 0L)
})

test_that("discover_studies finds all <studyID>/<studyID>.scminer.h5 bundles", {
  root <- tempfile("studies_")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  # Build two studies into the same root with different IDs.
  s1 <- synthetic_study(seed = 1)
  s1$meta$studyID    <- "1111"
  s1$meta$studyAbbr  <- "one"
  s1$meta$shortTitle <- "One"
  prepare_study_data(out_dir = root, meta = s1$meta,
                     cells = s1$cells, clusters = s1$clusters,
                     genes = s1$genes,
                     expression  = s1$expression,
                     activity_tf = s1$activity_tf,
                     activity_sig = s1$activity_sig,
                     network_tf = s1$network_tf,
                     network_sig = s1$network_sig,
                     emit = "bundle")

  s2 <- synthetic_study(seed = 2)
  s2$meta$studyID    <- "2222"
  s2$meta$studyAbbr  <- "two"
  s2$meta$shortTitle <- "Two"
  prepare_study_data(out_dir = root, meta = s2$meta,
                     cells = s2$cells, clusters = s2$clusters,
                     genes = s2$genes,
                     expression  = s2$expression,
                     emit = "bundle")

  df <- discover_studies(root)
  expect_equal(nrow(df), 2L)
  expect_setequal(df$studyID,   c("1111", "2222"))
  expect_setequal(df$studyAbbr, c("one",  "two"))
  expect_true(all(file.exists(df$bundle_path)))
  expect_true(all(df$n_cells > 0))
})

test_that("build_browser constructs a shiny app from a multi-study root", {
  skip_if_not_installed("bslib")
  skip_if_not_installed("plotly")
  skip_if_not_installed("DT")
  skip_if_not_installed("visNetwork")

  root <- tempfile("studies_")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  s <- synthetic_study(seed = 3)
  prepare_study_data(out_dir = root, meta = s$meta,
                     cells = s$cells, clusters = s$clusters,
                     genes = s$genes,
                     expression = s$expression,
                     emit = "bundle")

  app <- build_browser(root)
  expect_s3_class(app, "shiny.appobj")
})
