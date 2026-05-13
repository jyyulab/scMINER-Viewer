test_that("read_graph_study loads 2327 from the project data folder", {
  data_dir <- Sys.getenv("SCMINER_DATA_DIR", unset = "")
  if (!nzchar(data_dir) || !dir.exists(data_dir)) {
    candidates <- c(
      file.path(testthat::test_path(), "..", "..", "..", "data"),
      file.path(getwd(), "..", "..", "..", "data")
    )
    for (cand in candidates) {
      if (dir.exists(cand)) { data_dir <- cand; break }
    }
  }
  skip_if(!nzchar(data_dir) || !dir.exists(data_dir),
          "data/ folder not found; set SCMINER_DATA_DIR")

  s <- read_graph_study(data_dir, "2327")
  expect_equal(s$meta$studyID, "2327")
  expect_equal(s$meta$studyAbbr, "tex")
  expect_equal(s$meta$coordinate, "UMAP")
  expect_true(nrow(s$cells) > 0)
  expect_true(length(s$genes) > 0)
  expect_true(nrow(s$clusters) > 0)
  expect_true(all(s$clusters$count >= 0))
  expect_true(sum(s$clusters$count) <= nrow(s$cells))
  expect_true(!is.null(s$network_tf) && nrow(s$network_tf) > 0)
  expect_true(!is.null(s$network_sig) && nrow(s$network_sig) > 0)
  expect_true(all(c("source", "target", "cellType", "mi") %in%
                  colnames(s$network_tf)))
})

test_that("read_graph_study output writes and reads back through the bundle", {
  data_dir <- Sys.getenv("SCMINER_DATA_DIR", unset = "")
  if (!nzchar(data_dir) || !dir.exists(data_dir)) {
    candidates <- c(
      file.path(testthat::test_path(), "..", "..", "..", "data"),
      file.path(getwd(), "..", "..", "..", "data")
    )
    for (cand in candidates) {
      if (dir.exists(cand)) { data_dir <- cand; break }
    }
  }
  skip_if(!nzchar(data_dir) || !dir.exists(data_dir),
          "data/ folder not found; set SCMINER_DATA_DIR")

  s <- read_graph_study(data_dir, "2327")
  tmp <- tempfile(fileext = ".scminer.h5")
  on.exit(unlink(tmp), add = TRUE)

  write_bundle(
    bundle_path = tmp,
    meta        = s$meta,
    cells       = s$cells,
    clusters    = s$clusters,
    genes       = s$genes,
    network_tf  = s$network_tf,
    network_sig = s$network_sig
  )

  out <- load_study(tmp)
  expect_equal(out$meta$studyID, "2327")
  expect_equal(nrow(out$cells), nrow(s$cells))
  expect_equal(length(out$genes), length(s$genes))
  expect_equal(nrow(out$clusters), nrow(s$clusters))
  expect_equal(nrow(out$network_tf), nrow(s$network_tf))
  expect_equal(nrow(out$network_sig), nrow(s$network_sig))
  expect_null(out$expression)
})
