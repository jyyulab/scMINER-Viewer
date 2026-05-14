# Tests for fill_clusters() — auto-population of cluster colours
# (via ggsci) and label_1 / label_2 centroids (mean of cell coords).

test_that("fill_clusters builds clusters from cells when none provided", {
  cells <- data.frame(
    cellID    = paste0("c", 1:6),
    cellType  = rep(c("A", "B"), each = 3),
    cellGroup = rep(c("A", "B"), each = 3),
    coord1    = c(0, 1, 2, 10, 11, 12),
    coord2    = c(0, 0, 0, 5, 5, 5),
    stringsAsFactors = FALSE
  )
  clusters <- fill_clusters(cells)
  expect_equal(sort(clusters$cellType), c("A", "B"))
  expect_equal(clusters$count, c(3L, 3L))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", clusters$color)))
  # label_1 should be the per-cluster mean coord1
  expect_equal(clusters$label_1[clusters$cellType == "A"], 1)
  expect_equal(clusters$label_1[clusters$cellType == "B"], 11)
  expect_equal(clusters$label_2[clusters$cellType == "A"], 0)
  expect_equal(clusters$label_2[clusters$cellType == "B"], 5)
})

test_that("fill_clusters preserves existing colors and labels", {
  cells <- data.frame(
    cellID    = paste0("c", 1:4),
    cellType  = c("A", "A", "B", "B"),
    cellGroup = c("A", "A", "B", "B"),
    coord1    = c(0, 2, 8, 10),
    coord2    = c(0, 0, 5, 5),
    stringsAsFactors = FALSE
  )
  existing <- data.frame(
    cellType = c("A", "B"),
    count    = c(2L, 2L),
    color    = c("#123456", "#abcdef"),
    label_1  = c(99, 99),
    label_2  = c(-1, -1),
    stringsAsFactors = FALSE
  )
  out <- fill_clusters(cells, existing)
  expect_equal(out$color, c("#123456", "#abcdef"))
  expect_equal(out$label_1, c(99, 99))
  expect_equal(out$label_2, c(-1, -1))
})

test_that("fill_clusters tops up partially-present color / label cells", {
  cells <- data.frame(
    cellID    = paste0("c", 1:4),
    cellType  = c("A", "A", "B", "B"),
    cellGroup = c("A", "A", "B", "B"),
    coord1    = c(0, 2, 8, 10),
    coord2    = c(0, 0, 5, 5),
    stringsAsFactors = FALSE
  )
  partial <- data.frame(
    cellType = c("A", "B"),
    count    = c(2L, 2L),
    color    = c("#123456", NA_character_),
    label_1  = c(NA_real_, 99),
    label_2  = c(NA_real_, NA_real_),
    stringsAsFactors = FALSE
  )
  out <- fill_clusters(cells, partial)
  expect_equal(out$color[1], "#123456")        # preserved
  expect_match(out$color[2], "^#[0-9A-Fa-f]{6}$")  # filled
  expect_equal(out$label_1[1], 1)              # centroid for A
  expect_equal(out$label_1[2], 99)             # preserved for B
  expect_equal(out$label_2, c(0, 5))           # both auto-filled
})

test_that("fill_clusters respects the palette argument", {
  cells <- data.frame(
    cellID    = paste0("c", 1:4),
    cellType  = c("A", "B", "C", "D"),
    cellGroup = c("A", "B", "C", "D"),
    coord1    = 1:4, coord2 = 1:4,
    stringsAsFactors = FALSE
  )
  npg_cols  <- fill_clusters(cells, palette = "npg")$color
  jama_cols <- fill_clusters(cells, palette = "jama")$color
  if (requireNamespace("ggsci", quietly = TRUE)) {
    expect_false(identical(npg_cols, jama_cols))
  }
  # Always 7-char hex (no alpha)
  expect_true(all(nchar(npg_cols) == 7L))
})

test_that("fill_clusters warns on unknown palette and falls back to npg", {
  cells <- data.frame(
    cellID    = paste0("c", 1:3),
    cellType  = c("A", "B", "C"),
    cellGroup = c("A", "B", "C"),
    coord1    = 1:3, coord2 = 1:3,
    stringsAsFactors = FALSE
  )
  if (requireNamespace("ggsci", quietly = TRUE)) {
    expect_warning(
      out <- fill_clusters(cells, palette = "not_a_palette"),
      "Unknown cluster palette"
    )
    expect_match(out$color, "^#[0-9A-Fa-f]{6}$")
  } else {
    skip("ggsci not installed")
  }
})

test_that("prepare_study_data emits a study_meta CSV with real colors + labels", {
  s <- synthetic_study()
  out_dir <- tempfile("prep_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  # Pass clusters WITHOUT colors/labels to force auto-fill
  s$clusters$color   <- NULL
  s$clusters$label_1 <- NULL
  s$clusters$label_2 <- NULL

  res <- prepare_study_data(
    out_dir         = out_dir,
    meta            = s$meta,
    cells           = s$cells,
    clusters        = s$clusters,
    genes           = s$genes,
    cluster_palette = "jama",
    emit            = c("graph", "bundle")
  )
  sid <- s$meta$studyID
  meta_csv <- file.path(out_dir, sid, "study_meta",
                        paste0(sid, "_study_meta.csv"))
  expect_true(file.exists(meta_csv))
  df <- utils::read.csv(meta_csv, stringsAsFactors = FALSE)

  # Every colour is a real 7-char hex, not the old "#888888" default
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", df$Color)))
  expect_false(any(df$Color == "#888888"))

  # Labels are non-zero centroids drawn from the synthetic coord1/coord2
  expect_true(any(df$Label_1 != 0))

  # And the bundle's clusters round-trip the same values
  loaded <- load_study(res$bundle_path)
  expect_equal(sort(loaded$clusters$color), sort(df$Color))
  expect_equal(sort(loaded$clusters$label_1), sort(df$Label_1))
})
