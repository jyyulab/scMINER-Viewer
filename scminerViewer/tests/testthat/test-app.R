test_that("build_app constructs a shinyApp from a bundle", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("plotly")
  skip_if_not_installed("DT")
  skip_if_not_installed("visNetwork")

  s <- synthetic_study()
  tmp <- tempfile(fileext = ".scminer.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_bundle(
    bundle_path  = tmp,
    meta         = s$meta, cells = s$cells, clusters = s$clusters,
    genes        = s$genes, expression = s$expression,
    activity_tf  = s$activity_tf, activity_sig = s$activity_sig,
    network_tf   = s$network_tf, network_sig  = s$network_sig
  )

  app <- build_app(tmp)
  expect_s3_class(app, "shiny.appobj")
})

test_that("plot helpers produce plotly objects from a loaded study", {
  skip_if_not_installed("plotly")
  s <- synthetic_study()
  tmp <- tempfile(fileext = ".scminer.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_bundle(
    bundle_path  = tmp,
    meta         = s$meta, cells = s$cells, clusters = s$clusters,
    genes        = s$genes, expression = s$expression,
    activity_tf  = s$activity_tf, activity_sig = s$activity_sig,
    network_tf   = s$network_tf, network_sig  = s$network_sig
  )
  study <- load_study(tmp)

  p_cluster <- scminerViewer:::.cluster_plot(
    study, active_clusters = study$clusters$cellType,
    dot_size = 3, show_labels = TRUE
  )
  expect_s3_class(p_cluster, "plotly")

  g <- study$genes[1]
  p_feat <- scminerViewer:::.feature_plot(
    study, gene = g, relationship = "Express_normalized",
    active_clusters = study$clusters$cellType, dot_size = 3
  )
  expect_s3_class(p_feat, "plotly")

  p_vio <- scminerViewer:::.violin_plot(
    study, gene = g, relationship = "Express_normalized",
    active_clusters = study$clusters$cellType
  )
  expect_s3_class(p_vio, "plotly")

  p_heat <- scminerViewer:::.heatmap_plot(
    study, genes = study$genes[1:5],
    relationship = "Express_normalized",
    active_clusters = study$clusters$cellType
  )
  expect_s3_class(p_heat, "plotly")

  p_bub <- scminerViewer:::.bubble_plot(
    study, genes = study$genes[1:5],
    relationship = "Express_normalized",
    active_clusters = study$clusters$cellType
  )
  expect_s3_class(p_bub, "plotly")
})

test_that("network plot surfaces MI + Direction in hover and legend", {
  skip_if_not_installed("visNetwork")
  s <- synthetic_study()
  tmp <- tempfile(fileext = ".scminer.h5")
  on.exit(unlink(tmp), add = TRUE)
  write_bundle(
    bundle_path  = tmp,
    meta         = s$meta, cells = s$cells, clusters = s$clusters,
    genes        = s$genes, expression = s$expression,
    activity_tf  = s$activity_tf, activity_sig = s$activity_sig,
    network_tf   = s$network_tf, network_sig  = s$network_sig
  )
  study <- load_study(tmp)
  hub <- as.character(study$network_tf$source[1])

  plt <- scminerViewer:::.network_plot(
    study, gene = hub, network_type = "TF",
    active_clusters = study$clusters$cellType
  )
  expect_s3_class(plt, "visNetwork")

  # Focus gene pinned at the origin; neighbors on a circle around it.
  pn <- plt$x$nodes
  focus_row <- pn[pn$id == hub, , drop = FALSE]
  expect_equal(nrow(focus_row), 1L)
  expect_equal(focus_row$x, 0)
  expect_equal(focus_row$y, 0)
  expect_true(isTRUE(focus_row$fixed))
  neighbors <- pn[pn$id != hub, , drop = FALSE]
  if (nrow(neighbors) > 0L) {
    radii <- sqrt(neighbors$x^2 + neighbors$y^2)
    expect_true(all(abs(radii - radii[1]) < 1e-6))  # same radius for all
  }

  # Pull the underlying edge data from the htmlwidget payload and assert
  # the hover title and direction encoding are present.
  edges <- plt$x$edges
  expect_true(all(c("title", "color", "arrows") %in% colnames(edges)))

  # Arrows must land at the *neighbor* end. Visual edge orientation is
  # always focus -> neighbor regardless of underlying source/target, so
  # `from` should equal the focus gene for every edge.
  expect_true(all(edges$from == hub))
  expect_false(any(edges$to == hub))
  joined_title <- paste(edges$title, collapse = " | ")
  expect_match(joined_title, "MI:")
  expect_match(joined_title, "Direction:")
  expect_match(joined_title, "Spearman:|Pearson:")
  expect_true(any(grepl("Activator|Repressor|Unknown", edges$title)))
  expect_true(all(edges$arrows == "to"))

  # Legend definition carries the two direction entries.
  expect_true(!is.null(plt$x$legend))
})
