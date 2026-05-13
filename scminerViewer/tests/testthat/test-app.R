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
