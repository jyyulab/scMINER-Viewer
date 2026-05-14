# Plot helpers used by the Shiny server. Each accepts a list with the
# study and the current UI state, and returns a plotly object (or
# visNetwork object for the network panel).

# Shared color palette: rotates through a Tableau-like palette when the
# bundle doesn't provide cluster colors.
.default_colors <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"
)

.cluster_color_map <- function(clusters) {
  cols <- as.character(clusters$color)
  if (any(is.na(cols)) || any(!nzchar(cols))) {
    cols <- rep_len(.default_colors, nrow(clusters))
  }
  stats::setNames(cols, clusters$cellType)
}

.cells_in_active_clusters <- function(study, active_clusters,
                                       cell_mask = NULL) {
  base <- if (is.null(active_clusters) || length(active_clusters) == 0) {
    rep(TRUE, nrow(study$cells))
  } else {
    study$cells$cellType %in% active_clusters
  }
  if (!is.null(cell_mask)) base <- base & as.logical(cell_mask)
  base
}

# A no-data plotly figure that doesn't emit the "no trace type specified"
# warning that plotly::plotly_empty() does. We give plotly a real (but
# empty) scatter trace plus an annotation as the placeholder message.
.empty_plot <- function(title = NULL) {
  fig <- plotly::plot_ly(
    x = numeric(0), y = numeric(0),
    type = "scatter", mode = "markers", hoverinfo = "skip"
  )
  fig <- plotly::layout(
    fig,
    xaxis = list(visible = FALSE),
    yaxis = list(visible = FALSE),
    margin = list(l = 20, r = 20, t = 40, b = 20),
    showlegend = FALSE
  )
  if (!is.null(title) && nzchar(title)) {
    fig <- plotly::layout(fig, title = list(text = title, x = 0.5))
  }
  fig
}

.relationship_index <- function(study, relationship) {
  switch(relationship,
    "Express_normalized" = study$expression_index,
    "Activity_tf"        = study$activity_tf_index,
    "Activity_sig"       = study$activity_sig_index,
    NULL
  )
}

.gene_row_values <- function(study, gene, relationship) {
  gene_values(study, gene, relationship)
}

# --- Cluster plot -----------------------------------------------------------

.cluster_plot <- function(study, active_clusters, dot_size,
                          show_labels, cell_mask = NULL) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    return(.empty_plot())
  }
  mask <- .cells_in_active_clusters(study, active_clusters, cell_mask)
  cells <- study$cells[mask, , drop = FALSE]
  color_map <- .cluster_color_map(study$clusters)

  p <- plotly::plot_ly(
    cells,
    x = ~coord1, y = ~coord2,
    color = ~cellType, colors = color_map[unique(cells$cellType)],
    type = "scattergl", mode = "markers",
    marker = list(size = dot_size, opacity = 0.7),
    hovertemplate = paste0(
      "<b>%{customdata}</b><br>",
      study$meta$coordinate, "_1: %{x:.3f}<br>",
      study$meta$coordinate, "_2: %{y:.3f}<extra></extra>"
    ),
    customdata = ~cellID
  )

  if (isTRUE(show_labels) &&
      !is.null(study$clusters$label_1) &&
      !is.null(study$clusters$label_2)) {
    visible_labels <- study$clusters[
      study$clusters$cellType %in% (active_clusters %||% study$clusters$cellType), ,
      drop = FALSE
    ]
    if (nrow(visible_labels) > 0) {
      p <- plotly::add_annotations(
        p,
        x = visible_labels$label_1,
        y = visible_labels$label_2,
        text = visible_labels$cellType,
        showarrow = FALSE,
        font = list(size = 14, color = "#222"),
        bgcolor = "rgba(255,255,255,0.7)",
        bordercolor = "#666",
        borderpad = 3
      )
    }
  }

  plotly::layout(
    p,
    xaxis = list(title = paste0(study$meta$coordinate, "_1"),
                 zeroline = FALSE),
    yaxis = list(title = paste0(study$meta$coordinate, "_2"),
                 zeroline = FALSE,
                 scaleanchor = "x", scaleratio = 1),
    legend = list(title = list(text = "Cluster")),
    margin = list(l = 50, r = 20, t = 30, b = 50)
  )
}

# --- Feature plot -----------------------------------------------------------

.feature_plot <- function(study, gene, relationship, active_clusters,
                          dot_size, cell_mask = NULL) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    return(.empty_plot())
  }
  vals <- .gene_row_values(study, gene, relationship)
  if (is.null(vals)) {
    return(.empty_plot(paste0("No ", relationship, " data for ", gene)))
  }
  mask <- .cells_in_active_clusters(study, active_clusters, cell_mask)
  cells <- study$cells[mask, , drop = FALSE]
  v <- vals[mask]

  plotly::plot_ly(
    cells,
    x = ~coord1, y = ~coord2,
    type = "scattergl", mode = "markers",
    marker = list(
      size = dot_size,
      color = v,
      colorscale = list(c(0, "#dddddd"),
                         c(0.5, "#7c9fd1"),
                         c(1, "#1f3a72")),
      cmin = min(v, na.rm = TRUE),
      cmax = max(v, na.rm = TRUE),
      colorbar = list(title = list(text = gene)),
      opacity = 0.85
    ),
    hovertemplate = paste0(
      "<b>%{customdata}</b><br>",
      gene, ": %{marker.color:.3f}<extra></extra>"
    ),
    customdata = ~cellID
  ) |>
    plotly::layout(
      xaxis = list(title = paste0(study$meta$coordinate, "_1"),
                   zeroline = FALSE),
      yaxis = list(title = paste0(study$meta$coordinate, "_2"),
                   zeroline = FALSE,
                   scaleanchor = "x", scaleratio = 1),
      margin = list(l = 50, r = 20, t = 30, b = 50)
    )
}

# --- Violin plot ------------------------------------------------------------

.violin_plot <- function(study, gene, relationship, active_clusters,
                         cell_mask = NULL) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    return(.empty_plot())
  }
  vals <- .gene_row_values(study, gene, relationship)
  if (is.null(vals)) {
    return(.empty_plot(paste0("No ", relationship, " data for ", gene)))
  }
  mask <- .cells_in_active_clusters(study, active_clusters, cell_mask)
  df <- data.frame(
    cellType = study$cells$cellType[mask],
    value    = vals[mask],
    stringsAsFactors = FALSE
  )
  color_map <- .cluster_color_map(study$clusters)
  plotly::plot_ly(
    df,
    x = ~cellType, y = ~value, color = ~cellType,
    colors = color_map[unique(df$cellType)],
    type = "violin",
    box  = list(visible = TRUE),
    meanline = list(visible = TRUE),
    points = "outliers"
  ) |>
    plotly::layout(
      yaxis = list(title = gene),
      xaxis = list(title = ""),
      showlegend = FALSE,
      margin = list(l = 50, r = 20, t = 30, b = 80)
    )
}

# --- Heatmap (mean per cluster) --------------------------------------------

.cluster_aggregate_pair <- function(study, genes, relationship,
                                     active_clusters, cell_mask = NULL) {
  # Returns list(mean = matrix, pct = matrix), both gene x cluster.
  index <- .relationship_index(study, relationship)
  if (is.null(index) || length(genes) == 0) return(NULL)
  present <- genes[genes %in% index]
  if (length(present) == 0) return(NULL)
  ct <- study$cells$cellType
  sample_mask <- if (is.null(cell_mask)) rep(TRUE, length(ct))
                 else as.logical(cell_mask)
  clusters <- active_clusters %||% unique(ct)
  means <- matrix(0, nrow = length(present), ncol = length(clusters),
                  dimnames = list(present, clusters))
  pcts  <- matrix(0, nrow = length(present), ncol = length(clusters),
                  dimnames = list(present, clusters))
  for (i in seq_along(present)) {
    vals <- gene_values(study, present[i], relationship)
    if (is.null(vals)) next
    for (j in seq_along(clusters)) {
      mask <- (ct == clusters[j]) & sample_mask
      if (!any(mask)) next
      sub <- vals[mask]
      sub <- sub[!is.na(sub)]
      n <- length(sub)
      if (n == 0) next
      means[i, j] <- mean(sub)
      pcts[i, j]  <- sum(sub != 0) / n
    }
  }
  list(mean = means, pct = pcts)
}

.heatmap_plot <- function(study, genes, relationship, active_clusters,
                          cell_mask = NULL) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    return(.empty_plot())
  }
  agg <- .cluster_aggregate_pair(study, genes, relationship,
                                  active_clusters, cell_mask)
  if (is.null(agg)) {
    return(.empty_plot("Select gene(s) to build heatmap"))
  }
  m <- agg$mean
  # Blue-white-red diverging palette. zmid centres white at 0 so positive
  # values lean red and negative values (e.g. z-scored activity) lean blue.
  plotly::plot_ly(
    z = m,
    x = colnames(m), y = rownames(m),
    type = "heatmap",
    colorscale = list(c(0,   "#2166ac"),
                       c(0.5, "#ffffff"),
                       c(1,   "#b2182b")),
    zmid = 0,
    colorbar = list(title = list(text = "mean"))
  ) |>
    plotly::layout(
      xaxis = list(title = "Cluster"),
      yaxis = list(title = "Gene"),
      margin = list(l = 80, r = 20, t = 30, b = 80)
    )
}

# --- Bubble plot ------------------------------------------------------------

.bubble_plot <- function(study, genes, relationship, active_clusters,
                         cell_mask = NULL) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    return(.empty_plot())
  }
  agg <- .cluster_aggregate_pair(study, genes, relationship,
                                  active_clusters, cell_mask)
  if (is.null(agg)) {
    return(.empty_plot("Select gene(s) to build bubble plot"))
  }
  means <- agg$mean
  pcts  <- agg$pct
  rows <- rownames(means)
  cols <- colnames(means)
  grid <- expand.grid(gene = rows, cluster = cols,
                      stringsAsFactors = FALSE)
  grid$mean <- as.numeric(means[cbind(grid$gene, grid$cluster)])
  grid$pct  <- as.numeric(pcts[cbind(grid$gene, grid$cluster)])

  plotly::plot_ly(
    grid,
    x = ~cluster, y = ~gene,
    size = ~pct, sizes = c(2, 30),
    color = ~mean,
    colors = c("#f7f7f7", "#7c9fd1", "#1f3a72"),
    type = "scatter", mode = "markers",
    marker = list(sizemode = "area", line = list(width = 0)),
    hovertemplate = paste0(
      "<b>%{y}</b> in <b>%{x}</b><br>",
      "mean: %{marker.color:.3f}<br>",
      "pct expressing: %{marker.size:.1%}<extra></extra>"
    )
  ) |>
    plotly::layout(
      xaxis = list(title = "Cluster"),
      yaxis = list(title = "Gene", autorange = "reversed"),
      margin = list(l = 80, r = 20, t = 30, b = 80)
    )
}

# --- Network ----------------------------------------------------------------

.network_plot <- function(study, gene, network_type,
                          active_clusters, max_edges = 60) {
  edges_df <- switch(network_type,
    "TF"  = study$network_tf,
    "SIG" = study$network_sig,
    NULL
  )
  if (is.null(edges_df) || nrow(edges_df) == 0) return(NULL)
  edges_df <- edges_df[edges_df$source == gene | edges_df$target == gene, ,
                       drop = FALSE]
  if (!is.null(active_clusters) && length(active_clusters) > 0) {
    edges_df <- edges_df[edges_df$cellType %in% active_clusters, ,
                         drop = FALSE]
  }
  if (nrow(edges_df) == 0) return(NULL)
  edges_df <- edges_df[order(-abs(edges_df$mi)), , drop = FALSE]
  if (nrow(edges_df) > max_edges) {
    edges_df <- edges_df[seq_len(max_edges), , drop = FALSE]
  }

  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    return(edges_df)
  }
  nodes <- data.frame(
    id = unique(c(edges_df$source, edges_df$target)),
    stringsAsFactors = FALSE
  )
  nodes$label <- nodes$id
  nodes$color <- ifelse(nodes$id == gene, "#e15759", "#4e79a7")
  nodes$size  <- ifelse(nodes$id == gene, 32, 18)

  edges <- data.frame(
    from  = edges_df$source,
    to    = edges_df$target,
    value = pmax(abs(edges_df$mi), 1e-4),
    title = sprintf("mi=%.3f  pearson=%.3f  cellType=%s",
                    edges_df$mi, edges_df$pearson, edges_df$cellType),
    color = ifelse(edges_df$pearson >= 0, "#76b7b2", "#e15759"),
    stringsAsFactors = FALSE
  )
  visNetwork::visNetwork(nodes, edges, width = "100%", height = "640px") |>
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1)) |>
    visNetwork::visPhysics(stabilization = FALSE,
                            barnesHut = list(gravitationalConstant = -8000))
}
