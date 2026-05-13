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

.cells_in_active_clusters <- function(study, active_clusters) {
  if (is.null(active_clusters) || length(active_clusters) == 0) {
    return(rep(TRUE, nrow(study$cells)))
  }
  study$cells$cellType %in% active_clusters
}

.relationship_matrix <- function(study, relationship) {
  switch(relationship,
    "Express_normalized" = study$expression,
    "Activity_tf"        = study$activity_tf,
    "Activity_sig"       = study$activity_sig,
    NULL
  )
}

.gene_row_values <- function(study, gene, relationship) {
  mat <- .relationship_matrix(study, relationship)
  if (is.null(mat)) return(NULL)
  row_idx <- match(gene, study$genes)
  if (is.na(row_idx)) return(NULL)
  as.numeric(mat[row_idx, ])
}

# --- Cluster plot -----------------------------------------------------------

.cluster_plot <- function(study, active_clusters, dot_size,
                          show_labels) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    return(plotly::plotly_empty())
  }
  mask <- .cells_in_active_clusters(study, active_clusters)
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
                          dot_size) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    return(plotly::plotly_empty())
  }
  vals <- .gene_row_values(study, gene, relationship)
  if (is.null(vals)) {
    return(plotly::plotly_empty() %>%
             plotly::layout(title = paste0("No ", relationship,
                                           " data for ", gene)))
  }
  mask <- .cells_in_active_clusters(study, active_clusters)
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

.violin_plot <- function(study, gene, relationship, active_clusters) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    return(plotly::plotly_empty())
  }
  vals <- .gene_row_values(study, gene, relationship)
  if (is.null(vals)) {
    return(plotly::plotly_empty() %>%
             plotly::layout(title = paste0("No ", relationship,
                                           " data for ", gene)))
  }
  mask <- .cells_in_active_clusters(study, active_clusters)
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

.cluster_aggregate <- function(study, genes, relationship, active_clusters,
                                fun = mean) {
  mat <- .relationship_matrix(study, relationship)
  if (is.null(mat) || length(genes) == 0) return(NULL)
  rows <- match(genes, study$genes)
  ok <- !is.na(rows)
  if (!any(ok)) return(NULL)
  rows <- rows[ok]
  genes <- genes[ok]
  ct <- study$cells$cellType
  clusters <- active_clusters %||% unique(ct)
  out <- matrix(0, nrow = length(rows), ncol = length(clusters),
                dimnames = list(genes, clusters))
  for (j in seq_along(clusters)) {
    cell_mask <- ct == clusters[j]
    if (!any(cell_mask)) next
    sub <- mat[rows, cell_mask, drop = FALSE]
    out[, j] <- apply(sub, 1, fun)
  }
  out
}

.pct_expressing <- function(study, genes, relationship, active_clusters) {
  mat <- .relationship_matrix(study, relationship)
  if (is.null(mat) || length(genes) == 0) return(NULL)
  rows <- match(genes, study$genes)
  ok <- !is.na(rows)
  if (!any(ok)) return(NULL)
  rows <- rows[ok]
  genes <- genes[ok]
  ct <- study$cells$cellType
  clusters <- active_clusters %||% unique(ct)
  out <- matrix(0, nrow = length(rows), ncol = length(clusters),
                dimnames = list(genes, clusters))
  for (j in seq_along(clusters)) {
    cell_mask <- ct == clusters[j]
    n <- sum(cell_mask)
    if (n == 0) next
    sub <- mat[rows, cell_mask, drop = FALSE]
    out[, j] <- Matrix::rowSums(sub != 0) / n
  }
  out
}

.heatmap_plot <- function(study, genes, relationship, active_clusters) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    return(plotly::plotly_empty())
  }
  m <- .cluster_aggregate(study, genes, relationship, active_clusters)
  if (is.null(m) || nrow(m) == 0) {
    return(plotly::plotly_empty() %>%
             plotly::layout(title = "Select gene(s) to build heatmap"))
  }
  plotly::plot_ly(
    z = m,
    x = colnames(m), y = rownames(m),
    type = "heatmap",
    colorscale = list(c(0, "#f7f7f7"),
                       c(0.5, "#7c9fd1"),
                       c(1, "#1f3a72")),
    colorbar = list(title = list(text = "mean"))
  ) |>
    plotly::layout(
      xaxis = list(title = "Cluster"),
      yaxis = list(title = "Gene"),
      margin = list(l = 80, r = 20, t = 30, b = 80)
    )
}

# --- Bubble plot ------------------------------------------------------------

.bubble_plot <- function(study, genes, relationship, active_clusters) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    return(plotly::plotly_empty())
  }
  means <- .cluster_aggregate(study, genes, relationship, active_clusters)
  pcts  <- .pct_expressing(study, genes, relationship, active_clusters)
  if (is.null(means) || is.null(pcts)) {
    return(plotly::plotly_empty() %>%
             plotly::layout(title = "Select gene(s) to build bubble plot"))
  }
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
