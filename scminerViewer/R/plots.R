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

  # Direction encoding: positive correlation = activator (teal),
  # negative = repressor (red). We use spearman (rank-based, robust to
  # outliers — scMINER convention); fall back to pearson if spearman is
  # NA. Width = |MI| so dominant edges are visually heaviest.
  dir_score <- ifelse(is.na(edges_df$spearman),
                       edges_df$pearson, edges_df$spearman)
  direction <- ifelse(is.na(dir_score), "Unknown",
                       ifelse(dir_score >= 0, "Activator", "Repressor"))
  dir_color <- c(Activator = "#2E7D6A", Repressor = "#D7493A",
                  Unknown   = "#9e9e9e")

  # Per-neighbor dominant direction: the strongest |MI| edge to each
  # neighbor decides which side of the ring (and which color) it gets.
  # edges_df is already sorted by descending |MI|, so the first row that
  # mentions a neighbor carries that neighbor's dominant direction.
  edge_neighbor <- ifelse(edges_df$source == gene,
                           edges_df$target, edges_df$source)
  nb_seen <- !duplicated(edge_neighbor)
  nb_ids  <- edge_neighbor[nb_seen]
  nb_dir  <- direction[nb_seen]

  # Sort neighbors by (direction priority, descending MI). Activators
  # cluster on one arc, repressors on the opposite arc, with the
  # strongest of each group closest to the seam — reads cleanly at a
  # glance.
  dir_priority <- c(Activator = 0L, Unknown = 1L, Repressor = 2L)
  nb_mi <- abs(edges_df$mi[nb_seen])
  nb_order <- order(dir_priority[nb_dir], -nb_mi)
  nb_ids <- nb_ids[nb_order]

  nodes <- data.frame(
    id = c(gene, nb_ids),
    stringsAsFactors = FALSE
  )
  n_neighbors <- length(nb_ids)
  radius <- 260

  # Focus node: warm amber, slightly larger, white border ring.
  # Neighbor nodes: cool slate, smaller, white border.
  nodes$label <- nodes$id
  nodes$color <- c("#F2A33A", rep("#3F5775", n_neighbors))
  nodes$size  <- c(36L, rep(18L, n_neighbors))
  nodes$font.size  <- c(20L, rep(14L, n_neighbors))
  nodes$font.color <- "#1f1f1f"
  nodes$borderWidth <- c(3L, rep(1.5, n_neighbors))
  nodes$color.border <- "#ffffff"

  # Radial layout — focus at (0,0), neighbors evenly spaced starting at
  # 12 o'clock and walking clockwise so the sort above shows up visually.
  nodes$x <- 0
  nodes$y <- 0
  nodes$fixed <- c(TRUE, rep(FALSE, n_neighbors))
  if (n_neighbors >= 1L) {
    step <- 2 * pi / n_neighbors
    theta <- pi / 2 - (seq_len(n_neighbors) - 1L) * step
    nodes$x[-1L] <- radius * cos(theta)
    nodes$y[-1L] <- radius * sin(theta)
  }

  # Visual orientation: every edge is drawn as focus → neighbor so the
  # `to` arrow lands at the rim next to the non-focus gene regardless of
  # which side is `source` vs `target` in the underlying network. The
  # hover title preserves the original `source → target` so the true
  # regulation direction is still readable on hover.
  focus_is_source <- edges_df$source == gene
  vis_from <- ifelse(focus_is_source, edges_df$source, edges_df$target)
  vis_to   <- ifelse(focus_is_source, edges_df$target, edges_df$source)

  edges <- data.frame(
    from  = vis_from,
    to    = vis_to,
    value = pmax(abs(edges_df$mi), 1e-4),
    title = sprintf(
      paste0("<b>%s &rarr; %s</b><br>",
             "Direction: <b>%s</b><br>",
             "MI: %.3f<br>",
             "Spearman: %.3f<br>",
             "Pearson: %.3f<br>",
             "p-value: %.2g<br>",
             "Cell type: %s"),
      edges_df$source, edges_df$target,
      direction,
      edges_df$mi, edges_df$spearman, edges_df$pearson,
      edges_df$pvalue, edges_df$cellType
    ),
    color = unname(dir_color[direction]),
    arrows = "to",
    arrowStrikethrough = FALSE,
    stringsAsFactors = FALSE
  )

  ledges <- data.frame(
    label = c("Activator (+)", "Repressor (−)"),
    color = c(dir_color["Activator"], dir_color["Repressor"]),
    arrows = c("to", "to"),
    font.align = "top",
    stringsAsFactors = FALSE
  )

  main_html <- sprintf(
    paste0("<div style='text-align:center'>",
           "<span style='font-size:15px;font-weight:600;color:#1f1f1f'>",
           "%s network for %s</span>",
           "<br><span style='font-size:11px;color:#666'>",
           "%d edges &middot; width &prop; |MI| &middot; color = direction",
           "</span></div>"),
    network_type, gene, nrow(edges_df)
  )

  visNetwork::visNetwork(
      nodes, edges,
      width = "100%", height = "680px",
      main = list(text = main_html, style = "padding-top:6px;"),
      background = "#FBFBF9"
    ) |>
    visNetwork::visEdges(
      smooth = list(enabled = FALSE),       # straight lines so the arrow
                                             # sits cleanly at the target node
      # Width range 1.4 (smallest edge) -> 4 (largest edge). The min is
      # anchored so weak edges stay visible; the max is intentionally
      # restrained so heavy edges don't dominate the figure visually.
      scaling = list(min = 1.4, max = 4),
      # Small, restrained arrowheads — direction is clear from edge
      # orientation (source -> target) plus the per-edge color, so we
      # don't need a dominant arrow head. scaleFactor 0.7 keeps the
      # tip readable next to a node without overwhelming the plot.
      arrows = list(to = list(enabled = TRUE, scaleFactor = 0.7,
                               type = "arrow")),
      arrowStrikethrough = FALSE,            # arrow tip stops at the node border
      # Force the per-edge `color` column to apply to BOTH the line and
      # the arrow. Without `inherit = FALSE`, vis.js falls back to its
      # `inherit = "from"` default for arrow tinting, which can leave
      # the arrow a different shade than the line.
      color = list(inherit = FALSE, opacity = 0.78),
      hoverWidth = 0.5, selectionWidth = 1.2
    ) |>
    visNetwork::visNodes(
      shape = "dot",
      borderWidth = 2, borderWidthSelected = 3,
      shadow = list(enabled = TRUE, size = 10, x = 0, y = 2,
                     color = "rgba(0,0,0,0.15)"),
      font = list(face = "Helvetica, Arial, sans-serif")
    ) |>
    visNetwork::visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1,
                               algorithm = "hierarchical")
    ) |>
    # Physics disabled so the radial layout we computed for `nodes$x/y`
    # is exactly what gets drawn — focus gene stays pinned at center.
    visNetwork::visPhysics(enabled = FALSE) |>
    visNetwork::visInteraction(dragNodes = TRUE, dragView = TRUE,
                                zoomView = TRUE, tooltipDelay = 120,
                                hover = TRUE) |>
    visNetwork::visLegend(useGroups = FALSE, addEdges = ledges,
                           position = "right", main = "Direction",
                           width = 0.12, ncol = 1, zoom = FALSE)
}
