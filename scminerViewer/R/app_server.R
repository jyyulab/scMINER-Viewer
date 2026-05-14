#' Single-study server logic, callable with the standard Shiny signature.
#'
#' Returns a function suitable for `shiny::shinyApp(server = ...)`. The
#' browser app (`run_browser()`) reuses [`.app_server_logic()`] directly
#' so it can swap studies without rebuilding the whole app.
#'
#' @noRd
.app_server <- function(study) {
  function(input, output, session) {
    .app_server_logic(study, input, output, session)
  }
}

#' Body of the single-study server logic, decoupled from the closure.
#'
#' Splitting this out lets [run_browser()] call it for a study chosen
#' at session start.
#'
#' @noRd
.app_server_logic <- function(study, input, output, session) {
    if (!requireNamespace("shiny", quietly = TRUE)) return(invisible())

    # Server-side gene autocomplete keyed off the study's full gene list.
    # Default-load any genes the bundle marks as `defaults/genes`.
    defaults <- study$default_genes %||% character(0)
    shiny::updateSelectizeInput(
      session, "gene_select",
      choices = study$genes,
      selected = defaults,
      server = TRUE,
      options = list(maxOptions = 200,
                     placeholder = "Type to add gene(s)...")
    )

    # Reactive: which clusters are visible? Driven by DT row selection
    # (selected rows are the *visible* clusters).
    cluster_table_state <- shiny::reactiveValues(
      selected = study$clusters$cellType
    )

    # --- Downsampling ----------------------------------------------------
    # Shuffle each cluster's cell indices once per session; the reactive
    # mask then takes the first `pct%` from each shuffled vector. This is
    # deterministic for a given (cluster, pct) pair, so plots don't flicker
    # when the user switches tabs.
    .MAX_CELLS_DEFAULT <- 65000L
    set.seed(42L)
    sampling_order <- lapply(
      stats::setNames(study$clusters$cellType, study$clusters$cellType),
      function(ct) sample(which(study$cells$cellType == ct))
    )
    total_cells <- nrow(study$cells)
    initial_sampling_pct <- if (total_cells > .MAX_CELLS_DEFAULT) {
      max(1L, floor(.MAX_CELLS_DEFAULT / total_cells * 100))
    } else 100L
    if (initial_sampling_pct < 100L) {
      shiny::updateNumericInput(session, "sampling_percent",
                                value = initial_sampling_pct)
    }

    sampling_mask <- shiny::reactive({
      pct <- input$sampling_percent %||% 100
      pct <- max(1L, min(100L, as.integer(pct)))
      if (pct >= 100L) return(rep(TRUE, total_cells))
      out <- logical(total_cells)
      for (ord in sampling_order) {
        n_keep <- max(1L, round(length(ord) * pct / 100))
        out[ord[seq_len(n_keep)]] <- TRUE
      }
      out
    })

    output$cell_count <- shiny::renderText({
      mask <- .cells_in_active_clusters(study, cluster_table_state$selected,
                                         cell_mask = sampling_mask())
      sprintf("%s / %s",
              formatC(sum(mask),       big.mark = ",", format = "d"),
              formatC(nrow(study$cells), big.mark = ",", format = "d"))
    })

    # Per-cluster cell-index vectors. Precomputed once so the reactive
    # count only does a cheap mask-subset on each sampling change.
    cluster_indices <- lapply(
      as.character(study$clusters$cellType),
      function(ct) which(study$cells$cellType == ct)
    )
    names(cluster_indices) <- as.character(study$clusters$cellType)

    # Visible (sampled) count per cluster, reactive on the sampling mask.
    cluster_visible_counts <- shiny::reactive({
      mask <- sampling_mask()
      vapply(cluster_indices,
             function(idx) as.integer(sum(mask[idx])),
             integer(1))
    })

    # Build one <input type="checkbox"> that fires a Shiny event on click.
    # Cell type is read from a data-* attribute so we don't need to
    # escape the cellType into the JS string literal.
    .checkbox_html <- function(cell_type, checked) {
      safe <- gsub("[^A-Za-z0-9_]", "_", cell_type)
      shiny::HTML(sprintf(
        paste0(
          '<input type="checkbox" id="cluster_check_%s" ',
          'data-celltype="%s" %s ',
          'onclick="Shiny.setInputValue(',
          "'cluster_checkbox', ",
          '{ct: this.dataset.celltype, checked: this.checked, t: Date.now()}, ',
          "{priority: 'event'});\" />"
        ),
        safe,
        htmltools::htmlEscape(cell_type, attribute = TRUE),
        if (isTRUE(checked)) "checked" else ""
      ))
    }

    # Render the clusters as a plain HTML table via renderUI. This is
    # reactive on:
    #   - sampling_mask()             -> updates the Cells column
    #   - cluster_table_state$selected -> preserves checkbox `checked`
    #                                     state across re-renders
    output$clusters_table_ui <- shiny::renderUI({
      counts   <- cluster_visible_counts()
      selected <- cluster_table_state$selected
      totals   <- as.integer(study$clusters$count)
      colors   <- as.character(study$clusters$color)
      cts      <- as.character(study$clusters$cellType)

      rows <- lapply(seq_along(cts), function(i) {
        ct <- cts[i]
        cells_text <- sprintf(
          "%s / %s",
          formatC(counts[[ct]], big.mark = ",", format = "d"),
          formatC(totals[i],    big.mark = ",", format = "d")
        )
        swatch <- shiny::HTML(sprintf(
          paste0("<span style='display:inline-block;width:24px;",
                 "height:14px;background:%s;border:1px solid #ccc;'></span>"),
          colors[i]
        ))
        shiny::tags$tr(
          shiny::tags$td(.checkbox_html(ct, ct %in% selected),
                          class = "ct-show"),
          shiny::tags$td(ct, class = "ct-name"),
          shiny::tags$td(cells_text, class = "ct-count"),
          shiny::tags$td(swatch, class = "ct-color")
        )
      })

      shiny::tags$table(
        class = "cluster-table",
        shiny::tags$thead(shiny::tags$tr(
          shiny::tags$th("Show"),
          shiny::tags$th("Cluster"),
          shiny::tags$th("Cells"),
          shiny::tags$th("Color")
        )),
        shiny::tags$tbody(rows)
      )
    })

    # Listen for checkbox toggles fired by the inline onclick handler.
    shiny::observeEvent(input$cluster_checkbox, ignoreNULL = TRUE, {
      msg <- input$cluster_checkbox
      if (is.null(msg) || is.null(msg$ct)) return()
      cur <- cluster_table_state$selected
      if (isTRUE(msg$checked)) {
        cluster_table_state$selected <- unique(c(cur, msg$ct))
      } else {
        cluster_table_state$selected <- setdiff(cur, msg$ct)
      }
    })

    active_clusters <- shiny::reactive(cluster_table_state$selected)
    selected_genes  <- shiny::reactive(input$gene_select %||% character(0))

    # ---- Cluster plot ---------------------------------------------------
    output$cluster_plot <- plotly::renderPlotly({
      .cluster_plot(
        study           = study,
        active_clusters = active_clusters(),
        dot_size        = input$dot_size %||% 4,
        show_labels     = isTRUE(input$show_labels),
        cell_mask       = sampling_mask()
      )
    })

    # ---- Heatmap --------------------------------------------------------
    output$heatmap_plot <- plotly::renderPlotly({
      .heatmap_plot(
        study           = study,
        genes           = selected_genes(),
        relationship    = "Express_normalized",
        active_clusters = active_clusters(),
        cell_mask       = sampling_mask()
      )
    })

    # ---- Bubble ---------------------------------------------------------
    output$bubble_plot <- plotly::renderPlotly({
      .bubble_plot(
        study           = study,
        genes           = selected_genes(),
        relationship    = "Express_normalized",
        active_clusters = active_clusters(),
        cell_mask       = sampling_mask()
      )
    })

    # ---- Per-gene 3-level nested tabs (Gene -> CellType -> Relationship)
    # The innermost level is a `navset_card_pill` "slide card" toggling
    # between Expression / TF / SIG (or just TF / SIG for the Network
    # panel, which has no expression network).
    cell_types_all <- as.character(study$clusters$cellType)
    rel_labels_value <- c("Expression", "TF", "SIG")
    rel_labels_network <- c("TF", "SIG")
    .rel_key <- function(rel) {
      switch(rel,
        Expression = "Express_normalized",
        TF         = "Activity_tf",
        SIG        = "Activity_sig",
        rel
      )
    }
    .out_id <- function(plot_kind, gene, ct, rel) {
      sprintf("%s_%s__%s__%s", plot_kind,
              make.names(gene), make.names(ct), make.names(rel))
    }

    .three_level_tabs <- function(genes, plot_kind) {
      if (length(genes) == 0) {
        return(shiny::div(class = "no-data-msg",
                          "Add gene(s) to view this plot."))
      }
      cell_options <- c("All", cell_types_all)
      rels <- if (plot_kind == "network") rel_labels_network
              else rel_labels_value
      gene_panels <- lapply(genes, function(g) {
        ct_panels <- lapply(cell_options, function(ct) {
          rel_panels <- lapply(rels, function(rel) {
            oid <- .out_id(plot_kind, g, ct, rel)
            inner <- if (plot_kind == "network") {
              visNetwork::visNetworkOutput(oid, height = "620px")
            } else {
              plotly::plotlyOutput(oid, height = "620px")
            }
            bslib::nav_panel(rel, inner)
          })
          ct_id <- sprintf("%s_rel_%s_%s", plot_kind,
                           make.names(g), make.names(ct))
          bslib::nav_panel(
            ct,
            do.call(bslib::navset_card_pill,
                    c(list(id = ct_id), rel_panels))
          )
        })
        g_id <- sprintf("%s_ct_%s", plot_kind, make.names(g))
        bslib::nav_panel(
          g,
          do.call(bslib::navset_pill,
                  c(list(id = g_id), ct_panels))
        )
      })
      do.call(bslib::navset_tab,
              c(list(id = paste0(plot_kind, "_tabs")), gene_panels))
    }

    output$feature_tabs_ui <- shiny::renderUI({
      .three_level_tabs(selected_genes(), "feature")
    })
    output$violin_tabs_ui <- shiny::renderUI({
      .three_level_tabs(selected_genes(), "violin")
    })
    output$network_tabs_ui <- shiny::renderUI({
      .three_level_tabs(selected_genes(), "network")
    })

    # Effective clusters for a (global_selection, tab_ct) pair: "All" =
    # respect the global selection; specific cell type = only that one,
    # but empty if globally hidden.
    .effective_clusters <- function(global, ct) {
      if (identical(ct, "All")) return(global)
      if (ct %in% global) ct else character(0)
    }

    # Dynamically register an output for every (gene, cell_type, rel)
    # combination as the gene list changes. Shiny only evaluates
    # outputs that are actually visible, so even if the combinatorial
    # count is large, only the focused tab does work.
    shiny::observe({
      genes <- selected_genes()
      cell_options <- c("All", cell_types_all)
      for (g in genes) for (ct in cell_options) {
        for (rel in rel_labels_value) local({
          g_l <- g; ct_l <- ct; rel_l <- rel
          rel_k <- .rel_key(rel_l)
          output[[.out_id("feature", g_l, ct_l, rel_l)]] <-
            plotly::renderPlotly({
              eff <- .effective_clusters(active_clusters(), ct_l)
              .feature_plot(study = study, gene = g_l,
                            relationship = rel_k,
                            active_clusters = eff,
                            dot_size = input$dot_size %||% 4,
                            cell_mask = sampling_mask())
            })
          output[[.out_id("violin", g_l, ct_l, rel_l)]] <-
            plotly::renderPlotly({
              eff <- .effective_clusters(active_clusters(), ct_l)
              .violin_plot(study = study, gene = g_l,
                           relationship = rel_k,
                           active_clusters = eff,
                           cell_mask = sampling_mask())
            })
        })
        for (rel in rel_labels_network) local({
          g_l <- g; ct_l <- ct; rel_l <- rel
          output[[.out_id("network", g_l, ct_l, rel_l)]] <-
            visNetwork::renderVisNetwork({
              eff <- .effective_clusters(active_clusters(), ct_l)
              plt <- .network_plot(study = study, gene = g_l,
                                    network_type = rel_l,
                                    active_clusters = eff)
              if (is.null(plt) || is.data.frame(plt)) {
                visNetwork::visNetwork(
                  data.frame(id = integer(0)),
                  data.frame(from = integer(0), to = integer(0))
                )
              } else plt
            })
        })
      }
    })
}
