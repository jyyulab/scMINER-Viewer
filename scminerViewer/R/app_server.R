.app_server <- function(study) {
  function(input, output, session) {
    if (!requireNamespace("shiny", quietly = TRUE)) return(invisible())

    # Server-side gene autocomplete keyed off the study's full gene list.
    shiny::updateSelectizeInput(
      session, "gene_select",
      choices = study$genes,
      selected = character(0),
      server = TRUE,
      options = list(maxOptions = 200,
                     placeholder = "Type to add gene(s)...")
    )

    # Reactive: which clusters are visible? Driven by DT row selection
    # (selected rows are the *visible* clusters).
    cluster_table_state <- shiny::reactiveValues(
      selected = study$clusters$cellType
    )

    output$cell_count <- shiny::renderText({
      mask <- .cells_in_active_clusters(study, cluster_table_state$selected)
      sprintf("%s / %s",
              formatC(sum(mask),       big.mark = ",", format = "d"),
              formatC(nrow(study$cells), big.mark = ",", format = "d"))
    })

    output$clusters_table <- DT::renderDT({
      df <- data.frame(
        Cluster = study$clusters$cellType,
        Cells   = formatC(study$clusters$count, big.mark = ",", format = "d"),
        Color   = sprintf(
          "<span style='display:inline-block;width:24px;height:14px;background:%s;border:1px solid #ccc;'></span>",
          study$clusters$color
        ),
        stringsAsFactors = FALSE
      )
      DT::datatable(
        df,
        escape    = FALSE,
        rownames  = FALSE,
        selection = list(mode = "multiple", selected = seq_len(nrow(df))),
        options = list(
          dom = "t", paging = FALSE, searching = FALSE, info = FALSE,
          ordering = TRUE
        )
      )
    })

    shiny::observeEvent(input$clusters_table_rows_selected, ignoreNULL = FALSE, {
      rows <- input$clusters_table_rows_selected
      if (is.null(rows) || length(rows) == 0) {
        cluster_table_state$selected <- character(0)
      } else {
        cluster_table_state$selected <- study$clusters$cellType[rows]
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
        show_labels     = isTRUE(input$show_labels)
      )
    })

    # ---- Heatmap --------------------------------------------------------
    output$heatmap_plot <- plotly::renderPlotly({
      .heatmap_plot(
        study           = study,
        genes           = selected_genes(),
        relationship    = "Express_normalized",
        active_clusters = active_clusters()
      )
    })

    # ---- Bubble ---------------------------------------------------------
    output$bubble_plot <- plotly::renderPlotly({
      .bubble_plot(
        study           = study,
        genes           = selected_genes(),
        relationship    = "Express_normalized",
        active_clusters = active_clusters()
      )
    })

    # ---- Per-gene nested tabs (Feature / Violin / Network) --------------
    .nested_tabs <- function(genes, panel_fn, ns_prefix) {
      if (length(genes) == 0) {
        return(shiny::div(class = "no-data-msg",
                          "Add gene(s) to view this plot."))
      }
      tabs <- lapply(genes, function(g) {
        bslib::nav_panel(
          g,
          panel_fn(g)
        )
      })
      do.call(bslib::navset_tab, c(list(id = ns_prefix), tabs))
    }

    output$feature_tabs_ui <- shiny::renderUI({
      genes <- selected_genes()
      .nested_tabs(genes,
        panel_fn = function(g) {
          plotly::plotlyOutput(paste0("feature_", make.names(g)),
                               height = "640px")
        },
        ns_prefix = "feature_tabs"
      )
    })

    output$violin_tabs_ui <- shiny::renderUI({
      genes <- selected_genes()
      .nested_tabs(genes,
        panel_fn = function(g) {
          plotly::plotlyOutput(paste0("violin_", make.names(g)),
                               height = "640px")
        },
        ns_prefix = "violin_tabs"
      )
    })

    output$network_tabs_ui <- shiny::renderUI({
      genes <- selected_genes()
      .nested_tabs(genes,
        panel_fn = function(g) {
          shiny::tagList(
            shiny::radioButtons(paste0("netkind_", make.names(g)),
                                label = NULL,
                                choices = c("TF", "SIG"),
                                selected = "TF", inline = TRUE),
            visNetwork::visNetworkOutput(paste0("network_", make.names(g)),
                                          height = "640px")
          )
        },
        ns_prefix = "network_tabs"
      )
    })

    # Render per-gene plots dynamically. Each iteration generates an
    # observer that wires output$feature_<g> / violin_<g> / network_<g>.
    shiny::observe({
      genes <- selected_genes()
      for (g in genes) local({
        gene_local <- g
        oid <- make.names(gene_local)
        output[[paste0("feature_", oid)]] <- plotly::renderPlotly({
          .feature_plot(
            study           = study,
            gene            = gene_local,
            relationship    = "Express_normalized",
            active_clusters = active_clusters(),
            dot_size        = input$dot_size %||% 4
          )
        })
        output[[paste0("violin_", oid)]] <- plotly::renderPlotly({
          .violin_plot(
            study           = study,
            gene            = gene_local,
            relationship    = "Express_normalized",
            active_clusters = active_clusters()
          )
        })
        output[[paste0("network_", oid)]] <- visNetwork::renderVisNetwork({
          kind <- input[[paste0("netkind_", oid)]] %||% "TF"
          plt <- .network_plot(
            study           = study,
            gene            = gene_local,
            network_type    = kind,
            active_clusters = active_clusters()
          )
          if (is.null(plt) || is.data.frame(plt)) {
            visNetwork::visNetwork(
              data.frame(id = integer(0)),
              data.frame(from = integer(0), to = integer(0))
            )
          } else {
            plt
          }
        })
      })
    })
  }
}
