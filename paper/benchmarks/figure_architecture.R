#!/usr/bin/env Rscript
# paper/benchmarks/figure_architecture.R
#
# Renders the scMINER Viewer architecture as a single vector figure.
# Drives a base-graphics layout with rect()/arrows()/text() -- no
# extra package deps beyond what figures.R already loads.
#
# Writes:
#   paper/figures/architecture.pdf
#   paper/figures/architecture.png
#
# Run:
#   Rscript paper/benchmarks/figure_architecture.R

dir.create("paper/figures", recursive = TRUE, showWarnings = FALSE)

# --- Theme (NPG-aligned palette, matched to figure1) -----------------------

col_input    <- "#3C548822"   # very pale steel-blue (input panel)
col_input_ed <- "#3C5488"
col_prep     <- "#E64B3522"   # pale red (prepare_study)
col_prep_ed  <- "#E64B35"
col_artifact <- "#00A08722"   # pale teal (shared artifact)
col_artif_ed <- "#00A087"
col_viewer_r <- "#4DBBD522"   # pale blue (R viewer)
col_viewer_p <- "#F39B7F22"   # pale orange (Python viewer)
col_viewer_e <- "#3C5488"
col_webui    <- "#8491B422"   # pale slate (webui)
col_webui_ed <- "#8491B4"
col_arrow    <- "#444444"

draw_box <- function(x0, y0, x1, y1, fill, border,
                      title, body = NULL,
                      title_cex = 1.0, body_cex = 0.8,
                      title_font = 2, body_font = 1,
                      title_col = "black", body_col = "#222222",
                      pad_top = 0.4, line_gap = 0.5,
                      body_adj = c(0.5, 1)) {
  rect(x0, y0, x1, y1, col = fill, border = border, lwd = 1.4)
  cx <- (x0 + x1) / 2
  text(cx, y1 - pad_top, title,
       cex = title_cex, font = title_font, col = title_col)
  if (!is.null(body)) {
    # Accept either a vector of pre-split lines OR a single string
    # carrying embedded newlines. Both flatten to one char vector.
    lines <- unlist(strsplit(as.character(body), "\n", fixed = TRUE),
                    use.names = FALSE)
    y_start <- y1 - pad_top - 0.65
    for (k in seq_along(lines)) {
      text(cx, y_start - (k - 1) * line_gap,
           lines[k], cex = body_cex, font = body_font, col = body_col,
           adj = body_adj)
    }
  }
}

draw_arrow <- function(x0, y0, x1, y1, lwd = 1.6, length = 0.10) {
  arrows(x0, y0, x1, y1, lwd = lwd, length = length, col = col_arrow,
         angle = 25)
}

# --- Layout coordinates (12 wide x 20 tall) ---------------------------------

W <- 12; H <- 20
render_arch <- function() {
  par(mar = c(0.4, 0.4, 0.4, 0.4))
  plot.new()
  plot.window(xlim = c(0, W), ylim = c(0, H), asp = NA)

  # --- Tier 1: Inputs (top row) -----------------------------------------
  in_top <- H - 0.4
  in_bot <- H - 3.0
  inputs <- list(
    list(t = "expression.rds",
         body = "Biobase ExpressionSet\n(exprs + pData + fData)"),
    list(t = "activity.rds",
         body = "ExpressionSet\nrows _TF / _SIG"),
    list(t = "networks.txt",
         body = "SJARACNe TSV\nsource, target, MI..."),
    list(t = "config.yaml",
         body = "column names, paths,\npalette")
  )
  n_inputs <- length(inputs)
  pad <- 0.25
  total_w <- W - 2 * pad
  box_w <- total_w / n_inputs - 0.20
  x_left <- pad
  for (i in seq_along(inputs)) {
    x0 <- x_left + (i - 1) * (box_w + 0.20)
    x1 <- x0 + box_w
    draw_box(x0, in_bot, x1, in_top,
             fill = col_input, border = col_input_ed,
             title = inputs[[i]]$t,
             body  = inputs[[i]]$body,
             title_cex = 0.95, body_cex = 0.75,
             pad_top = 0.42, line_gap = 0.46)
  }
  text(W / 2, in_top + 0.15, "INPUTS  (per study)",
       cex = 0.78, font = 3, col = col_input_ed, adj = c(0.5, 0))

  # Arrows from each input box down to the prepare_study box.
  prep_top <- in_bot - 1.5
  prep_bot <- prep_top - 2.6
  prep_x0  <- pad + 0.6
  prep_x1  <- W - pad - 0.6
  prep_xc  <- (prep_x0 + prep_x1) / 2
  for (i in seq_along(inputs)) {
    x0 <- x_left + (i - 1) * (box_w + 0.20) + box_w / 2
    draw_arrow(x0, in_bot - 0.05, prep_xc + (x0 - W / 2) * 0.2, prep_top + 0.05)
  }

  # --- Tier 2: prepare_study() -----------------------------------------
  draw_box(prep_x0, prep_bot, prep_x1, prep_top,
           fill = col_prep, border = col_prep_ed,
           title = "scminerViewer::prepare_study(config_path)",
           body = c(
             "1.  load_study_config()        YAML + defaults",
             "2.  extract_{cells, genes, expression, activity, networks}()",
             "3.  .write_graph_shards()      row-by-row gzip CSVs",
             "4.  write_bundle()             HDF5 (no matrix values)"
           ),
           title_cex = 1.00, body_cex = 0.80,
           pad_top = 0.45, line_gap = 0.48)

  # Arrow down to shared artifact
  art_top <- prep_bot - 1.0
  art_bot <- art_top - 3.1
  art_xc  <- W / 2
  draw_arrow(prep_xc, prep_bot - 0.05, art_xc, art_top + 0.05)

  # --- Tier 3: Shared lazy artifact -----------------------------------------
  art_x0 <- pad + 0.3
  art_x1 <- W - pad - 0.3
  draw_box(art_x0, art_bot, art_x1, art_top,
           fill = col_artifact, border = col_artif_ed,
           title = "SHARED LAZY ARTIFACT  (bundleVersion = 1)",
           body = c(
             "<output>/<studyID>/",
             "    <studyID>.scminer.h5            metadata + indexes",
             "    expression_files/<sid>/<letter>/<gene>.csv.gz",
             "    activity_files/<sid>/{TF,SIG}/<letter>/<gene>.csv.gz"
           ),
           title_cex = 1.00, body_cex = 0.78,
           pad_top = 0.45, line_gap = 0.48)

  # Arrows fanning out to the two viewer boxes
  view_top <- art_bot - 1.2
  view_bot <- view_top - 3.2
  vR_x0 <- pad + 0.6
  vR_x1 <- W / 2 - 0.4
  vP_x0 <- W / 2 + 0.4
  vP_x1 <- W - pad - 0.6
  draw_arrow(art_xc - 1.2, art_bot - 0.05, (vR_x0 + vR_x1) / 2, view_top + 0.05)
  draw_arrow(art_xc + 1.2, art_bot - 0.05, (vP_x0 + vP_x1) / 2, view_top + 0.05)

  # --- Tier 4: R / Python viewers -----------------------------------------
  draw_box(vR_x0, view_bot, vR_x1, view_top,
           fill = col_viewer_r, border = col_viewer_e,
           title = "R viewer  (scminerViewer)",
           body = c(
             "load_study()",
             "gene_values()",
             "run_app()",
             "run_browser()"
           ),
           title_cex = 0.95, body_cex = 0.80,
           pad_top = 0.45, line_gap = 0.48)
  draw_box(vP_x0, view_bot, vP_x1, view_top,
           fill = col_viewer_p, border = col_viewer_e,
           title = "Python viewer  (scminer_viewer)",
           body = c(
             "load_study()",
             "gene_values()",
             "run_app()",
             "scminer-viewer browse <root>"
           ),
           title_cex = 0.95, body_cex = 0.80,
           pad_top = 0.45, line_gap = 0.48)

  # --- Tier 5: Shiny webui ------------------------------------------------
  web_top <- view_bot - 1.2
  web_bot <- web_top - 2.2
  web_x0 <- W / 2 - 3.8
  web_x1 <- W / 2 + 3.8
  draw_arrow((vR_x0 + vR_x1) / 2, view_bot - 0.05, W / 2 - 0.5, web_top + 0.05)
  draw_arrow((vP_x0 + vP_x1) / 2, view_bot - 0.05, W / 2 + 0.5, web_top + 0.05)
  draw_box(web_x0, web_bot, web_x1, web_top,
           fill = col_webui, border = col_webui_ed,
           title = "Shiny / Shiny-for-Python webui",
           body = c(
             "3 nested tabs:",
             "Gene  ->  CellType  ->  { Expression, TF, SIG }"
           ),
           title_cex = 0.95, body_cex = 0.85,
           pad_top = 0.55, line_gap = 0.60)

  # Footer note (sits below the webui box with breathing room)
  text(W / 2, max(0.25, web_bot - 0.55),
       "Single writer.  Twin readers.  216 tests (194 R + 22 Python) lock the contract.",
       cex = 0.78, font = 3, col = "#555555")
}

# --- Write PDF + PNG -------------------------------------------------------

w_in <- 7.5
h_in <- 9.0

pdf("paper/figures/architecture.pdf", width = w_in, height = h_in)
render_arch()
dev.off()

png("paper/figures/architecture.png", width = w_in, height = h_in,
    units = "in", res = 300)
render_arch()
dev.off()

message("Wrote paper/figures/architecture.pdf")
message("Wrote paper/figures/architecture.png")
