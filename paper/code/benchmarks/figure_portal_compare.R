#!/usr/bin/env Rscript
# paper/code/benchmarks/figure_portal_compare.R
#
# Render two figure sets from the joined wide table produced by
# paper/code/portal/portal_compare.R (paper/metrics/portal_studies/portal_studies_compare.tsv).
# Figure 3 (the full-mode 6-panel figure for the manuscript) is rendered
# by paper/code/benchmarks/figure_portal.R from the same compare data; this
# script handles the supplemental and paired-comparison views only:
#
#   1. Supplemental expression-only: one row per study with a valid
#      expression-only run. Shows what the "no TF/sig, no activity"
#      baseline looks like on its own (sits in the supplement).
#   2. Paired compare:               the TF/sig-eligible studies, paired
#      full vs expression-only side by side (geom_segment + two
#      geom_point markers) so the per-study cost of TF/sig is visible.
#
# Outputs (under paper/figures/compare/):
#   figureS_expr_only_{A..D}.{pdf,png}  four standalone supplemental panels
#   figureS_expr_only.{pdf,png}         combined supplemental figure
#   compare_{A..D}.{pdf,png}            four standalone paired panels
#   compare.{pdf,png}                   combined paired comparison
#
# Run from the project root (after the compare TSV exists):
#   Rscript paper/code/benchmarks/figure_portal_compare.R
#
# CLI overrides:
#   --compare-tsv <path>  default paper/metrics/portal_studies/portal_studies_compare.tsv
#   --figures-dir <dir>   default paper/figures/compare

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

.cli <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  i <- which(args == name)
  if (length(i) == 0L || i == length(args)) return(default)
  args[i + 1L]
}
compare_tsv <- .cli("--compare-tsv",
                     "paper/metrics/portal_studies/portal_studies_compare.tsv")
figures_dir <- .cli("--figures-dir", "paper/figures/compare")

# NPG palette pulled from ggsci so we don't pick up a dependency:
#   E64B35 red, 4DBBD5 cyan, 00A087 green, 3C5488 navy, F39B7F peach,
#   91D1C2 mint, 8491B4 slate.
npg <- c(red    = "#E64B35", cyan = "#4DBBD5", green = "#00A087",
         navy   = "#3C5488", peach = "#F39B7F",
         mint   = "#91D1C2", slate = "#8491B4")

theme_paper <- function() {
  theme_classic(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 11,
                                       margin = margin(b = 4)),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position  = "top",
      legend.key.size  = unit(10, "pt"),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 9),
      axis.text        = element_text(size = 9),
      axis.title       = element_text(size = 10)
    )
}

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

.cairo_ok <- isTRUE(capabilities("cairo")) &&
  is.null(tryCatch(
    withCallingHandlers(
      {
        grDevices::cairo_pdf(tempfile(fileext = ".pdf"),
                              width = 1, height = 1)
        try(grDevices::dev.off(), silent = TRUE)
        NULL
      },
      warning = function(w) {
        if (grepl("cairo", conditionMessage(w), ignore.case = TRUE)) {
          stop("cairo unavailable")
        }
        invokeRestart("muffleWarning")
      }),
    error = function(e) e))
pdf_device <- if (.cairo_ok) cairo_pdf else grDevices::pdf
if (!.cairo_ok) {
  message("(cairo PDF device unavailable; using default pdf() instead)")
}

save_panel <- function(plot, name, width, height) {
  pdf_path <- file.path(figures_dir, paste0(name, ".pdf"))
  png_path <- file.path(figures_dir, paste0(name, ".png"))
  ggsave(pdf_path, plot, width = width, height = height,
         units = "in", device = pdf_device)
  ggsave(png_path, plot, width = width, height = height,
         units = "in", dpi = 300)
  message("  ", pdf_path, " + ", basename(png_path))
}

# ---- 1. Load + slice -------------------------------------------------------

if (!file.exists(compare_tsv)) {
  stop("Compare TSV not found: ", compare_tsv,
       "\nRun paper/code/portal/portal_compare.R first.")
}
df <- utils::read.table(compare_tsv, sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE,
                         check.names = FALSE,
                         na.strings = c("", "NA"))
message(sprintf("Loaded %d rows from %s", nrow(df), compare_tsv))

# Coalesce shape columns: prefer the full-mode value, fall back to
# expression-only when the full row is missing (e.g. 2322 ran only
# expression-only). Without this, those studies have NA on the x-axis
# and silently disappear from the plot.
.coalesce <- function(d, base) {
  a <- paste0(base, "_full")
  b <- paste0(base, "_expr_only")
  full <- if (a %in% names(d)) as.numeric(d[[a]]) else rep(NA_real_, nrow(d))
  expr <- if (b %in% names(d)) as.numeric(d[[b]]) else rep(NA_real_, nrow(d))
  ifelse(is.na(full), expr, full)
}
df$n_cells_any <- .coalesce(df, "n_cells")
df$n_genes_any <- .coalesce(df, "n_genes")

# Studies whose full-mode run actually carried TF/sig edges. A few studies
# in the manifest have no networks file on disk; their full-mode run is
# byte-identical to expression-only and would confuse the comparison
# plots, so we drop them per the user's instruction.
df$has_tfsig <- ((!is.na(df$net_tf_edges_full)  & df$net_tf_edges_full  > 0L) |
                 (!is.na(df$net_sig_edges_full) & df$net_sig_edges_full > 0L))

# Helper: a row is valid for the supplemental expression-only figure
# set if BOTH modes produced metrics. Studies where the full-mode run
# is missing (e.g. 2322) are excluded so the supplemental figure only
# plots studies whose comparison benchmark ran cleanly end to end.
expr_ok <- !is.na(df$status_expr_only) & df$status_expr_only == "ok" &
  !is.na(df$prepare_seconds_expr_only) &
  !is.na(df$status_full) & df$status_full == "ok"
df_expr <- df[expr_ok, , drop = FALSE]

# A row is valid for the paired "compare" figure set if BOTH runs
# produced metrics AND TF/sig was actually present in the full run.
full_ok <- !is.na(df$status_full) & df$status_full == "ok" &
  !is.na(df$prepare_seconds_full) & df$has_tfsig
both_ok <- expr_ok & full_ok
df_cmp <- df[both_ok, , drop = FALSE]

message(sprintf("  expression-only valid rows: %d", nrow(df_expr)))
message(sprintf("  paired (both modes) rows:   %d", nrow(df_cmp)))

# Derive plot-friendly columns up front so each panel stays one-liner-ish.
add_derived <- function(d, suffix) {
  bb <- d[[paste0("bundle_bytes_", suffix)]]
  d[[paste0("bundle_mb_", suffix)]] <- as.numeric(bb) / 1024^2
  d
}
df_expr <- add_derived(df_expr, "expr_only")
df_cmp  <- add_derived(df_cmp,  "expr_only")
df_cmp  <- add_derived(df_cmp,  "full")

# Reused axis aesthetics
sx_log_cells <- list(
  scale_x_log10(labels = label_log()),
  labs(x = "Cells (log10)")
)

# ---- 2. Standalone-panel renderer ------------------------------------------

# Single-mode panel set (used for both expression-only and full).
# `d` is a data frame already filtered to that mode's eligible rows.
# `suffix` is "expr_only" or "full", used to pick the right columns.
# `mode_label` is the human-readable title decoration.
render_single_mode <- function(d, suffix, mode_label, name_prefix,
                                colour) {
  prep_col <- paste0("prepare_seconds_",  suffix)
  peak_col <- paste0("prepare_peak_mb_",  suffix)
  load_col <- paste0("load_seconds_",     suffix)
  bnd_col  <- paste0("bundle_mb_",        suffix)

  panels <- list()

  panels$A <- ggplot(d, aes(.data$n_cells_any, .data[[prep_col]])) +
    geom_point(size = 2.0, colour = colour, alpha = 0.85) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                linewidth = 0.5, linetype = "dashed", colour = colour) +
    sx_log_cells +
    scale_y_log10(labels = label_log()) +
    labs(title = sprintf("prepare_study() wall time  [%s]", mode_label),
         y = "Seconds (log10)") +
    theme_paper()

  panels$B <- ggplot(d, aes(.data$n_cells_any, .data[[peak_col]])) +
    geom_point(size = 2.0, colour = colour, alpha = 0.85) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                linewidth = 0.5, linetype = "dashed", colour = colour) +
    sx_log_cells +
    scale_y_log10(labels = label_log()) +
    labs(title = sprintf("prepare_study() peak memory  [%s]", mode_label),
         y = "Peak MB (log10)") +
    theme_paper()

  panels$C <- ggplot(d, aes(.data[[bnd_col]], .data[[load_col]])) +
    geom_point(size = 2.0, colour = colour, alpha = 0.85) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                linewidth = 0.5, linetype = "dashed", colour = colour) +
    scale_x_log10(labels = label_log()) +
    labs(title = sprintf("load_study() cold-start latency  [%s]",
                          mode_label),
         x = "Bundle MB (log10)", y = "Seconds") +
    theme_paper()

  panels$D <- ggplot(d, aes(.data$n_cells_any, .data[[bnd_col]])) +
    geom_point(size = 2.0, colour = colour, alpha = 0.85) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                linewidth = 0.5, linetype = "dashed", colour = colour) +
    sx_log_cells +
    scale_y_log10(labels = label_log()) +
    labs(title = sprintf("Bundle size  [%s]", mode_label),
         y = "Bundle MB (log10)") +
    theme_paper()

  for (k in names(panels)) {
    save_panel(panels[[k]],
               paste0(name_prefix, "_", k),
               width = 5.5, height = 4.0)
  }
  if (requireNamespace("patchwork", quietly = TRUE)) {
    pw <- patchwork::wrap_plots(panels, ncol = 2) +
      patchwork::plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(face = "bold"))
    save_panel(pw, name_prefix, width = 11.5, height = 8.0)
  }
}

# Paired "compare" panel: one slope-style segment per study showing
# expression-only -> full for a given metric. Studies are sorted by
# the full-mode value so the y-axis stays interpretable.
render_compare_panel <- function(d, prep, expr_col, full_col,
                                  y_label, log_y = TRUE,
                                  name) {
  # Long-format frame: two rows per study (one per mode).
  long <- data.frame(
    studyID = rep(as.character(d$studyID), 2L),
    mode    = c(rep("expression-only", nrow(d)),
                rep("full",            nrow(d))),
    value   = c(as.numeric(d[[expr_col]]),
                as.numeric(d[[full_col]])),
    stringsAsFactors = FALSE
  )
  # Order studies left-to-right by full-mode value so the slope direction
  # reads consistently.
  ord <- order(as.numeric(d[[full_col]]))
  long$studyID <- factor(long$studyID,
                          levels = as.character(d$studyID)[ord])
  segs <- data.frame(
    studyID = factor(as.character(d$studyID),
                      levels = as.character(d$studyID)[ord]),
    y_expr  = as.numeric(d[[expr_col]]),
    y_full  = as.numeric(d[[full_col]]),
    stringsAsFactors = FALSE
  )
  p <- ggplot() +
    geom_segment(data = segs,
                  aes(x = .data$studyID, xend = .data$studyID,
                      y = .data$y_expr,  yend = .data$y_full),
                  colour = "grey60", linewidth = 0.4) +
    geom_point(data = long,
                aes(.data$studyID, .data$value,
                    colour = .data$mode, shape = .data$mode),
                size = 2.4, alpha = 0.95) +
    scale_colour_manual(values = c(`expression-only` = npg[["mint"]],
                                    `full`            = npg[["navy"]]),
                        name = NULL) +
    scale_shape_manual(values = c(`expression-only` = 16, `full` = 17),
                       name = NULL) +
    labs(title = prep, x = NULL, y = y_label) +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))
  if (isTRUE(log_y)) {
    p <- p + scale_y_log10(labels = label_log())
  }
  save_panel(p, name, width = 7.5, height = 4.5)
  p
}

# ---- 3. Render: supplemental expression-only set ---------------------------

message("\nRendering supplemental expression-only panels:")
render_single_mode(df_expr, suffix = "expr_only",
                   mode_label  = "expression-only",
                   name_prefix = "figureS_expr_only",
                   colour      = npg[["mint"]])

# ---- 4. Render: paired compare set -----------------------------------------

message("\nRendering paired compare panels:")
cmp_panels <- list()
if (nrow(df_cmp) > 0L) {
  cmp_panels$A <- render_compare_panel(
    df_cmp,
    prep     = "prepare_study() wall time: full vs expression-only",
    expr_col = "prepare_seconds_expr_only",
    full_col = "prepare_seconds_full",
    y_label  = "Seconds (log10)",
    log_y    = TRUE,
    name     = "compare_A_prepare_seconds")
  cmp_panels$B <- render_compare_panel(
    df_cmp,
    prep     = "prepare_study() peak memory: full vs expression-only",
    expr_col = "prepare_peak_mb_expr_only",
    full_col = "prepare_peak_mb_full",
    y_label  = "Peak MB (log10)",
    log_y    = TRUE,
    name     = "compare_B_prepare_peak_mb")
  cmp_panels$C <- render_compare_panel(
    df_cmp,
    prep     = "load_study() cold-start: full vs expression-only",
    expr_col = "load_seconds_expr_only",
    full_col = "load_seconds_full",
    y_label  = "Seconds (log10)",
    log_y    = TRUE,
    name     = "compare_C_load_seconds")
  cmp_panels$D <- render_compare_panel(
    df_cmp,
    prep     = "Bundle size: full vs expression-only",
    expr_col = "bundle_mb_expr_only",
    full_col = "bundle_mb_full",
    y_label  = "Bundle MB (log10)",
    log_y    = TRUE,
    name     = "compare_D_bundle_mb")
  if (requireNamespace("patchwork", quietly = TRUE)) {
    pw <- patchwork::wrap_plots(cmp_panels, ncol = 2) +
      patchwork::plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(face = "bold"))
    save_panel(pw, "compare", width = 15.0, height = 9.0)
  }
} else {
  message("  (no rows with both modes -- skipping compare panels)")
}

message(sprintf(
  "\nDone. Wrote %d figure sets under %s",
  sum(c(nrow(df_expr) > 0L, nrow(df_cmp) > 0L)),
  figures_dir))
