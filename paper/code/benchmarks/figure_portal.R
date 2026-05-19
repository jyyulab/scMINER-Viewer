#!/usr/bin/env Rscript
# paper/code/benchmarks/figure_portal.R
#
# Reads the aggregated per-study metrics produced by the portal-study
# benchmark (paper/metrics/portal_studies/portal_studies.tsv -- written either by
# paper/code/portal/portal_merge.R after a job-array run, or by the
# single-bsub wrapper directly) and renders six standalone figures
# summarising how each metric scales across the real scMINER Portal
# studies. Each panel is its own file so figures can be placed
# independently in the manuscript (or used in supplementary slides).
#
# Outputs (six pairs of {pdf,png} under paper/figures/):
#   figure3_A_size_vs_cells.{pdf,png}   bundle + shard tree size
#   figure3_B_prepare_time.{pdf,png}    wall time of prepare_study()
#   figure3_C_peak_memory.{pdf,png}     R peak memory vs input
#   figure3_D_load_latency.{pdf,png}    cold load_study() latency
#   figure3_E_fetch_latency.{pdf,png}   gene_values() median ms
#   figure3_F_size_ratio.{pdf,png}      output:input ratio per study
# plus:
#   paper/metrics/portal_studies/portal_studies_summary.tsv  one row per status bucket
#
# Run from the project root:
#   Rscript paper/code/benchmarks/figure_portal.R
#
# Falls back to concat'ing every paper/metrics/portal_studies/portal_studies_*.tsv it
# can find if the merged TSV is absent.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

# NPG palette: prefer ggsci's scale_colour_npg() when installed; fall
# back to a manual scale using the same hex codes so the script runs on
# machines without ggsci.
.have_ggsci <- requireNamespace("ggsci", quietly = TRUE)
if (.have_ggsci) {
  suppressPackageStartupMessages(library(ggsci))
}
npg_hex <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
              "#F39B7F", "#91D1C2", "#8491B4")
scale_colour_npg_safe <- function() {
  if (.have_ggsci) ggsci::scale_color_npg()
  else scale_colour_manual(values = npg_hex)
}

# CLI overrides:
#   --metrics-dir <dir>  read TSVs from <dir>/ instead of paper/metrics/
#   --figures-dir <dir>  write figures into <dir>/ instead of paper/figures/
#   --summary-out <tsv>  override the status-count TSV path
.cli <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  i <- which(args == name)
  if (length(i) == 0L || i == length(args)) return(default)
  args[i + 1L]
}
metrics_dir   <- .cli("--metrics-dir", "paper/metrics/portal_studies")
figures_dir   <- .cli("--figures-dir", "paper/figures")
# Figure 3 is the real-portal benchmark restricted to studies that
# actually exercise the full pipeline (expression + activity + TF/SIG
# networks). Source priority:
#   1. Per-mode full-mode TSVs directly under <metrics_dir>/ (produced
#      by portal_studies_compare.sh). Rows are filtered to
#      net_tf_edges + net_sig_edges > 0 -- studies whose YAML lacked
#      networks data have a full row that's byte-identical to
#      expression-only and don't belong in a TF/sig figure.
#   2. <metrics_dir>/portal_studies.tsv (legacy single-mode merge), used
#      only as a fallback when the comparison run hasn't been done yet.
compare_glob  <- file.path(metrics_dir, "portal_studies_*_full.tsv")
merged_path   <- file.path(metrics_dir, "portal_studies.tsv")
# --metrics-tsv lets callers point at any TSV with the per-study schema
# (overrides both the comparison glob and the legacy merged TSV).
override_tsv  <- .cli("--metrics-tsv", NULL)

# ---- 1. Load + sanity-check ------------------------------------------------

read_one <- function(p) {
  utils::read.table(p, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE,
                    check.names = FALSE,
                    na.strings = c("", "NA"))
}

bind_rows_fill <- function(rows) {
  all_cols <- unique(unlist(lapply(rows, colnames)))
  rows <- lapply(rows, function(d) {
    miss <- setdiff(all_cols, colnames(d))
    for (m in miss) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  })
  do.call(rbind, rows)
}

compare_files <- Sys.glob(compare_glob)
df <- if (!is.null(override_tsv)) {
  message("Reading override TSV: ", override_tsv)
  read_one(override_tsv)
} else if (length(compare_files) > 0L) {
  message(sprintf("Reading %d per-mode full TSVs from %s",
                  length(compare_files), dirname(compare_glob)))
  bind_rows_fill(lapply(compare_files, read_one))
} else if (file.exists(merged_path)) {
  message("Reading legacy merged TSV: ", merged_path)
  read_one(merged_path)
} else {
  stop("No metrics found. Looked under:\n",
       "  ", compare_glob, "\n",
       "  ", merged_path, "\n",
       "Run portal_studies_compare.sh + portal_compare.R first, ",
       "or pass --metrics-tsv <path>.")
}

message(sprintf("Loaded %d rows  (status: %s)",
                nrow(df),
                paste(sprintf("%s=%d", names(table(df$status)),
                              as.integer(table(df$status))),
                      collapse = ", ")))

# Local %||% (figure_portal.R is its own script; no global helpers)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

ok <- df[!is.na(df$status) & df$status == "ok", , drop = FALSE]
if (nrow(ok) == 0L) stop("No status=='ok' rows to plot.")

# TFsig filter: when the source is the per-mode comparison glob (or any
# TSV that carries net_tf_edges / net_sig_edges columns), restrict
# Figure 3 to studies whose full-mode run actually consumed a networks
# file. Studies whose YAML lacks networks have a full-mode row that's
# byte-identical to expression-only and would clutter the plot.
has_tfsig_cols <- all(c("net_tf_edges", "net_sig_edges") %in% names(ok))
if (has_tfsig_cols) {
  tfsig_n <- as.numeric(ok$net_tf_edges %||% 0) +
    as.numeric(ok$net_sig_edges %||% 0)
  keep_tfsig <- !is.na(tfsig_n) & tfsig_n > 0
  if (any(keep_tfsig)) {
    message(sprintf("  filtering to %d TF/sig-eligible studies (dropping %d)",
                    sum(keep_tfsig), sum(!keep_tfsig)))
    ok <- ok[keep_tfsig, , drop = FALSE]
  } else {
    message("  no TF/sig edges in any row; keeping all status=ok rows")
  }
}

# Derived helpers
ok$studyID      <- as.character(ok$studyID)
ok$cells_genes  <- as.numeric(ok$n_cells) * as.numeric(ok$n_genes)
ok$bundle_mb    <- ok$bundle_bytes       / 1024^2
ok$shard_mb     <- ok$shard_bytes        / 1024^2
ok$in_mb        <- ok$total_input_bytes  / 1024^2
ok$out_mb       <- ok$total_output_bytes / 1024^2
ok$fetch_ms     <- ok$fetch_median * 1000
ok$ratio_out_in <- ok$total_output_bytes / pmax(1, ok$total_input_bytes)
ok$has_activity <- ok$act_input_bytes > 0
ok <- ok[order(ok$n_cells), , drop = FALSE]

# Per-status summary TSV (cited from the manuscript text)
status_summary <- as.data.frame(table(df$status, useNA = "ifany"),
                                  stringsAsFactors = FALSE)
names(status_summary) <- c("status", "n")
utils::write.table(status_summary,
                   file.path(metrics_dir, "portal_studies_summary.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote ", file.path(metrics_dir, "portal_studies_summary.tsv"))

# ---- 2. Shared theme + saver ----------------------------------------------

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

# Probe whether cairo_pdf actually works at this R install. Some Mac
# installs report capabilities("cairo") == TRUE but emit a "failed to
# load cairo DLL" warning at device creation, which silently truncates
# the file. Fall back to the default pdf() device when that happens.
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

# ---- 3. Panels (one save_panel() call each) -------------------------------

message("\nRendering individual panels:")

# A. Bundle + shard size vs cells
df_size <- rbind(
  data.frame(studyID = ok$studyID, n_cells = ok$n_cells,
             mb = ok$bundle_mb,
             kind = "Bundle (.scminer.h5)",
             stringsAsFactors = FALSE),
  data.frame(studyID = ok$studyID, n_cells = ok$n_cells,
             mb = ok$shard_mb,
             kind = "Shard tree (.csv.gz)",
             stringsAsFactors = FALSE)
)
p_A <- ggplot(df_size, aes(n_cells, mb, colour = kind, shape = kind)) +
  geom_point(size = 2.0, alpha = 0.85) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.5, linetype = "dashed", alpha = 0.6) +
  scale_x_log10(labels = label_log()) +
  scale_y_log10(labels = label_log()) +
  scale_colour_npg_safe() +
  labs(title = "Bundle + shard tree size vs cell count",
       x = "Cells (log10)", y = "MB (log10)",
       colour = NULL, shape = NULL) +
  theme_paper()
save_panel(p_A, "figure3_A_size_vs_cells",   width = 5.5, height = 4.0)

# B. prepare time vs cells x genes
p_B <- ggplot(ok, aes(cells_genes, prepare_seconds)) +
  geom_point(size = 2.0, colour = "#E64B35", alpha = 0.85) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.5, linetype = "dashed", colour = "#E64B35") +
  scale_x_log10(labels = label_log()) +
  scale_y_log10(labels = label_log()) +
  labs(title = "prepare_study() wall time",
       x = "n_cells x n_genes (log10)", y = "Seconds (log10)") +
  theme_paper()
save_panel(p_B, "figure3_B_prepare_time",    width = 5.5, height = 4.0)

# C. peak memory vs input size
p_C <- ggplot(ok, aes(in_mb, prepare_peak_mb)) +
  geom_point(size = 2.0, colour = "#00A087", alpha = 0.85) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              colour = "grey50", linewidth = 0.4) +
  scale_x_log10(labels = label_log()) +
  scale_y_log10(labels = label_log()) +
  labs(title = "Peak memory vs total input size",
       x = "Input MB (log10)", y = "Peak Mb (log10)") +
  theme_paper()
save_panel(p_C, "figure3_C_peak_memory",     width = 5.5, height = 4.0)

# D. cold load latency vs bundle size
p_D <- ggplot(ok, aes(bundle_mb, load_seconds)) +
  geom_point(size = 2.0, colour = "#3C5488", alpha = 0.85) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.5, linetype = "dashed", colour = "#3C5488") +
  scale_x_log10(labels = label_log()) +
  labs(title = "load_study() cold-start latency",
       x = "Bundle MB (log10)", y = "Seconds") +
  theme_paper()
save_panel(p_D, "figure3_D_load_latency",    width = 5.5, height = 4.0)

# E. gene-fetch latency vs gene count
p_E <- ggplot(ok, aes(n_genes, fetch_ms)) +
  geom_point(size = 2.0, colour = "#F39B7F", alpha = 0.85) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.5, linetype = "dashed", colour = "#F39B7F") +
  scale_x_log10(labels = label_log()) +
  scale_y_log10(labels = label_log()) +
  labs(title = "gene_values() fetch latency (median)",
       x = "Genes (log10)", y = "Milliseconds (log10)") +
  theme_paper()
save_panel(p_E, "figure3_E_fetch_latency",   width = 5.5, height = 4.0)

# F. output : input ratio per study
ok$id_label <- factor(ok$studyID,
                       levels = ok$studyID[order(ok$ratio_out_in)])
p_F <- ggplot(ok, aes(id_label, ratio_out_in, fill = has_activity)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dotted",
             colour = "grey50", linewidth = 0.4) +
  # Name the label vector explicitly so the legend stays correct even
  # when has_activity has only one level (e.g. Figure 3's TF/sig-only
  # subset where every study has activity).
  scale_fill_manual(values = c(`FALSE` = "#91D1C2",
                                `TRUE`  = "#3C5488"),
                    labels = c(`FALSE` = "expression only",
                                `TRUE`  = "expression + activity"),
                    name = NULL,
                    guide = if (length(unique(ok$has_activity)) > 1L)
                              guide_legend() else "none") +
  labs(title = "Output : input size ratio per study",
       x = NULL, y = "total_output / total_input") +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))
save_panel(p_F, "figure3_F_size_ratio",      width = 7.5, height = 4.0)

# Combined 3x2 patchwork grid -- matches the manuscript's "Figure 3"
# (real portal sweep). Standalone panels above remain available for
# any case where only one panel is needed.
if (requireNamespace("patchwork", quietly = TRUE)) {
  suppressPackageStartupMessages(library(patchwork))
  fig3 <- (p_A | p_B) / (p_C | p_D) / (p_E | p_F) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold"))
  combo_pdf <- file.path(figures_dir, "figure3.pdf")
  combo_png <- file.path(figures_dir, "figure3.png")
  ggsave(combo_pdf, fig3,
         width = 11.5, height = 12.0, units = "in", device = pdf_device)
  ggsave(combo_png, fig3,
         width = 11.5, height = 12.0, units = "in", dpi = 300)
  message("  ", combo_pdf, " + figure3.png (combined 3x2)")
} else {
  message("(patchwork not installed; skipping combined figure3)")
}

# ---- 4. Console summary ---------------------------------------------------

message(sprintf("\nRendered %d standalone panels + 1 combined under %s/",
                6L, figures_dir))
message("\n=== Portal benchmark summary (n=", nrow(ok), " OK rows) ===")
disp <- ok[, c("studyID", "n_cells", "n_genes", "n_clusters",
               "in_mb", "out_mb", "bundle_mb", "shard_mb",
               "prepare_seconds", "prepare_peak_mb",
               "load_seconds", "fetch_ms")]
disp$in_mb     <- round(disp$in_mb,    1)
disp$out_mb    <- round(disp$out_mb,   1)
disp$bundle_mb <- round(disp$bundle_mb, 1)
disp$shard_mb  <- round(disp$shard_mb,  1)
disp$prepare_seconds <- round(disp$prepare_seconds, 1)
disp$prepare_peak_mb <- round(disp$prepare_peak_mb, 0)
disp$load_seconds    <- round(disp$load_seconds, 3)
disp$fetch_ms        <- round(disp$fetch_ms, 1)
print(disp, row.names = FALSE)

message(sprintf(
  paste0("\nTotals: input %s-%s MB, output %s-%s MB, ",
         "prepare %.0f-%.0f s, peak %s-%s MB, ",
         "load %.2f-%.2f s, fetch %s-%s ms"),
  format(round(min(ok$in_mb), 1),    big.mark = ","),
  format(round(max(ok$in_mb), 1),    big.mark = ","),
  format(round(min(ok$out_mb), 1),   big.mark = ","),
  format(round(max(ok$out_mb), 1),   big.mark = ","),
  min(ok$prepare_seconds),  max(ok$prepare_seconds),
  format(round(min(ok$prepare_peak_mb), 0), big.mark = ","),
  format(round(max(ok$prepare_peak_mb), 0), big.mark = ","),
  min(ok$load_seconds),     max(ok$load_seconds),
  format(round(min(ok$fetch_ms), 1), big.mark = ","),
  format(round(max(ok$fetch_ms), 1), big.mark = ",")
))
