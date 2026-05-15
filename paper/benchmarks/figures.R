#!/usr/bin/env Rscript
# paper/benchmarks/figures.R
#
# Runs the synthetic-sweep benchmarks defined in paper/benchmarks/methods.R
# against:
#   1. a 7x4 grid of synthetic studies (cells x genes scaling),
#   2. the real 2327 (Tex) study at data/2327/2327.scminer.h5 (optional),
#   3. a discover_studies() scan over 1..32 multi-study roots,
# and writes:
#   paper/figures/figure1_A_size_vs_cells.{pdf,png}    bundle + shards
#   paper/figures/figure1_B_load_latency.{pdf,png}     load_study cold start
#   paper/figures/figure1_C_fetch_latency.{pdf,png}    gene_values median + max
#   paper/figures/figure1_D_discover_scaling.{pdf,png} discover_studies()
#   paper/figures/figure1_E_prepare_time.{pdf,png}     prepare_study_data wall
#   paper/figures/figure1_F_peak_memory.{pdf,png}      prepare_study_data peak
#   paper/metrics/bundle_scaling.tsv                   per-config rows
#   paper/metrics/discover_scaling.tsv
#   paper/metrics/real_study.tsv                       (if bundle present)
#
# Run from the project root:
#   Rscript paper/benchmarks/figures.R
#
# Knobs: edit SCALING_GRID below to grow / shrink the synthetic sweep.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(ggsci)
  library(dplyr)
})

source("paper/benchmarks/methods.R")

# Each row in SCALING_GRID becomes one synthetic study x N_REPS metric
# rows. With the 7x4 = 28 grid and N_REPS = 5 below, that's 140 bench
# runs at ~ 10-60 s each -- expect ~ 60-90 min on a laptop and less on
# HPC. Shrink SCALING_GRID or drop N_REPS to iterate faster.
N_REPS <- 5L
SCALING_GRID <- expand.grid(
  n_cells = c(500L, 1000L, 2000L, 4000L, 6000L, 8000L, 10000L),
  n_genes = c(2000L, 5000L, 8000L, 10000L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

DISCOVER_GRID <- c(1L, 2L, 4L, 8L, 16L, 32L)

# ---------------------------------------------------------------- bench loop

run_scaling <- function(grid, scratch_root, n_reps = N_REPS) {
  message(sprintf("Scaling benchmark -- %d configurations x %d reps:",
                  nrow(grid), n_reps))
  out <- list()
  for (i in seq_len(nrow(grid))) {
    row <- grid[i, , drop = FALSE]
    for (rep in seq_len(n_reps)) {
      label <- sprintf("[%d/%d] cells=%d genes=%d  rep %d/%d",
                       i, nrow(grid), row$n_cells, row$n_genes,
                       rep, n_reps)
      message("  ", label)
      sub_root <- file.path(scratch_root,
                              sprintf("c%d_g%d_r%d",
                                      row$n_cells, row$n_genes, rep))
      dir.create(sub_root, recursive = TRUE, showWarnings = FALSE)
      s <- make_synthetic_study(
        n_cells = row$n_cells, n_genes = row$n_genes,
        n_clusters = 4L, density = 0.10,
        # Distinct seed per (config, rep) so replicates are
        # independent draws of the synthetic generator.
        seed = 1000L + i * 100L + rep,
        # Scaling sweep focuses on the expression index + network
        # rows; the real 2327 study (run separately below) covers
        # the full-featured case with activity matrices too.
        with_activity = FALSE, with_networks = TRUE
      )
      b <- as.data.frame(bench_bundle(s, sub_root))
      b$replicate <- rep
      out[[length(out) + 1L]] <- b
    }
  }
  do.call(rbind, out)
}

run_discover <- function(grid, n_reps = N_REPS) {
  message(sprintf("\ndiscover_studies() benchmark x %d reps:", n_reps))
  out <- list()
  for (i in seq_along(grid)) {
    for (rep in seq_len(n_reps)) {
      message(sprintf("  [%d/%d] %d studies  rep %d/%d",
                      i, length(grid), grid[i], rep, n_reps))
      r <- as.data.frame(bench_discover(grid[i],
                                         cells_per = 200L,
                                         genes_per = 500L))
      r$replicate <- rep
      out[[length(out) + 1L]] <- r
    }
  }
  do.call(rbind, out)
}

# Mean + standard error helpers. SE = sd / sqrt(n).
.summarise_reps <- function(df, group_cols, value_cols) {
  keys <- unique(df[, group_cols, drop = FALSE])
  rows <- lapply(seq_len(nrow(keys)), function(k) {
    key <- keys[k, , drop = FALSE]
    mask <- rep(TRUE, nrow(df))
    for (g in group_cols) mask <- mask & df[[g]] == key[[g]]
    sub <- df[mask, , drop = FALSE]
    out <- key
    out$n_reps <- nrow(sub)
    for (v in value_cols) {
      xs <- as.numeric(sub[[v]])
      out[[paste0(v, "_mean")]] <- mean(xs, na.rm = TRUE)
      out[[paste0(v, "_sd")  ]] <- stats::sd(xs, na.rm = TRUE)
      out[[paste0(v, "_se")  ]] <- stats::sd(xs, na.rm = TRUE) /
                                       sqrt(sum(!is.na(xs)))
    }
    out
  })
  do.call(rbind, rows)
}

# ----------------------------------------------------------------------- run

scratch <- tempfile("paper_bench_")
dir.create(scratch, recursive = TRUE)
on.exit(unlink(scratch, recursive = TRUE), add = TRUE)

scaling  <- run_scaling(SCALING_GRID, scratch)
discover <- run_discover(DISCOVER_GRID)
real     <- bench_real_study()

# Aggregate replicates -> mean / sd / SE per config.
scaling_summary <- .summarise_reps(
  scaling,
  group_cols = c("n_cells", "n_genes"),
  value_cols = c("bundle_bytes", "shard_bytes",
                  "prepare_seconds", "prepare_peak_mb",
                  "load_seconds", "fetch_median", "fetch_max")
)
discover_summary <- .summarise_reps(
  discover,
  group_cols = "n_studies",
  value_cols = c("discover_seconds", "n_found")
)

# Persist raw + summary metrics -------------------------------------------
dir.create("paper/metrics", recursive = TRUE, showWarnings = FALSE)
write.table(scaling, "paper/metrics/bundle_scaling.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(scaling_summary, "paper/metrics/bundle_scaling_summary.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(discover, "paper/metrics/discover_scaling.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(discover_summary, "paper/metrics/discover_scaling_summary.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
if (!is.null(real)) {
  write.table(as.data.frame(real, stringsAsFactors = FALSE),
              "paper/metrics/real_study.tsv",
              sep = "\t", row.names = FALSE, quote = FALSE)
}
message(sprintf("\nWrote metrics (%d rep'd rows, %d configs) to paper/metrics/",
                nrow(scaling), nrow(scaling_summary)))

# ---------------------------------------------------------- shared figure deps

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

dir.create("paper/figures", recursive = TRUE, showWarnings = FALSE)

save_panel <- function(plot, name, width, height) {
  pdf_path <- file.path("paper/figures", paste0(name, ".pdf"))
  png_path <- file.path("paper/figures", paste0(name, ".png"))
  ggsave(pdf_path, plot, width = width, height = height,
         units = "in", device = cairo_pdf)
  ggsave(png_path, plot, width = width, height = height,
         units = "in", dpi = 300)
  message("  ", pdf_path, " + ", basename(png_path))
}

# Derived columns on the per-config summary (one row per cells x genes).
scaling_summary$mb_bundle_mean <- scaling_summary$bundle_bytes_mean / 1024^2
scaling_summary$mb_bundle_se   <- scaling_summary$bundle_bytes_se   / 1024^2
scaling_summary$mb_shard_mean  <- scaling_summary$shard_bytes_mean  / 1024^2
scaling_summary$mb_shard_se    <- scaling_summary$shard_bytes_se    / 1024^2
scaling_summary$gene_label <- factor(
  paste0(scaling_summary$n_genes / 1000, "K genes"),
  levels = paste0(sort(unique(scaling_summary$n_genes)) / 1000, "K genes")
)

# ---------------------------------------------------------------- panels

message("\nRendering individual panels (mean +/- SE over ",
        N_REPS, " replicates):")

# A. Bundle vs shard tree size, two series via mk_half().
mk_half <- function(df, mean_col, se_col, kind_label) {
  out <- df[, c("n_cells", "n_genes", "gene_label", mean_col, se_col)]
  names(out)[(ncol(out) - 1L):ncol(out)] <- c("mb", "se")
  out$kind <- kind_label
  out
}
df_long <- rbind(
  mk_half(scaling_summary, "mb_bundle_mean", "mb_bundle_se", "Bundle"),
  mk_half(scaling_summary, "mb_shard_mean",  "mb_shard_se",  "Shard tree")
)
p_A <- ggplot(df_long,
              aes(n_cells, mb, colour = gene_label, linetype = kind,
                  shape = kind)) +
  geom_line(linewidth = 0.5) +
  geom_errorbar(aes(ymin = pmax(mb - se, .Machine$double.eps),
                     ymax = mb + se),
                width = 0, linewidth = 0.4, show.legend = FALSE) +
  geom_point(size = 1.9) +
  scale_x_continuous(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  scale_colour_npg() +
  labs(title = "Bundle size vs shard-tree size",
       x = "Number of cells", y = "MB (log)",
       colour = NULL, linetype = NULL, shape = NULL) +
  theme_paper() +
  guides(colour = guide_legend(order = 1, nrow = 1),
         linetype = guide_legend(order = 2, nrow = 1),
         shape    = guide_legend(order = 2, nrow = 1))
save_panel(p_A, "figure1_A_size_vs_cells",   width = 6.5, height = 4.2)

# B. cold-load latency
p_B <- ggplot(scaling_summary,
              aes(n_cells, load_seconds_mean,
                  colour = gene_label, group = gene_label)) +
  geom_line(linewidth = 0.5) +
  geom_errorbar(aes(ymin = load_seconds_mean - load_seconds_se,
                     ymax = load_seconds_mean + load_seconds_se),
                width = 0, linewidth = 0.4) +
  geom_point(size = 1.9) +
  scale_x_continuous(labels = label_comma()) +
  scale_colour_npg() +
  labs(title = "load_study() cold-start latency",
       x = "Number of cells", y = "Seconds (mean +/- SE)",
       colour = NULL) +
  theme_paper()
save_panel(p_B, "figure1_B_load_latency",    width = 6.0, height = 4.2)

# C. gene_values() median fetch latency (mean +/- SE across replicates).
scaling_summary$fetch_median_ms_mean <- scaling_summary$fetch_median_mean * 1000
scaling_summary$fetch_median_ms_se   <- scaling_summary$fetch_median_se   * 1000
p_C <- ggplot(scaling_summary,
              aes(n_cells, fetch_median_ms_mean,
                  colour = gene_label, group = gene_label)) +
  geom_line(linewidth = 0.5) +
  geom_errorbar(aes(ymin = fetch_median_ms_mean - fetch_median_ms_se,
                     ymax = fetch_median_ms_mean + fetch_median_ms_se),
                width = 0, linewidth = 0.4) +
  geom_point(size = 1.9) +
  scale_x_continuous(labels = label_comma()) +
  scale_colour_npg() +
  labs(title = "gene_values() median fetch latency",
       x = "Number of cells", y = "Milliseconds (mean +/- SE)",
       colour = NULL) +
  theme_paper()
save_panel(p_C, "figure1_C_fetch_latency",   width = 6.0, height = 4.2)

# D. discover_studies() scaling
p_D <- ggplot(discover_summary,
              aes(n_studies, discover_seconds_mean)) +
  geom_line(linewidth = 0.5, colour = "#3C5488") +
  geom_errorbar(aes(ymin = discover_seconds_mean - discover_seconds_se,
                     ymax = discover_seconds_mean + discover_seconds_se),
                width = 0, linewidth = 0.4, colour = "#3C5488") +
  geom_point(size = 1.9, colour = "#3C5488") +
  scale_x_continuous(breaks = DISCOVER_GRID) +
  labs(title = "discover_studies() scaling",
       x = "Studies under root",
       y = "Seconds (mean +/- SE)") +
  theme_paper()
save_panel(p_D, "figure1_D_discover_scaling", width = 6.0, height = 4.2)

# E. prepare_study_data() wall time
p_E <- ggplot(scaling_summary,
              aes(n_cells, prepare_seconds_mean,
                  colour = gene_label, group = gene_label)) +
  geom_line(linewidth = 0.5) +
  geom_errorbar(aes(ymin = prepare_seconds_mean - prepare_seconds_se,
                     ymax = prepare_seconds_mean + prepare_seconds_se),
                width = 0, linewidth = 0.4) +
  geom_point(size = 1.9) +
  scale_x_continuous(labels = label_comma()) +
  scale_colour_npg() +
  labs(title = "prepare_study_data() wall time",
       x = "Number of cells", y = "Seconds (mean +/- SE)",
       colour = NULL) +
  theme_paper()
save_panel(p_E, "figure1_E_prepare_time",    width = 6.0, height = 4.2)

# F. prepare_study_data() peak memory
p_F <- ggplot(scaling_summary,
              aes(n_cells, prepare_peak_mb_mean,
                  colour = gene_label, group = gene_label)) +
  geom_line(linewidth = 0.5) +
  geom_errorbar(aes(ymin = prepare_peak_mb_mean - prepare_peak_mb_se,
                     ymax = prepare_peak_mb_mean + prepare_peak_mb_se),
                width = 0, linewidth = 0.4) +
  geom_point(size = 1.9) +
  scale_x_continuous(labels = label_comma()) +
  scale_colour_npg() +
  labs(title = "prepare_study_data() peak memory",
       x = "Number of cells", y = "Peak Mb (mean +/- SE)",
       colour = NULL) +
  theme_paper()
save_panel(p_F, "figure1_F_peak_memory",     width = 6.0, height = 4.2)

# Clean up any stale combined-grid artifact from earlier renders.
for (legacy in c("paper/figures/figure1.pdf",
                  "paper/figures/figure1.png")) {
  if (file.exists(legacy)) file.remove(legacy)
}

message(sprintf("\nRendered %d standalone panels under paper/figures/", 6L))

# ------------------------------------------------------- summary on stdout

if (!is.null(real)) {
  message("\n=== Real study (2327 / Tex) ===")
  message(sprintf("  bundle:        %s (%s cells x %s genes, %d clusters)",
                  human_bytes(real$bundle_bytes),
                  format(real$n_cells, big.mark = ","),
                  format(real$n_genes, big.mark = ","),
                  real$n_clusters))
  message(sprintf("  load_seconds:  %.3f", real$load_seconds))
  message(sprintf("  fetch_median:  %.1f ms (%d genes fetched)",
                  real$fetch_median * 1000, real$n_fetched))
  message(sprintf("  network edges: TF=%d, SIG=%d",
                  real$net_tf, real$net_sig))
}

message(sprintf(
  "\n=== Synthetic scaling summary (mean across %d reps) ===", N_REPS))
disp <- scaling_summary[, c(
  "n_cells", "n_genes", "n_reps",
  "bundle_bytes_mean", "shard_bytes_mean",
  "load_seconds_mean", "fetch_median_mean",
  "prepare_seconds_mean", "prepare_peak_mb_mean"
)]
names(disp) <- c("n_cells", "n_genes", "n_reps",
                  "bundle_B", "shard_B",
                  "load_s", "fetch_s",
                  "prepare_s", "peak_Mb")
print(disp, row.names = FALSE)
