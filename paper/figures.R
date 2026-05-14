#!/usr/bin/env Rscript
# paper/figures.R
#
# Runs the benchmarks defined in paper/methods.R against
#   1. a sweep of synthetic studies (cells × genes scaling),
#   2. the real 2327 (Tex) study at data/2327/2327.scminer.h5,
#   3. a discover_studies() scan over 1..32 multi-study roots,
# and writes:
#   paper/figures/figure1.pdf        single multi-panel figure
#   paper/figures/figure1.png        300 dpi rendering
#   paper/metrics/bundle_scaling.tsv per-study numbers
#   paper/metrics/discover_scaling.tsv
#   paper/metrics/real_study.tsv
#
# Run from the project root:
#   Rscript paper/figures.R
#
# Knobs: edit SCALING_GRID below to grow / shrink the synthetic sweep.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(ggsci)
  library(dplyr)
})

source("paper/methods.R")

# Each row in SCALING_GRID becomes one synthetic study + one row of metrics.
# Keep the sweep modest so the benchmark finishes in <2 min on a laptop.
SCALING_GRID <- expand.grid(
  n_cells = c(500L, 1000L, 2000L, 4000L),
  n_genes = c(2000L, 5000L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

DISCOVER_GRID <- c(1L, 2L, 4L, 8L, 16L, 32L)

# ---------------------------------------------------------------- bench loop

run_scaling <- function(grid, scratch_root) {
  message(sprintf("Scaling benchmark — %d configurations:", nrow(grid)))
  out <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    row <- grid[i, , drop = FALSE]
    label <- sprintf("[%d/%d] cells=%d genes=%d",
                     i, nrow(grid), row$n_cells, row$n_genes)
    message("  ", label)
    sub_root <- file.path(scratch_root, sprintf("c%d_g%d",
                                                 row$n_cells, row$n_genes))
    dir.create(sub_root, recursive = TRUE, showWarnings = FALSE)
    s <- make_synthetic_study(
      n_cells = row$n_cells, n_genes = row$n_genes,
      n_clusters = 4L, density = 0.10,
      seed = 1000L + i,
      # Scaling sweep focuses on the expression index + network rows.
      # The real 2327 study (run separately below) covers the
      # full-featured case with activity matrices too.
      with_activity = FALSE, with_networks = TRUE
    )
    b <- bench_bundle(s, sub_root)
    out[[i]] <- as.data.frame(b)
  }
  do.call(rbind, out)
}

run_discover <- function(grid) {
  message("\ndiscover_studies() benchmark:")
  out <- vector("list", length(grid))
  for (i in seq_along(grid)) {
    message(sprintf("  [%d/%d] %d studies", i, length(grid), grid[i]))
    out[[i]] <- as.data.frame(bench_discover(grid[i],
                                              cells_per = 200L,
                                              genes_per = 500L))
  }
  do.call(rbind, out)
}

# ----------------------------------------------------------------------- run

scratch <- tempfile("paper_bench_")
dir.create(scratch, recursive = TRUE)
on.exit(unlink(scratch, recursive = TRUE), add = TRUE)

scaling <- run_scaling(SCALING_GRID, scratch)
discover <- run_discover(DISCOVER_GRID)
real <- bench_real_study()

# Persist raw metrics ------------------------------------------------------
dir.create("paper/metrics", recursive = TRUE, showWarnings = FALSE)
write.table(scaling, "paper/metrics/bundle_scaling.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(discover, "paper/metrics/discover_scaling.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
if (!is.null(real)) {
  write.table(as.data.frame(real, stringsAsFactors = FALSE),
              "paper/metrics/real_study.tsv",
              sep = "\t", row.names = FALSE, quote = FALSE)
}

message("\nWrote metrics to paper/metrics/")

# ------------------------------------------------------------------- figure

theme_paper <- function() {
  theme_classic(base_size = 9) +
    theme(
      plot.title       = element_text(face = "bold", size = 9),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      legend.key.size  = unit(8, "pt"),
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 8),
      axis.text        = element_text(size = 8),
      axis.title       = element_text(size = 8),
      plot.tag         = element_text(face = "bold", size = 11),
      plot.tag.position= c(0.02, 0.98)
    )
}

scaling$mb_bundle <- scaling$bundle_bytes / 1024^2
scaling$mb_shard  <- scaling$shard_bytes  / 1024^2
scaling$total_cells_genes <- scaling$n_cells * scaling$n_genes
scaling$gene_label <- factor(paste0(scaling$n_genes / 1000, "K genes"),
                              levels = paste0(sort(unique(scaling$n_genes)) / 1000,
                                              "K genes"))

# Panel A — bundle vs shard tree size. Build the two halves with
# matched column names so rbind / pivot is happy.
mk_half <- function(df, value_col, kind_label) {
  out <- df[, c("n_cells", "n_genes", "gene_label", value_col)]
  names(out)[ncol(out)] <- "mb"
  out$kind <- kind_label
  out
}
df_long <- rbind(
  mk_half(scaling, "mb_bundle", "Bundle (metadata + indexes)"),
  mk_half(scaling, "mb_shard",  "Shard tree (per-gene .csv.gz)")
)

p_A <- ggplot(df_long,
              aes(n_cells, mb, colour = gene_label, linetype = kind,
                  shape = kind)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.6) +
  scale_x_continuous(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  scale_colour_npg() +
  labs(title = "A. Bundle size vs shard-tree size",
       x = "Number of cells", y = "MB (log)",
       colour = NULL, linetype = NULL, shape = NULL) +
  theme_paper()

# Panel B — cold-load latency
p_B <- ggplot(scaling, aes(n_cells, load_seconds,
                            colour = gene_label, group = gene_label)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.6) +
  scale_x_continuous(labels = label_comma()) +
  scale_colour_npg() +
  labs(title = "B. load_study() cold-start latency",
       x = "Number of cells", y = "Seconds",
       colour = NULL) +
  theme_paper()

# Panel C — first-gene fetch latency
p_C <- ggplot(scaling, aes(n_cells, fetch_median * 1000,
                            colour = gene_label, group = gene_label)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.6) +
  geom_errorbar(aes(ymin = fetch_median * 1000,
                     ymax = fetch_max    * 1000),
                width = 0, alpha = 0.4) +
  scale_x_continuous(labels = label_comma()) +
  scale_colour_npg() +
  labs(title = "C. gene_values() fetch latency (median, bar to max)",
       x = "Number of cells", y = "Milliseconds per gene",
       colour = NULL) +
  theme_paper()

# Panel D — discover_studies() scaling
p_D <- ggplot(discover, aes(n_studies, discover_seconds)) +
  geom_line(linewidth = 0.5, colour = "#3C5488") +
  geom_point(size = 1.6, colour = "#3C5488") +
  scale_x_continuous(breaks = DISCOVER_GRID) +
  labs(title = "D. discover_studies() scaling",
       x = "Studies under root", y = "Seconds") +
  theme_paper()

# Panel E — prepare_study_data() wall time
p_E <- ggplot(scaling, aes(n_cells, prepare_seconds,
                            colour = gene_label, group = gene_label)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.6) +
  scale_x_continuous(labels = label_comma()) +
  scale_colour_npg() +
  labs(title = "E. prepare_study_data() wall time",
       x = "Number of cells", y = "Seconds",
       colour = NULL) +
  theme_paper()

# Panel F — prepare_study_data() peak resident memory
p_F <- ggplot(scaling, aes(n_cells, prepare_peak_mb,
                            colour = gene_label, group = gene_label)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.6) +
  scale_x_continuous(labels = label_comma()) +
  scale_colour_npg() +
  labs(title = "F. prepare_study_data() peak memory",
       x = "Number of cells", y = "Peak Mb (gc-reported)",
       colour = NULL) +
  theme_paper()

fig <- (p_A | p_B) / (p_C | p_D) / (p_E | p_F) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

dir.create("paper/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("paper/figures/figure1.pdf", fig,
       width = 7.0, height = 8.1, units = "in", device = cairo_pdf)
ggsave("paper/figures/figure1.png", fig,
       width = 7.0, height = 8.1, units = "in", dpi = 300)

message("\nFigure written to paper/figures/figure1.{pdf,png}")

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

message("\n=== Synthetic scaling summary ===")
print(scaling[, c("n_cells", "n_genes", "bundle_bytes", "shard_bytes",
                    "load_seconds", "fetch_median")])
