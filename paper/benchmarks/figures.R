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

# Each row in SCALING_GRID becomes one synthetic study + one row of metrics.
# Expanded 7x4 grid: covers 500-10K cells x 2-10K genes (28 configs).
# At 10K x 10K density 0.10 the in-memory dgCMatrix is ~80 M nnz, well
# under any limit; total wall time on a laptop is ~ 15-25 min for the
# full sweep. Shrink the vectors below for a faster iteration cycle.
SCALING_GRID <- expand.grid(
  n_cells = c(500L, 1000L, 2000L, 4000L, 6000L, 8000L, 10000L),
  n_genes = c(2000L, 5000L, 8000L, 10000L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

DISCOVER_GRID <- c(1L, 2L, 4L, 8L, 16L, 32L)

# ---------------------------------------------------------------- bench loop

run_scaling <- function(grid, scratch_root) {
  message(sprintf("Scaling benchmark -- %d configurations:", nrow(grid)))
  out <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    row <- grid[i, , drop = FALSE]
    label <- sprintf("[%d/%d] cells=%d genes=%d",
                     i, nrow(grid), row$n_cells, row$n_genes)
    message("  ", label)
    sub_root <- file.path(scratch_root,
                            sprintf("c%d_g%d", row$n_cells, row$n_genes))
    dir.create(sub_root, recursive = TRUE, showWarnings = FALSE)
    s <- make_synthetic_study(
      n_cells = row$n_cells, n_genes = row$n_genes,
      n_clusters = 4L, density = 0.10,
      seed = 1000L + i,
      # Scaling sweep focuses on the expression index + network rows;
      # the real 2327 study (run separately below) covers the
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

# Derived columns
scaling$mb_bundle  <- scaling$bundle_bytes / 1024^2
scaling$mb_shard   <- scaling$shard_bytes  / 1024^2
scaling$gene_label <- factor(
  paste0(scaling$n_genes / 1000, "K genes"),
  levels = paste0(sort(unique(scaling$n_genes)) / 1000, "K genes")
)

# ---------------------------------------------------------------- panels

message("\nRendering individual panels:")

# A. Bundle vs shard tree size, two series via mk_half().
mk_half <- function(df, value_col, kind_label) {
  out <- df[, c("n_cells", "n_genes", "gene_label", value_col)]
  names(out)[ncol(out)] <- "mb"
  out$kind <- kind_label
  out
}
df_long <- rbind(
  mk_half(scaling, "mb_bundle", "Bundle"),
  mk_half(scaling, "mb_shard",  "Shard tree")
)
p_A <- ggplot(df_long,
              aes(n_cells, mb, colour = gene_label, linetype = kind,
                  shape = kind)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.9) +
  scale_x_continuous(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  scale_colour_npg() +
  labs(title = "Bundle size vs shard-tree size",
       x = "Number of cells", y = "MB (log)",
       colour = NULL, linetype = NULL, shape = NULL) +
  theme_paper() +
  guides(colour = guide_legend(order = 1),
         linetype = guide_legend(order = 2),
         shape    = guide_legend(order = 2))
save_panel(p_A, "figure1_A_size_vs_cells",   width = 6.5, height = 4.2)

# B. cold-load latency
p_B <- ggplot(scaling, aes(n_cells, load_seconds,
                            colour = gene_label, group = gene_label)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.9) +
  scale_x_continuous(labels = label_comma()) +
  scale_colour_npg() +
  labs(title = "load_study() cold-start latency",
       x = "Number of cells", y = "Seconds",
       colour = NULL) +
  theme_paper()
save_panel(p_B, "figure1_B_load_latency",    width = 6.0, height = 4.2)

# C. gene_values() median + max fetch latency
p_C <- ggplot(scaling, aes(n_cells, fetch_median * 1000,
                            colour = gene_label, group = gene_label)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.9) +
  geom_errorbar(aes(ymin = fetch_median * 1000,
                     ymax = fetch_max    * 1000),
                width = 0, alpha = 0.45) +
  scale_x_continuous(labels = label_comma()) +
  scale_colour_npg() +
  labs(title = "gene_values() fetch latency (median, bar to max)",
       x = "Number of cells", y = "Milliseconds per gene",
       colour = NULL) +
  theme_paper()
save_panel(p_C, "figure1_C_fetch_latency",   width = 6.0, height = 4.2)

# D. discover_studies() scaling
p_D <- ggplot(discover, aes(n_studies, discover_seconds)) +
  geom_line(linewidth = 0.5, colour = "#3C5488") +
  geom_point(size = 1.9, colour = "#3C5488") +
  scale_x_continuous(breaks = DISCOVER_GRID) +
  labs(title = "discover_studies() scaling",
       x = "Studies under root", y = "Seconds") +
  theme_paper()
save_panel(p_D, "figure1_D_discover_scaling", width = 6.0, height = 4.2)

# E. prepare_study_data() wall time
p_E <- ggplot(scaling, aes(n_cells, prepare_seconds,
                            colour = gene_label, group = gene_label)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.9) +
  scale_x_continuous(labels = label_comma()) +
  scale_colour_npg() +
  labs(title = "prepare_study_data() wall time",
       x = "Number of cells", y = "Seconds",
       colour = NULL) +
  theme_paper()
save_panel(p_E, "figure1_E_prepare_time",    width = 6.0, height = 4.2)

# F. prepare_study_data() peak memory
p_F <- ggplot(scaling, aes(n_cells, prepare_peak_mb,
                            colour = gene_label, group = gene_label)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.9) +
  scale_x_continuous(labels = label_comma()) +
  scale_colour_npg() +
  labs(title = "prepare_study_data() peak memory",
       x = "Number of cells", y = "Peak Mb (gc-reported)",
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

message("\n=== Synthetic scaling summary ===")
print(scaling[, c("n_cells", "n_genes", "bundle_bytes", "shard_bytes",
                    "load_seconds", "fetch_median",
                    "prepare_seconds", "prepare_peak_mb")])
