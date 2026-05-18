#!/usr/bin/env Rscript
# paper/benchmarks/tables.R
#
# Generate publication-ready tables from the benchmark TSVs:
#
#   paper/metrics/bundle_scaling.tsv      (synthetic sweep -- figure 1)
#   paper/metrics/portal_studies.tsv      (real 26+ portal studies)
#
# Writes for each:
#   paper/tables/<name>.md     markdown (pandoc-renderable in the manuscript)
#   paper/tables/<name>.tsv    clean TSV (machine-readable, paper-friendly)
#
# Run from the project root:
#   Rscript paper/benchmarks/tables.R
#
# Idempotent: regenerate whenever the source TSVs change.

dir.create("paper/tables", recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(p) {
  if (!file.exists(p)) {
    stop("Required TSV missing: ", p)
  }
  utils::read.table(p, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE,
                    check.names = FALSE,
                    na.strings = c("", "NA"))
}

# Local %||% helper (R < 4.4 doesn't ship one). Defined here so it's
# in scope for all downstream column-fallback uses.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Same loader as figure_portal.R: prefer the per-mode full TSVs from the
# compare benchmark (filtered to TF/sig-eligible studies); fall back to
# the legacy single-mode merged TSV when the comparison run isn't there.
bind_rows_fill <- function(rows) {
  all_cols <- unique(unlist(lapply(rows, colnames)))
  rows <- lapply(rows, function(d) {
    miss <- setdiff(all_cols, colnames(d))
    for (m in miss) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  })
  do.call(rbind, rows)
}

load_portal_metrics <- function(metrics_dir = "paper/metrics") {
  compare_glob <- file.path(metrics_dir, "comparison",
                             "portal_studies_*_full.tsv")
  merged_path  <- file.path(metrics_dir, "portal_studies.tsv")
  compare_files <- Sys.glob(compare_glob)
  if (length(compare_files) > 0L) {
    message(sprintf("Reading %d per-mode full TSVs for Table 1",
                    length(compare_files)))
    bind_rows_fill(lapply(compare_files, read_tsv))
  } else {
    message("No comparison/*_full.tsv; falling back to ", merged_path)
    read_tsv(merged_path)
  }
}

# ---- Markdown helpers -----------------------------------------------------

# Build a GFM-style markdown table. `align` accepts "l", "r", or "c"
# per column; default right-align for numeric columns.
md_table <- function(df, align = NULL) {
  cols <- colnames(df)
  if (is.null(align)) {
    align <- ifelse(vapply(df, is.numeric, logical(1)), "r", "l")
  }
  sep_for <- function(a) switch(a, l = ":---", r = "---:", c = ":---:", "---")
  # Render each cell -- numerics get pretty formatting via the caller.
  render_row <- function(row) {
    paste0("| ", paste(as.character(row), collapse = " | "), " |")
  }
  header <- render_row(cols)
  sep    <- paste0("| ", paste(vapply(align, sep_for, ""),
                                collapse = " | "), " |")
  body   <- vapply(seq_len(nrow(df)), function(i) {
    render_row(df[i, , drop = TRUE])
  }, character(1))
  c(header, sep, body)
}

fmt_n     <- function(x, big = ",") formatC(x, format = "d", big.mark = big)
fmt_mb    <- function(x, d = 1)
  ifelse(is.na(x), "—", formatC(x, format = "f", digits = d, big.mark = ","))
fmt_secs  <- function(x, d = 1)
  ifelse(is.na(x), "—", formatC(x, format = "f", digits = d, big.mark = ","))
fmt_ms    <- function(x, d = 1)
  ifelse(is.na(x), "—", formatC(x, format = "f", digits = d, big.mark = ","))

# ---- 1. Portal studies (real data) -----------------------------------------

portal <- load_portal_metrics("paper/metrics")
portal_ok <- portal[!is.na(portal$status) & portal$status == "ok",
                     , drop = FALSE]
# Table 1 mirrors Figure 3: restrict to TF/sig-eligible studies whenever
# the source carries network-edge columns. Falls through unchanged on
# legacy merged TSVs that don't have those columns.
if (all(c("net_tf_edges", "net_sig_edges") %in% names(portal_ok))) {
  tfsig_n <- as.numeric(portal_ok$net_tf_edges %||% 0) +
    as.numeric(portal_ok$net_sig_edges %||% 0)
  keep <- !is.na(tfsig_n) & tfsig_n > 0
  if (any(keep)) {
    message(sprintf("Table 1: restricting to %d TF/sig-eligible studies (dropping %d)",
                    sum(keep), sum(!keep)))
    portal_ok <- portal_ok[keep, , drop = FALSE]
  }
}
portal_ok <- portal_ok[order(portal_ok$n_cells), , drop = FALSE]

portal_tbl <- data.frame(
  studyID        = as.character(portal_ok$studyID),
  n_cells        = fmt_n(portal_ok$n_cells),
  n_genes        = fmt_n(portal_ok$n_genes),
  n_clusters     = fmt_n(portal_ok$n_clusters),
  `input MB`     = fmt_mb(portal_ok$total_input_bytes  / 1024^2),
  `bundle MB`    = fmt_mb(portal_ok$bundle_bytes       / 1024^2),
  `shard MB`     = fmt_mb(portal_ok$shard_bytes        / 1024^2),
  `output MB`    = fmt_mb(portal_ok$total_output_bytes / 1024^2),
  `prepare s`    = fmt_secs(portal_ok$prepare_seconds, 0),
  `peak Mb`      = fmt_mb(portal_ok$prepare_peak_mb, 0),
  `load s`       = fmt_secs(portal_ok$load_seconds,    2),
  `fetch ms`     = fmt_ms  (portal_ok$fetch_median * 1000, 0),
  `TF edges`     = fmt_n(portal_ok$net_tf_edges %||% 0),
  `SIG edges`    = fmt_n(portal_ok$net_sig_edges %||% 0),
  check.names    = FALSE,
  stringsAsFactors = FALSE
)
# Mask zero edge counts so the table doesn't read as 0 vs missing.
portal_tbl$`TF edges`  <- gsub("^0$", "—", portal_tbl$`TF edges`)
portal_tbl$`SIG edges` <- gsub("^0$", "—", portal_tbl$`SIG edges`)

# Write TSV
write.table(portal_tbl, "paper/tables/figure3_portal_studies.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote paper/tables/figure3_portal_studies.tsv")

# Write Markdown
md_lines <- c(
  sprintf("**Table 1.** Per-study metrics for the %d real scMINER ",
          nrow(portal_ok)),
  "Portal studies, sorted ascending by `n_cells`. `input MB` is the",
  "sum of `expression.rds + activity.rds + networks.txt`. `bundle MB`",
  "is the lazy `.scminer.h5`; `shard MB` is the per-gene gzipped CSV tree.",
  "`prepare s` is wall time of `prepare_study_from_eset()`; `peak Mb`",
  "is the R-reported peak working-set memory during that call.",
  "`fetch ms` is the median time of `gene_values()` over 25 randomly",
  "sampled genes per study. Studies with empty edge cells lack a",
  "networks input.",
  "",
  md_table(portal_tbl)
)
writeLines(md_lines, "paper/tables/figure3_portal_studies.md")
message("Wrote paper/tables/figure3_portal_studies.md")

# ---- 2. Figure 1 synthetic-sweep table -------------------------------------

# Prefer the per-config summary (one row per cells x genes with mean +
# SE across replicates). Fall back to the raw per-rep TSV when the
# summary hasn't been generated yet (single-rep mode).
summary_path <- "paper/metrics/bundle_scaling_summary.tsv"
if (file.exists(summary_path)) {
  s <- read_tsv(summary_path)
  s <- s[order(s$n_cells, s$n_genes), , drop = FALSE]

  # Cell formatter: "mean ± SE" with units already baked in.
  fmt_pm <- function(mean, se, digits = 2, big = ",") {
    ifelse(is.na(mean), "—",
           sprintf("%s ± %s",
                   formatC(mean, format = "f", digits = digits,
                            big.mark = big),
                   formatC(se,   format = "f", digits = digits,
                            big.mark = big)))
  }

  n_reps_value <- if ("n_reps" %in% colnames(s)) max(s$n_reps) else 1L
  fig1_tbl <- data.frame(
    n_cells       = fmt_n(s$n_cells),
    n_genes       = fmt_n(s$n_genes),
    n_reps        = fmt_n(s$n_reps %||% n_reps_value),
    `bundle MB`   = fmt_pm(s$bundle_bytes_mean / 1024^2,
                            s$bundle_bytes_se   / 1024^2, digits = 2),
    `shard MB`    = fmt_pm(s$shard_bytes_mean  / 1024^2,
                            s$shard_bytes_se    / 1024^2, digits = 1),
    `prepare s`   = fmt_pm(s$prepare_seconds_mean,
                            s$prepare_seconds_se, digits = 1),
    `peak Mb`     = fmt_pm(s$prepare_peak_mb_mean,
                            s$prepare_peak_mb_se, digits = 0),
    `load s`      = fmt_pm(s$load_seconds_mean,
                            s$load_seconds_se, digits = 3),
    `fetch median ms` = fmt_pm(s$fetch_median_mean * 1000,
                                  s$fetch_median_se   * 1000, digits = 1),
    `fetch max ms`    = fmt_pm(s$fetch_max_mean    * 1000,
                                  s$fetch_max_se      * 1000, digits = 1),
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  source_label <- sprintf("(mean ± SE across %d replicates)",
                          n_reps_value)
} else {
  message("No summary TSV; falling back to per-rep table at ",
          "paper/metrics/bundle_scaling.tsv")
  scaling <- read_tsv("paper/metrics/bundle_scaling.tsv")
  scaling <- scaling[order(scaling$n_cells, scaling$n_genes), ,
                       drop = FALSE]
  fig1_tbl <- data.frame(
    n_cells       = fmt_n(scaling$n_cells),
    n_genes       = fmt_n(scaling$n_genes),
    `bundle MB`   = fmt_mb(scaling$bundle_bytes / 1024^2, 2),
    `shard MB`    = fmt_mb(scaling$shard_bytes  / 1024^2, 1),
    `prepare s`   = fmt_secs(scaling$prepare_seconds, 1),
    `peak Mb`     = fmt_mb(scaling$prepare_peak_mb, 0),
    `load s`      = fmt_secs(scaling$load_seconds, 3),
    `fetch median ms` = fmt_ms(scaling$fetch_median * 1000, 1),
    `fetch max ms`    = fmt_ms(scaling$fetch_max    * 1000, 1),
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  source_label <- "(single run per configuration)"
}
write.table(fig1_tbl, "paper/tables/figure2_scaling.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote paper/tables/figure2_scaling.tsv")

md_lines2 <- c(
  sprintf("**Table 2.** Synthetic-sweep benchmark (%d configurations, %s) ",
          nrow(fig1_tbl), source_label),
  "underlying figure 2. Each configuration is a synthetic study generated",
  "by `make_synthetic_study(n_cells, n_genes, n_clusters = 4, density = 0.10)`",
  "and benchmarked end-to-end via `bench_bundle()`: bundle + shard write,",
  "`load_study()` cold-start, and 25 random `gene_values()` fetches.",
  "Each replicate uses an independent random seed.",
  "Bundles stay near-constant in size while the shard tree grows linearly",
  "with `n_cells × n_genes`.",
  "",
  md_table(fig1_tbl)
)
writeLines(md_lines2, "paper/tables/figure2_scaling.md")
message("Wrote paper/tables/figure2_scaling.md")

# ---- 3. Console summary ---------------------------------------------------

message(sprintf("\nGenerated tables:"))
message(sprintf("  Table 1 (portal): %d rows", nrow(portal_tbl)))
message(sprintf("  Table 2 (fig 1):  %d rows", nrow(fig1_tbl)))
