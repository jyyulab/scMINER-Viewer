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

portal <- read_tsv("paper/metrics/portal_studies.tsv")
portal_ok <- portal[!is.na(portal$status) & portal$status == "ok",
                     , drop = FALSE]
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

# Local %||% helper (R < 4.4 doesn't ship one)
`%||%` <- function(a, b) if (is.null(a)) b else a

# Write TSV
write.table(portal_tbl, "paper/tables/portal_studies.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote paper/tables/portal_studies.tsv")

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
writeLines(md_lines, "paper/tables/portal_studies.md")
message("Wrote paper/tables/portal_studies.md")

# ---- 2. Figure 1 synthetic-sweep table -------------------------------------

scaling <- read_tsv("paper/metrics/bundle_scaling.tsv")
scaling <- scaling[order(scaling$n_cells, scaling$n_genes), , drop = FALSE]

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
  check.names   = FALSE,
  stringsAsFactors = FALSE
)
write.table(fig1_tbl, "paper/tables/figure1_scaling.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote paper/tables/figure1_scaling.tsv")

md_lines2 <- c(
  sprintf("**Table 2.** Synthetic-sweep benchmark (%d configurations) ",
          nrow(scaling)),
  "underlying figure 1. Each row is a single synthetic study generated",
  "by `make_synthetic_study(n_cells, n_genes, n_clusters = 4, density = 0.10)`",
  "and benchmarked end-to-end via `bench_bundle()`: bundle + shard write,",
  "`load_study()` cold-start, and 25 random `gene_values()` fetches.",
  "Bundles stay near-constant in size while the shard tree grows linearly",
  "with `n_cells × n_genes`.",
  "",
  md_table(fig1_tbl)
)
writeLines(md_lines2, "paper/tables/figure1_scaling.md")
message("Wrote paper/tables/figure1_scaling.md")

# ---- 3. Console summary ---------------------------------------------------

message(sprintf("\nGenerated tables:"))
message(sprintf("  Table 1 (portal): %d rows", nrow(portal_tbl)))
message(sprintf("  Table 2 (fig 1):  %d rows", nrow(fig1_tbl)))
