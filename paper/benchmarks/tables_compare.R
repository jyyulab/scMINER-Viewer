#!/usr/bin/env Rscript
# paper/benchmarks/tables_compare.R
#
# Generate two tables from the joined wide TSV produced by
# paper/portal/portal_compare.R. Table 1 (the manuscript table for
# the TF/sig-eligible studies) is built by paper/benchmarks/tables.R
# from the same compare data, so it isn't duplicated here.
#
#   1. paper/tables/tableS_expr_only.{md,tsv}      (Supplemental table)
#      One row per study with a valid expression-only run. Columns:
#      n_cells, n_genes, prepare_seconds, prepare_peak_mb, load_seconds,
#      bundle_mb. The "no TF/sig, no activity" baseline by itself.
#
#   2. paper/tables/compare_delta.{md,tsv}         (Paired comparison)
#      The TFsig-eligible studies, side by side: each metric appears
#      as (full, expr_only, ratio = full / expr_only).
#
# Run from the project root (after the compare TSV exists):
#   Rscript paper/benchmarks/tables_compare.R

dir.create("paper/tables", recursive = TRUE, showWarnings = FALSE)

.cli <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  i <- which(args == name)
  if (length(i) == 0L || i == length(args)) return(default)
  args[i + 1L]
}
compare_tsv <- .cli("--compare-tsv",
                     "paper/metrics/portal_studies_compare.tsv")
tables_dir  <- .cli("--tables-dir", "paper/tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(compare_tsv)) {
  stop("Compare TSV not found: ", compare_tsv,
       "\nRun paper/portal/portal_compare.R first.")
}
df <- utils::read.table(compare_tsv, sep = "\t", header = TRUE,
                         stringsAsFactors = FALSE,
                         check.names = FALSE,
                         na.strings = c("", "NA"))
message(sprintf("Loaded %d rows from %s", nrow(df), compare_tsv))

# Coalesce shape columns: prefer the full-mode value, fall back to
# expression-only when full is missing (studies where only one mode
# ran -- e.g. 2322 has no full row, so n_cells_full is NA but
# n_cells_expr_only is populated). Without this, those rows show up
# in the table with empty shape cells even though the run succeeded.
.coalesce <- function(d, base) {
  a <- paste0(base, "_full")
  b <- paste0(base, "_expr_only")
  full <- if (a %in% names(d)) as.numeric(d[[a]]) else rep(NA_real_, nrow(d))
  expr <- if (b %in% names(d)) as.numeric(d[[b]]) else rep(NA_real_, nrow(d))
  ifelse(is.na(full), expr, full)
}
df$n_cells_any    <- .coalesce(df, "n_cells")
df$n_genes_any    <- .coalesce(df, "n_genes")
df$n_clusters_any <- .coalesce(df, "n_clusters")

# ---- Slicing ---------------------------------------------------------------

# The supplemental set requires BOTH modes to have produced metrics --
# studies where the full-mode run is missing (e.g. 2322) are excluded
# so the supplemental table only reports studies where the comparison
# benchmark ran cleanly end to end.
expr_ok <- !is.na(df$status_expr_only) & df$status_expr_only == "ok" &
  !is.na(df$prepare_seconds_expr_only) &
  !is.na(df$status_full) & df$status_full == "ok"

has_tfsig <- ((!is.na(df$net_tf_edges_full)  & df$net_tf_edges_full  > 0L) |
              (!is.na(df$net_sig_edges_full) & df$net_sig_edges_full > 0L))

full_ok <- !is.na(df$status_full) & df$status_full == "ok" &
  !is.na(df$prepare_seconds_full) & has_tfsig

both_ok <- expr_ok & full_ok

df_expr <- df[expr_ok, , drop = FALSE]
df_cmp  <- df[both_ok, , drop = FALSE]

# Sort by cell count so the tables read smallest -> largest. Use the
# coalesced n_cells_any so studies with only an expression-only run
# (2322) sort into the right place instead of getting dumped at NA.
df_expr <- df_expr[order(as.numeric(df_expr$n_cells_any)), , drop = FALSE]
df_cmp  <- df_cmp[order(as.numeric(df_cmp$n_cells_any)),  , drop = FALSE]

message(sprintf("  expression-only rows: %d", nrow(df_expr)))
message(sprintf("  paired delta rows:    %d", nrow(df_cmp)))

# ---- Formatters + Markdown helpers -----------------------------------------

fmt_n    <- function(x, big = ",") formatC(x, format = "d", big.mark = big)
fmt_secs <- function(x, d = 1) {
  ifelse(is.na(x), "-",
         formatC(x, format = "f", digits = d, big.mark = ","))
}
fmt_mb   <- function(x, d = 1) {
  ifelse(is.na(x), "-",
         formatC(x, format = "f", digits = d, big.mark = ","))
}
fmt_ratio <- function(x, d = 2) {
  ifelse(is.na(x) | !is.finite(x), "-",
         formatC(x, format = "f", digits = d, big.mark = ","))
}

md_table <- function(df, align = NULL) {
  cols <- colnames(df)
  if (is.null(align)) {
    align <- ifelse(vapply(df, is.numeric, logical(1)), "r", "l")
  }
  sep_for <- function(a) switch(a, l = ":---", r = "---:", c = ":---:", "---")
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

write_tsv <- function(tbl, path) {
  utils::write.table(tbl, path, sep = "\t",
                     row.names = FALSE, quote = FALSE)
  message("Wrote ", path)
}

# ---- Table 1: expression-only -----------------------------------------------

t1 <- data.frame(
  studyID            = as.character(df_expr$studyID),
  n_cells            = fmt_n(as.integer(df_expr$n_cells_any)),
  n_genes            = fmt_n(as.integer(df_expr$n_genes_any)),
  `prepare s`        = fmt_secs(df_expr$prepare_seconds_expr_only, 0),
  `peak MB`          = fmt_mb  (df_expr$prepare_peak_mb_expr_only, 0),
  `load s (cold)`    = fmt_secs(df_expr$load_seconds_expr_only,    2),
  `load s (warm)`    = fmt_secs(df_expr$load_seconds_warm_expr_only, 2),
  `bundle MB`        = fmt_mb(
    as.numeric(df_expr$bundle_bytes_expr_only) / 1024^2, 1),
  check.names    = FALSE,
  stringsAsFactors = FALSE
)
write_tsv(t1, file.path(tables_dir, "tableS_expr_only.tsv"))
md1 <- c(
  sprintf(paste0("**Supplemental Table.** Expression-only benchmark ",
                  "(no activity, no networks), %d studies. `prepare s` ",
                  "is wall time of `prepare_study_from_eset()`; `peak ",
                  "MB` is R-reported peak working-set memory during ",
                  "that call. `load s (cold)` is the first ",
                  "`load_study()` after a fresh R session; `load s ",
                  "(warm)` is a second call in the same session (page ",
                  "cache pre-warmed). Studies are sorted by `n_cells`. ",
                  "For the manuscript Table 1 (full pipeline with ",
                  "activity + TF/sig networks) see ",
                  "`figure3_portal_studies.{md,tsv}`, generated by ",
                  "`paper/benchmarks/tables.R`."),
          nrow(t1)),
  "",
  md_table(t1))
writeLines(md1, file.path(tables_dir, "tableS_expr_only.md"))
message("Wrote ", file.path(tables_dir, "tableS_expr_only.md"))

# ---- Paired delta table ----------------------------------------------------

bnd_full <- as.numeric(df_cmp$bundle_bytes_full)      / 1024^2
bnd_expr <- as.numeric(df_cmp$bundle_bytes_expr_only) / 1024^2

t3 <- data.frame(
  studyID = as.character(df_cmp$studyID),
  n_cells = fmt_n(as.integer(df_cmp$n_cells_any)),
  # prepare seconds
  `prep s (expr)`   = fmt_secs(df_cmp$prepare_seconds_expr_only, 0),
  `prep s (full)`   = fmt_secs(df_cmp$prepare_seconds_full,      0),
  `prep s ratio`    = fmt_ratio(df_cmp$prepare_seconds_full /
                                 pmax(1e-9,
                                      df_cmp$prepare_seconds_expr_only),
                                 1),
  # peak memory
  `peak MB (expr)`  = fmt_mb(df_cmp$prepare_peak_mb_expr_only, 0),
  `peak MB (full)`  = fmt_mb(df_cmp$prepare_peak_mb_full,      0),
  `peak MB ratio`   = fmt_ratio(df_cmp$prepare_peak_mb_full /
                                 pmax(1e-9,
                                      df_cmp$prepare_peak_mb_expr_only),
                                 1),
  # cold load
  `load s (expr)`   = fmt_secs(df_cmp$load_seconds_expr_only, 2),
  `load s (full)`   = fmt_secs(df_cmp$load_seconds_full,      2),
  `load s ratio`    = fmt_ratio(df_cmp$load_seconds_full /
                                 pmax(1e-9, df_cmp$load_seconds_expr_only),
                                 1),
  # bundle MB
  `bundle MB (expr)` = fmt_mb(bnd_expr, 1),
  `bundle MB (full)` = fmt_mb(bnd_full, 1),
  `bundle MB ratio`  = fmt_ratio(bnd_full / pmax(1e-9, bnd_expr), 1),
  check.names      = FALSE,
  stringsAsFactors = FALSE
)
write_tsv(t3, file.path(tables_dir, "compare_delta.tsv"))
md3 <- c(
  sprintf(paste0("**Paired delta table.** With-vs-without-TF/sig comparison, %d ",
                  "studies. For each metric, the table shows the ",
                  "expression-only value, the full value, and the ratio ",
                  "`full / expression-only` (so 1.0 = TF/sig adds no ",
                  "cost; > 1 = TF/sig costs more). `prep` rows come ",
                  "from `prepare_study_from_eset()`; `load s` rows are ",
                  "cold `load_study()` latency; `bundle MB` is the ",
                  "lazy `.scminer.h5` on disk. Sorted by ",
                  "`n_cells`."), nrow(t3)),
  "",
  md_table(t3))
writeLines(md3, file.path(tables_dir, "compare_delta.md"))
message("Wrote ", file.path(tables_dir, "compare_delta.md"))

# ---- Console summary -------------------------------------------------------

message(sprintf("\nGenerated tables under %s/:", tables_dir))
message(sprintf("  tableS_expr_only (Supplemental) : %d rows", nrow(t1)))
message(sprintf("  compare_delta    (paired)       : %d rows", nrow(t3)))
