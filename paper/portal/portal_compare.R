#!/usr/bin/env Rscript
# paper/portal/portal_compare.R
#
# Join the per-study TSVs produced by portal_studies.R --mode runs into
# one wide comparison table, one row per study, with the columns the
# benchmark cares about for the with-vs-without-TFsig comparison:
#
#   prepare_seconds_full, prepare_seconds_expr_only, delta_prepare_seconds
#   prepare_peak_mb_full, prepare_peak_mb_expr_only, delta_prepare_peak_mb
#   load_seconds_full,    load_seconds_expr_only,    delta_load_seconds
#   load_seconds_warm_full, load_seconds_warm_expr_only, delta_load_seconds_warm
#   bundle_bytes_full,    bundle_bytes_expr_only,    delta_bundle_bytes
#   total_output_bytes_full, total_output_bytes_expr_only, delta_total_output_bytes
#   net_tf_edges_full,    net_sig_edges_full        (sanity check: > 0 in full)
#
# Usage:
#   Rscript paper/portal/portal_compare.R \
#       --expr-only-glob "paper/metrics/portal_studies_*_expression-only.tsv" \
#       --full-glob      "paper/metrics/portal_studies_*_full.tsv" \
#       --out            "paper/metrics/portal_studies_compare.tsv"

.libPaths("~/R_libs")

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- which(args == name)
  if (length(i) == 0L || i == length(args)) return(default)
  args[i + 1L]
}

expr_glob <- get_arg("--expr-only-glob",
                      "paper/metrics/portal_studies_*_expression-only.tsv")
full_glob <- get_arg("--full-glob",
                      "paper/metrics/portal_studies_*_full.tsv")
out_path  <- get_arg("--out",
                      "paper/metrics/portal_studies_compare.tsv")

load_mode <- function(glob, label) {
  files <- Sys.glob(glob)
  if (length(files) == 0L) {
    stop("No files matched glob for mode '", label, "': ", glob)
  }
  rows <- lapply(files, function(f) {
    dt <- tryCatch(
      data.table::fread(f, sep = "\t", header = TRUE),
      error = function(e) {
        message("  skip (read error): ", f, " -- ", conditionMessage(e))
        NULL
      })
    if (is.null(dt) || nrow(dt) == 0L) return(NULL)
    dt
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) stop("All TSVs empty for mode '", label, "'")
  dt <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!"studyID" %in% names(dt)) {
    stop("Mode '", label, "' TSVs missing studyID column")
  }
  message(sprintf("[%s] %d files, %d rows", label, length(files), nrow(dt)))
  dt
}

# Pull the metrics we want from each mode and tag column names with the
# mode label so they don't collide after the join.
columns_kept <- c(
  "studyID",
  "n_cells", "n_genes", "n_clusters",
  "prepare_seconds", "prepare_peak_mb",
  "load_seconds", "load_seconds_warm",
  "bundle_bytes", "total_output_bytes",
  "net_tf_edges", "net_sig_edges",
  "status", "note"
)

narrow <- function(dt, mode_label) {
  miss <- setdiff(columns_kept, names(dt))
  for (m in miss) dt[[m]] <- NA
  dt <- dt[, columns_kept, with = FALSE]
  # Rename every non-key column with the mode suffix.
  rename_cols <- setdiff(columns_kept, "studyID")
  data.table::setnames(dt, rename_cols,
                       paste0(rename_cols, "_", mode_label))
  dt
}

dt_full      <- narrow(load_mode(full_glob,      "full"),      "full")
dt_expr_only <- narrow(load_mode(expr_glob, "expression-only"), "expr_only")

joined <- merge(dt_full, dt_expr_only, by = "studyID", all = TRUE)

# Compute deltas (full - expr_only). Positive = TF/sig adds cost.
.delta <- function(dt, base, suffix_a = "_full", suffix_b = "_expr_only") {
  a <- paste0(base, suffix_a)
  b <- paste0(base, suffix_b)
  if (!a %in% names(dt) || !b %in% names(dt)) return(invisible())
  dt[[paste0("delta_", base)]] <-
    as.numeric(dt[[a]]) - as.numeric(dt[[b]])
}
for (b in c("prepare_seconds", "prepare_peak_mb",
            "load_seconds", "load_seconds_warm",
            "bundle_bytes", "total_output_bytes")) {
  .delta(joined, b)
}

# Order columns: identity, then per-metric (full, expr_only, delta) triples.
metric_groups <- c("prepare_seconds", "prepare_peak_mb",
                   "load_seconds", "load_seconds_warm",
                   "bundle_bytes", "total_output_bytes")
metric_cols <- unlist(lapply(metric_groups, function(b) {
  c(paste0(b, "_full"), paste0(b, "_expr_only"), paste0("delta_", b))
}))
identity_cols <- c("studyID",
                   "n_cells_full", "n_genes_full", "n_clusters_full",
                   "net_tf_edges_full", "net_sig_edges_full",
                   "status_full", "status_expr_only",
                   "note_full", "note_expr_only")
keep <- c(identity_cols, metric_cols)
keep <- keep[keep %in% names(joined)]
joined <- joined[, keep, with = FALSE]

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
data.table::fwrite(joined, out_path, sep = "\t", na = "NA")

message(sprintf("Wrote %d-row comparison TSV: %s", nrow(joined), out_path))

# Console summary --------------------------------------------------------------
disp <- joined[, .(
  studyID,
  n_cells       = n_cells_full,
  prep_s_full   = round(as.numeric(prepare_seconds_full),      1),
  prep_s_xonly  = round(as.numeric(prepare_seconds_expr_only), 1),
  d_prep_s      = round(as.numeric(delta_prepare_seconds),     1),
  peak_mb_full  = round(as.numeric(prepare_peak_mb_full),      0),
  peak_mb_xonly = round(as.numeric(prepare_peak_mb_expr_only), 0),
  d_peak_mb     = round(as.numeric(delta_prepare_peak_mb),     0),
  load_s_full   = round(as.numeric(load_seconds_full),         3),
  load_s_xonly  = round(as.numeric(load_seconds_expr_only),    3),
  d_load_s      = round(as.numeric(delta_load_seconds),        3),
  bundle_mb_full  = round(as.numeric(bundle_bytes_full)      / 1024^2, 1),
  bundle_mb_xonly = round(as.numeric(bundle_bytes_expr_only) / 1024^2, 1)
)]
print(disp, row.names = FALSE)
