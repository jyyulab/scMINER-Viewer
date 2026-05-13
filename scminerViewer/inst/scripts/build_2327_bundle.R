#!/usr/bin/env Rscript
# Build data-bundles/<studyID>.scminer.h5 from the graph-import layout in
# data/ plus the per-gene shard tree at <shard_dir>/<studyID>/{expression,
# activity}_files/. Missing shards / manifests are tolerated; the bundle
# will simply skip the corresponding group.
#
# Usage:
#   Rscript scminerViewer/inst/scripts/build_2327_bundle.R \
#     [data_dir] [out_path] [study_id] [shard_dir]
#
# Defaults: data, data-bundles/2327.scminer.h5, 2327, <data_dir>

suppressPackageStartupMessages({
  library(scminerViewer)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir  <- if (length(args) >= 1) args[1] else "data"
out_path  <- if (length(args) >= 2) args[2] else "data-bundles/2327.scminer.h5"
study_id  <- if (length(args) >= 3) args[3] else "2327"
shard_dir <- if (length(args) >= 4) args[4] else data_dir

data_dir  <- normalizePath(data_dir,  mustWork = TRUE)
shard_dir <- normalizePath(shard_dir, mustWork = TRUE)

message("data_dir:  ", data_dir)
message("shard_dir: ", shard_dir)
message("study_id:  ", study_id)
message("out_path:  ", out_path)

study <- read_graph_study(
  data_dir          = data_dir,
  study_id          = study_id,
  shard_dir         = shard_dir,
  load_expression   = TRUE,
  load_activity_tf  = TRUE,
  load_activity_sig = TRUE,
  verbose           = TRUE
)

mat_info <- function(m) {
  if (is.null(m)) return("(absent)")
  sprintf("%dx%d, nnz=%d", nrow(m), ncol(m), length(m@x))
}
message(sprintf("  cells=%d  genes=%d  clusters=%d",
                nrow(study$cells), length(study$genes),
                nrow(study$clusters)))
message("  expression:   ", mat_info(study$expression))
message("  activity_tf:  ", mat_info(study$activity_tf))
message("  activity_sig: ", mat_info(study$activity_sig))
message(sprintf("  network_tf=%d  network_sig=%d",
                if (is.null(study$network_tf))  0L else nrow(study$network_tf),
                if (is.null(study$network_sig)) 0L else nrow(study$network_sig)))

message("\nWriting bundle...")
write_bundle(
  bundle_path  = out_path,
  meta         = study$meta,
  cells        = study$cells,
  clusters     = study$clusters,
  genes        = study$genes,
  expression   = study$expression,
  activity_tf  = study$activity_tf,
  activity_sig = study$activity_sig,
  network_tf   = study$network_tf,
  network_sig  = study$network_sig,
  overwrite    = TRUE
)

message("\nReloading bundle for sanity check...")
print(load_study(out_path))
