#!/usr/bin/env Rscript
# Build the bundle for one study from an existing on-disk graph layout.
# Writes a small lazy-mode bundle next to its shard tree.
#
# Usage:
#   Rscript scminerViewer/inst/scripts/build_2327_bundle.R \
#     [data_dir] [out_path] [study_id]
#
# Defaults: data_dir = data, study_id = 2327,
#           out_path = <data_dir>/<study_id>/<study_id>.scminer.h5
#
# The `<studyID>/` subfolder convention matches what prepare_study_data()
# writes, so run_browser(<data_dir>) finds the bundle without extra
# configuration.

suppressPackageStartupMessages({
  library(scminerViewer)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data"
study_id <- if (length(args) >= 3) args[3] else "2327"
out_path <- if (length(args) >= 2) args[2] else {
  file.path(data_dir, study_id, paste0(study_id, ".scminer.h5"))
}

data_dir <- normalizePath(data_dir, mustWork = TRUE)

message("data_dir:  ", data_dir)
message("study_id:  ", study_id)
message("out_path:  ", out_path)

# read_graph_study autodetects the wrapped <data_dir>/<sid>/ layout if
# present, otherwise reads the flat layout directly under data_dir.
study <- read_graph_study(data_dir = data_dir, study_id = study_id)

count_or_zero <- function(v) if (is.null(v)) 0L else length(v)
message(sprintf("  cells=%d  genes=%d  clusters=%d",
                nrow(study$cells), length(study$genes),
                nrow(study$clusters)))
message(sprintf("  expression_index=%d genes", count_or_zero(study$expression_genes)))
message(sprintf("  activity_tf_index=%d genes", count_or_zero(study$activity_tf_genes)))
message(sprintf("  activity_sig_index=%d genes", count_or_zero(study$activity_sig_genes)))
message(sprintf("  network_tf=%d  network_sig=%d",
                if (is.null(study$network_tf))  0L else nrow(study$network_tf),
                if (is.null(study$network_sig)) 0L else nrow(study$network_sig)))

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

message("\nWriting bundle...")
write_bundle(
  bundle_path        = out_path,
  meta               = study$meta,
  cells              = study$cells,
  clusters           = study$clusters,
  genes              = study$genes,
  expression_genes   = study$expression_genes,
  activity_tf_genes  = study$activity_tf_genes,
  activity_sig_genes = study$activity_sig_genes,
  network_tf         = study$network_tf,
  network_sig        = study$network_sig,
  overwrite          = TRUE
)

message("\nReloading bundle for sanity check...")
# The shards live alongside the source data, not the (new) bundle path.
print(load_study(out_path, shard_dir = data_dir))
