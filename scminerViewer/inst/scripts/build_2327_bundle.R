#!/usr/bin/env Rscript
# Build the bundle for one study by reading the graph-import layout in
# data/. The bundle stores only metadata + per-matrix gene indexes;
# expression / activity values stay in the on-disk shard tree and are
# read lazily by `gene_values()`.
#
# Usage:
#   Rscript scminerViewer/inst/scripts/build_2327_bundle.R \
#     [data_dir] [out_path] [study_id]
#
# Defaults: data, data/2327.scminer.h5, 2327
# The bundle is written into data/ by default so it is co-located with
# the shard tree (load_study auto-discovers shards via dirname(bundle)).

suppressPackageStartupMessages({
  library(scminerViewer)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data"
study_id <- if (length(args) >= 3) args[3] else "2327"
out_path <- if (length(args) >= 2) args[2] else {
  file.path(data_dir, paste0(study_id, ".scminer.h5"))
}

data_dir <- normalizePath(data_dir, mustWork = TRUE)

message("data_dir:  ", data_dir)
message("study_id:  ", study_id)
message("out_path:  ", out_path)

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
print(load_study(out_path))
