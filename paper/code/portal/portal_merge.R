#!/usr/bin/env Rscript
# paper/code/portal/portal_merge.R
#
# Concatenate the per-study TSVs that the bsub job array produced into
# the final paper/metrics/portal_studies/portal_studies.tsv. Called by the merge job
# scheduled in paper/code/portal/portal_studies_hpc.sh.
#
# Usage:
#   Rscript paper/code/portal/portal_merge.R <out.tsv> <in1.tsv> <in2.tsv> ...

.libPaths("~/R_libs")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: portal_merge.R <out.tsv> <in1.tsv> [<in2.tsv> ...]")
}
out_path <- args[1L]
in_paths <- Sys.glob(args[-1L])
in_paths <- unique(in_paths[file.exists(in_paths)])

if (length(in_paths) == 0L) {
  stop("No input TSVs matched: ", paste(args[-1L], collapse = ", "))
}

frames <- lapply(in_paths, function(p) {
  utils::read.table(p, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE, check.names = FALSE)
})

# Column-align via union; missing columns get NA.
all_cols <- unique(unlist(lapply(frames, colnames)))
frames <- lapply(frames, function(df) {
  miss <- setdiff(all_cols, colnames(df))
  for (m in miss) df[[m]] <- NA
  df[, all_cols, drop = FALSE]
})
merged <- do.call(rbind, frames)

# Stable per-study order
if ("studyID" %in% colnames(merged)) {
  merged <- merged[order(merged$studyID), , drop = FALSE]
}

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
utils::write.table(merged, out_path, sep = "\t",
                   row.names = FALSE, quote = FALSE)
message(sprintf("Wrote %d-row TSV to %s (from %d inputs)",
                nrow(merged), out_path, length(in_paths)))
