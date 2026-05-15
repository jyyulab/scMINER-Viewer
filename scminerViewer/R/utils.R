`%||%` <- function(a, b) if (is.null(a)) b else a

.bundle_version <- 1L

# Shared by graph_read.R (when finding shards) and write_graph.R (when
# emitting them). Lowercase first letter, or "nm" for non-alphabetic.
.shard_letter <- function(gene) {
  first <- tolower(substr(gene, 1, 1))
  if (grepl("^[a-z]$", first)) first else "nm"
}
