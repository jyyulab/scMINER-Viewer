#!/usr/bin/env Rscript
# paper/portal/sparseify_eset.R
#
# Convert a Biobase::ExpressionSet whose `exprs()` slot is stored as a
# dense matrix (`dgeMatrix` or base `matrix`) into one backed by a
# `dgCMatrix`. Use this one-time on HPC for any portal study whose
# loaded eset blows past the available RAM (a typical symptom is
# `TERM_MEMLIMIT` during readRDS for 100k+ cell studies, like
# 2317 / Covid650k).
#
# Once the sparse rds is written, update the corresponding YAML's
# `input.expression` (and/or `input.activity`) to point at the new
# path; subsequent `prepare_study_from_eset()` runs use ~5-10x less
# memory and finish proportionally faster.
#
# Usage:
#   Rscript paper/portal/sparseify_eset.R \
#       --in  /research/.../Covid650k/expression.rds \
#       --out /research/.../Covid650k/expression.sparse.rds
#
#   # Or convert in place (writes <stem>.sparse.rds next to the source):
#   Rscript paper/portal/sparseify_eset.R \
#       --in  /research/.../Covid650k/expression.rds
#
#   # Verify only (no write); reports class + dims + nnz of the source:
#   Rscript paper/portal/sparseify_eset.R \
#       --in  /research/.../Covid650k/expression.rds \
#       --verify-only
#
# Flags:
#   --in  <path>      (required)  source .rds (Biobase ExpressionSet)
#   --out <path>      (optional)  destination .rds; defaults to
#                                  <stem>.sparse.rds in source dir
#   --force           overwrite --out if it exists (default: skip)
#   --verify-only     load + report only; do not write
#   --quiet           suppress progress messages

.libPaths("~/R_libs")

suppressPackageStartupMessages({
  library(Biobase)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- which(args == name)
  if (length(i) == 0L || i == length(args)) return(default)
  args[i + 1L]
}
src_path     <- get_arg("--in",  NULL)
dst_path     <- get_arg("--out", NULL)
force_over   <- "--force"       %in% args
verify_only  <- "--verify-only" %in% args
verbose      <- !("--quiet" %in% args)

if (is.null(src_path)) {
  stop("Missing --in <path>; pass the source .rds to convert.")
}
if (!file.exists(src_path)) {
  stop("--in does not exist: ", src_path)
}

if (is.null(dst_path) && !verify_only) {
  stem <- tools::file_path_sans_ext(src_path)
  dst_path <- paste0(stem, ".sparse.rds")
}

log_msg <- function(...) if (isTRUE(verbose)) message(...)

# ---- Load source -----------------------------------------------------------

t_load <- Sys.time()
log_msg(sprintf("Reading %s (%s bytes on disk) ...",
                src_path,
                format(file.info(src_path)$size, big.mark = ",")))

eset <- readRDS(src_path)
if (!is(eset, "ExpressionSet")) {
  stop("Source is not an ExpressionSet (got class ", class(eset)[1L], ")")
}

m <- Biobase::exprs(eset)
nr <- nrow(m); nc <- ncol(m)
src_class <- class(m)[1L]
already_sparse <- inherits(m, "sparseMatrix")

if (already_sparse) {
  src_nnz <- length(m@x)
} else {
  src_nnz <- sum(as.numeric(m != 0), na.rm = TRUE)
}
load_secs <- as.numeric(difftime(Sys.time(), t_load, units = "secs"))

log_msg(sprintf(
  "  load:     %.1f s",                load_secs))
log_msg(sprintf(
  "  class:    %s",                    src_class))
log_msg(sprintf(
  "  dims:     %d genes x %d cells",   nr, nc))
log_msg(sprintf(
  "  nnz:      %s  (density %.2f%%)",
  format(src_nnz, big.mark = ","),
  100 * src_nnz / (as.numeric(nr) * nc)))

if (already_sparse) {
  log_msg("  status:   exprs() is already sparse -- nothing to convert.")
  if (!verify_only && !is.null(dst_path) &&
      normalizePath(src_path, mustWork = FALSE) !=
      normalizePath(dst_path, mustWork = FALSE)) {
    log_msg("  copying source to ", dst_path, " unchanged ...")
    if (file.exists(dst_path) && !force_over) {
      stop("Destination exists; pass --force to overwrite: ", dst_path)
    }
    file.copy(src_path, dst_path, overwrite = TRUE)
    log_msg("Done.")
  }
  quit(status = 0)
}

# Cap check before allocating ----------------------------------------------
.nnz_cap <- .Machine$integer.max  # 2^31 - 1
if (src_nnz > .nnz_cap) {
  stop(sprintf(paste0(
    "Source nnz=%.3g exceeds 2^31-1 dgCMatrix cap. dgCMatrix cannot ",
    "represent this matrix; consider a different storage backend ",
    "(HDF5Array, BPCells) or downsampling."),
    src_nnz))
}

if (isTRUE(verify_only)) {
  log_msg("  --verify-only set; not writing.")
  quit(status = 0)
}

if (is.null(dst_path)) {
  stop("--out resolved to NULL; bug?")
}

if (file.exists(dst_path) && !force_over) {
  stop("Destination already exists (pass --force to overwrite): ",
       dst_path)
}

# ---- Convert + write -------------------------------------------------------

log_msg(sprintf("Converting %s -> dgCMatrix ...", src_class))
t_conv <- Sys.time()
# Biobase only registers `exprs<-` for signature (ExpressionSet, matrix);
# a dgCMatrix isn't a base matrix, so dispatch fails for every ExpressionSet
# subclass (including scMINER's SparseExpressionSet). Write to assayData
# directly to bypass the restriction; dimnames carry over from `m`.
sp <- methods::as(m, "CsparseMatrix")
rm(m); invisible(gc(verbose = FALSE))
eset <- Biobase::assayDataElementReplace(eset, "exprs", sp, validate = FALSE)
rm(sp); invisible(gc(verbose = FALSE))
conv_secs <- as.numeric(difftime(Sys.time(), t_conv, units = "secs"))
log_msg(sprintf("  convert:  %.1f s", conv_secs))

new_m <- Biobase::exprs(eset)
log_msg(sprintf("  new class: %s   nnz=%s",
                class(new_m)[1L],
                format(length(new_m@x), big.mark = ",")))

log_msg(sprintf("Writing %s ...", dst_path))
t_save <- Sys.time()
saveRDS(eset, dst_path, compress = TRUE)
save_secs <- as.numeric(difftime(Sys.time(), t_save, units = "secs"))
log_msg(sprintf("  save:     %.1f s", save_secs))

src_sz <- file.info(src_path)$size
dst_sz <- file.info(dst_path)$size
log_msg(sprintf(
  "Done. Source %s bytes  ->  sparse %s bytes  (%.1fx %s)",
  format(src_sz, big.mark = ","),
  format(dst_sz, big.mark = ","),
  abs(dst_sz / src_sz),
  if (dst_sz < src_sz) "smaller" else "bigger"))

log_msg("Next: update the relevant YAML(s) under paper/configs/ to point")
log_msg(sprintf("      input.expression at %s", dst_path))
