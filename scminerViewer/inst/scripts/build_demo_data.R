#!/usr/bin/env Rscript
# Materialize the demo bundle that ships in inst/extdata/.
#
# Copies the 2327 (Tex) h5 bundle plus a curated 200-gene shard subset
# (~67 expression / ~67 SIG activity / ~66 TF activity) into
# scminerViewer/inst/extdata/. Re-run after the source data under
# data/ changes.
#
# Usage (from the repo root):
#   Rscript scminerViewer/inst/scripts/build_demo_data.R
#
# Source paths (relative to repo root):
#   data/2327/2327.scminer.h5
#   data/example/expression_files/2327/{meta.csv, <letter>/<gene>.csv.gz}
#   data/example/activity_files/2327/{meta.csv, {TF,SIG}/<letter>/<gene>.csv.gz}
#
# Destination: scminerViewer/inst/extdata/ (mirrors the source layout so
# load_study(demo_bundle_path()) finds shards alongside the bundle).

suppressPackageStartupMessages({
  library(hdf5r)
})

STUDY_ID <- "2327"
N_EXP    <- 67L
N_SIG    <- 67L
N_TF     <- 66L

# Curated T-cell / exhaustion markers — picked first so the demo shows
# biologically meaningful plots. Genes not present in a modality's index
# are skipped; remaining slots are filled from the head of that index.
CURATED <- c(
  "Cd8a", "Cd8b1", "Cd4", "Cd3e", "Cd3d", "Cd3g", "Cd2", "Cd5", "Cd7",
  "Pdcd1", "Lag3", "Havcr2", "Tigit", "Tox", "Tcf7", "Eomes", "Tbx21",
  "Ifng", "Gzmb", "Gzmk", "Prf1", "Il2", "Il7r", "Il2ra", "Il2rb",
  "Foxp3", "Ctla4", "Cd28", "Icos", "Cd44", "Sell", "Klrg1", "Cd69",
  "Cxcr3", "Cxcr5", "Ccr7", "Ccr5", "Nkg7", "Klf2", "Bach2", "Bcl6",
  "Runx3", "Stat1", "Stat3", "Stat4", "Stat5a", "Irf4", "Irf8", "Batf",
  "Nr4a1", "Nr4a2", "Nr4a3", "Ezh2", "Foxo1", "Myc", "Hif1a", "Mki67",
  "Top2a", "Lef1", "Itga4", "Itgae", "S1pr1", "S1pr5", "Bcl2", "Bax",
  "Gata3", "Rora", "Rorc", "Maf", "Ikzf2", "Tnfrsf9", "Tnfrsf4", "Tnfrsf18"
)

script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- grep("^--file=", args, value = TRUE)
  if (length(hit)) return(sub("^--file=", "", hit[1]))
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(ofile)
  stop("Could not determine script path; run via Rscript or source().")
})()
repo_root <- normalizePath(file.path(dirname(script_path),
                                     "..", "..", ".."),
                           mustWork = TRUE)
setwd(repo_root)

src_h5  <- file.path("data", STUDY_ID, paste0(STUDY_ID, ".scminer.h5"))
src_exp <- file.path("data", "example", "expression_files", STUDY_ID)
src_act <- file.path("data", "example", "activity_files",   STUDY_ID)

dst_root <- file.path("scminerViewer", "inst", "extdata")
dst_h5   <- file.path(dst_root, basename(src_h5))
dst_exp  <- file.path(dst_root, "expression_files", STUDY_ID)
dst_act  <- file.path(dst_root, "activity_files",   STUDY_ID)

stopifnot(file.exists(src_h5),
          dir.exists(src_exp),
          dir.exists(src_act))

# --- Read modality indexes from the bundle so we only pick real genes --
f <- H5File$new(src_h5, mode = "r")
exp_idx <- f[["index/expression"]]$read()
tf_idx  <- f[["index/activity_tf"]]$read()
sig_idx <- f[["index/activity_sig"]]$read()
f$close_all()

pick_genes <- function(index, n) {
  hits <- intersect(CURATED, index)
  fill <- setdiff(index, hits)
  c(hits, head(fill, max(0L, n - length(hits))))
}

genes_exp <- pick_genes(exp_idx, N_EXP)
genes_sig <- pick_genes(sig_idx, N_SIG)
genes_tf  <- pick_genes(tf_idx,  N_TF)

message(sprintf("expression: %d genes (%d curated, %d fill)",
                length(genes_exp), length(intersect(CURATED, genes_exp)),
                length(setdiff(genes_exp, CURATED))))
message(sprintf("SIG activity: %d genes (%d curated, %d fill)",
                length(genes_sig), length(intersect(CURATED, genes_sig)),
                length(setdiff(genes_sig, CURATED))))
message(sprintf("TF activity:  %d genes (%d curated, %d fill)",
                length(genes_tf), length(intersect(CURATED, genes_tf)),
                length(setdiff(genes_tf, CURATED))))

# --- Reset destination tree (so removed genes don't linger) ------------
unlink(dst_h5, force = TRUE)
unlink(dst_exp, recursive = TRUE, force = TRUE)
unlink(dst_act, recursive = TRUE, force = TRUE)
dir.create(dst_root, recursive = TRUE, showWarnings = FALSE)
dir.create(dst_exp,  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(dst_act, "SIG"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(dst_act, "TF"),  recursive = TRUE, showWarnings = FALSE)

shard_letter <- function(g) {
  ch <- tolower(substr(g, 1L, 1L))
  if (grepl("^[a-z]$", ch)) ch else "nm"
}

copy_shards <- function(genes, src_dir, dst_dir) {
  copied <- 0L; missed <- character()
  for (g in genes) {
    letter <- shard_letter(g)
    src <- file.path(src_dir, letter, paste0(g, ".csv.gz"))
    if (!file.exists(src)) { missed <- c(missed, g); next }
    dst_letter <- file.path(dst_dir, letter)
    if (!dir.exists(dst_letter)) {
      dir.create(dst_letter, recursive = TRUE, showWarnings = FALSE)
    }
    file.copy(src, file.path(dst_letter, basename(src)), overwrite = TRUE)
    copied <- copied + 1L
  }
  list(copied = copied, missed = missed)
}

message("\nCopying expression shards...")
res_exp <- copy_shards(genes_exp, src_exp, dst_exp)
message(sprintf("  copied %d, missing %d", res_exp$copied, length(res_exp$missed)))

message("Copying SIG activity shards...")
res_sig <- copy_shards(genes_sig, file.path(src_act, "SIG"),
                       file.path(dst_act, "SIG"))
message(sprintf("  copied %d, missing %d", res_sig$copied, length(res_sig$missed)))

message("Copying TF activity shards...")
res_tf <- copy_shards(genes_tf, file.path(src_act, "TF"),
                      file.path(dst_act, "TF"))
message(sprintf("  copied %d, missing %d", res_tf$copied, length(res_tf$missed)))

message("\nCopying meta.csv files...")
file.copy(file.path(src_exp, "meta.csv"),
          file.path(dst_exp, "meta.csv"), overwrite = TRUE)
file.copy(file.path(src_act, "meta.csv"),
          file.path(dst_act, "meta.csv"), overwrite = TRUE)

message("Copying h5 bundle...")
file.copy(src_h5, dst_h5, overwrite = TRUE)

# Manifest so users (and demo_bundle_path() docs) know which genes work.
manifest <- rbind(
  data.frame(modality = "expression",  gene = setdiff(genes_exp, res_exp$missed),
             stringsAsFactors = FALSE),
  data.frame(modality = "activity_sig", gene = setdiff(genes_sig, res_sig$missed),
             stringsAsFactors = FALSE),
  data.frame(modality = "activity_tf",  gene = setdiff(genes_tf,  res_tf$missed),
             stringsAsFactors = FALSE)
)
write.csv(manifest, file.path(dst_root, "demo_genes.csv"),
          row.names = FALSE)

message(sprintf("\nDone. inst/extdata/ now contains %d gene shards across %d modalities.",
                nrow(manifest), length(unique(manifest$modality))))
message("Bundle: ", normalizePath(dst_h5))
