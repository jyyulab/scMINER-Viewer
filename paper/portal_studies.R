#!/usr/bin/env Rscript
# paper/portal_studies.R
#
# Benchmark the full prepare_study_from_eset -> load_study -> gene_values
# pipeline against the 21 real scMINER Portal studies whose source data
# lives on HPC.
#
# Three input layouts are accepted:
#
# A. `--configs-dir <dir>` (preferred for the 21-study sweep).
#    A single folder of YAML config files (one per study). Each YAML
#    references the actual data files via absolute paths, or via paths
#    relative to its `input_root` key, or relative to the YAML's own
#    parent dir.
#
# B. `--config <yaml>` (preferred for an array task: one YAML at a time).
#    Process exactly one YAML.
#
# C. `--studies-root <dir>` (legacy). Each subfolder of <dir> holds
#    config.yaml plus its data files in the same folder, OR contains
#    bare expression.rds / activity.rds / networks.txt with column
#    names matching the 2327 sample (CellID, CellGroup, UMAP_1/2,
#    GeneSymbol).
#
# Sample YAML lives at data/input/2327/config.yaml.
#
# Each study's bundle, expression_files/, activity_files/, and graph
# files are written to a per-study subdirectory named after study.ID
# inside the resolved output root:
#   1. --output-root <dir>   CLI override (a single root for all studies)
#   2. cfg$output             from the YAML (per-study, default for HPC)
#   3. --scratch <dir>        (default: tempfile() -- ephemeral)
#
# So a YAML with `output: /research/.../scMINERViewerMetrics` and
# `study.ID: 2327` produces:
#   /research/.../scMINERViewerMetrics/2327/2327.scminer.h5
#   /research/.../scMINERViewerMetrics/2327/expression_files/...
#   /research/.../scMINERViewerMetrics/2327/activity_files/...
#   /research/.../scMINERViewerMetrics/2327/graph_files/...   (etc.)
#
# Per-study metrics written to paper/metrics/portal_studies.tsv:
#
#   Inputs (from the YAML's resolved paths):
#     expr_input_bytes    size of expression.rds on disk
#     act_input_bytes     size of activity.rds on disk (0 if absent)
#     net_input_bytes     size of networks.txt on disk (0 if absent)
#     total_input_bytes   sum of the three above
#   Outputs (under <scratch>/<studyID>/):
#     bundle_bytes        .scminer.h5 size
#     shard_bytes         total size of *.csv.gz shards
#     graph_bytes         everything else (graph files + meta TSVs)
#     total_output_bytes  every file under the per-study out dir
#   Wall time / memory / fetch:
#     prepare_seconds     wall time of prepare_study_from_eset()
#     prepare_peak_mb     R peak working-set memory during prepare
#     load_seconds        cold load_study() latency
#     fetch_median        median time for gene_values() over 25 random genes
#     fetch_max           worst case in that 25-gene sample
#
# Usage (from project root):
#   Rscript paper/portal_studies.R --studies-root data/input
#   Rscript paper/portal_studies.R --studies-root /hpc/.../studies
#   Rscript paper/portal_studies.R --studies-root <path> --only 2327

.libPaths("~/R_libs")

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(Matrix)
  library(scMINER)
  library(tidyverse)
  library(stringr)
  library(data.table)
  library(R.utils)
  library(scminerViewer)
})

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
}

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  i <- which(args == name)
  if (length(i) == 0L || i == length(args)) return(default)
  args[i + 1L]
}

configs_dir  <- get_arg("--configs-dir", NULL)
config_one   <- get_arg("--config",      NULL)
studies_root <- get_arg("--studies-root", NULL)
out_tsv      <- get_arg("--out", "paper/metrics/portal_studies.tsv")
only_csv     <- get_arg("--only", "")
# Where the produced bundle / expression_files / activity_files / graph
# files for each study should land. Resolution order (highest -> lowest):
#   1. --output-root <dir> CLI flag (forces a single root for every study)
#   2. cfg$output from the YAML (per-study, resolved relative to YAML dir)
#   3. --scratch <dir> (defaults to a tempfile -- ephemeral)
output_root  <- get_arg("--output-root", NULL)
scratch_root <- get_arg("--scratch",     tempfile("portal_bundle_"))
verbose      <- !("--quiet" %in% args)

modes_set <- sum(!is.null(configs_dir), !is.null(config_one),
                  !is.null(studies_root))
if (modes_set == 0L) {
  stop("Specify exactly one of --configs-dir <dir>, --config <yaml>, ",
       "or --studies-root <dir>.")
}
if (modes_set > 1L) {
  stop("--configs-dir, --config, and --studies-root are mutually exclusive.")
}
only_set <- if (nzchar(only_csv)) {
  trimws(strsplit(only_csv, ",", fixed = TRUE)[[1L]])
} else character(0)

dir.create(dirname(out_tsv), showWarnings = FALSE, recursive = TRUE)
dir.create(scratch_root,    showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) if (isTRUE(verbose)) message(...)

# ---- Discover study folders ------------------------------------------------

# Discovery for --configs-dir mode: scan for *.yaml / *.yml.
discover_configs <- function(dir_path) {
  files <- list.files(dir_path,
                      pattern = "\\.ya?ml$",
                      full.names = TRUE,
                      ignore.case = TRUE,
                      recursive = FALSE)
  sort(files)
}

# Discovery for --studies-root mode: each subdir is a study folder
# that holds either config.yaml or bare expression.rds.
discover_real_studies <- function(root) {
  dirs <- list.dirs(root, recursive = FALSE)
  bundleable <- vapply(dirs, function(d) {
    file.exists(file.path(d, "config.yaml")) ||
      file.exists(file.path(d, "expression.rds"))
  }, logical(1))
  dirs[bundleable]
}

# Networks file may be called network.txt or networks.txt (both seen
# in the wild).
find_networks_path <- function(study_dir) {
  cand <- file.path(study_dir,
                     c("networks.txt", "network.txt",
                       "networks.tsv", "network.tsv"))
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0L) NULL else hit[1L]
}

# Resolve a path stored in a YAML config. Absolute paths are returned
# as-is; relative paths are resolved against `base_dir` (typically
# the YAML's `input_root` if set, else the YAML's parent dir).
.resolve_path <- function(p, base_dir) {
  if (is.null(p) || !nzchar(p)) return(NULL)
  if (.Platform$OS.type == "windows" && grepl("^[A-Za-z]:", p)) return(p)
  if (startsWith(p, "/") || startsWith(p, "~")) {
    return(normalizePath(p, mustWork = FALSE))
  }
  normalizePath(file.path(base_dir, p), mustWork = FALSE)
}

# Build a normalised study spec from one YAML config path. Relative
# paths in `input.*` resolve against `cfg$input_root` if present
# (itself resolvable relative to the YAML's parent dir), else the
# YAML's own directory.
spec_from_yaml <- function(yaml_path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("yaml package required; install.packages('yaml').")
  }
  cfg <- yaml::yaml.load_file(yaml_path)
  yaml_dir <- dirname(normalizePath(yaml_path))
  input_base <- if (!is.null(cfg$input_root)) {
    .resolve_path(cfg$input_root, yaml_dir)
  } else yaml_dir
  # Trim accidental whitespace in IDs (seen in older portal exports).
  study_id <- trimws(as.character(cfg$study$ID %||%
                            tools::file_path_sans_ext(basename(yaml_path))))
  # `input.genes` flags the legacy raw-matrix + genes.csv layout
  # (KKYan studies) which prepare_study_from_eset() cannot consume.
  # Return a "skipped" spec so the run loop emits a placeholder row
  # and the array task still exits 0 (otherwise the merge job's
  # done(array) dependency would never be satisfied).
  if (!is.null(cfg$input$genes)) {
    return(list(
      source    = "skipped-legacy-genes-csv",
      yaml_path = yaml_path,
      study_id  = study_id,
      skip_reason = "input.genes (raw matrix + genes.csv layout)"
    ))
  }
  # `output:` from the YAML names a directory on HPC where the bundle
  # and the gene-shard tree should land. Relative paths resolve
  # against the YAML's parent dir (same rule as input.*).
  output_root <- .resolve_path(cfg$output, yaml_dir)
  list(
    source          = "yaml",
    yaml_path       = yaml_path,
    study_id        = study_id,
    output_root     = output_root,
    expr_path       = .resolve_path(cfg$input$expression, input_base),
    act_path        = .resolve_path(cfg$input$activity,   input_base),
    net_path        = .resolve_path(cfg$input$networks,   input_base),
    meta = list(
      studyID    = study_id,
      studyAbbr  = as.character(cfg$study$studyAbbr  %||% study_id),
      shortTitle = as.character(cfg$study$shortTitle %||% study_id),
      longTitle  = as.character(cfg$study$longTitle  %||% study_id),
      species    = as.character(cfg$species          %||% "unknown"),
      coordinate = as.character(cfg$coordinate       %||% "UMAP")
    ),
    # Defaults match scminerViewer / scMINER conventions:
    #   cellID     -> "cellID"
    #   cellType   -> "cellGroup"
    #   cellGroup  -> falls back to cellType (or its default)
    #   geneSymbol -> "geneSymbol"
    #   coordinate -> "UMAP"
    cell_id_col     = as.character(cfg$cellID     %||% "cellID"),
    cell_type_col   = as.character(cfg$cellType   %||% "cellGroup"),
    cell_group_col  = as.character(cfg$cellGroup  %||%
                                      cfg$cellType %||% "cellGroup"),
    coordinate_col  = as.character(cfg$coordinate %||% "UMAP"),
    gene_symbol_col = as.character(cfg$geneSymbol %||% "geneSymbol"),
    cluster_palette = as.character(cfg$cluster_palette %||% "npg")
  )
}

# Build a spec from a study folder. Prefers an in-folder config.yaml;
# falls back to convention-based defaults matching the 2327 sample.
spec_from_dir <- function(study_dir) {
  yaml_path <- file.path(study_dir, "config.yaml")
  if (file.exists(yaml_path)) return(spec_from_yaml(yaml_path))
  list(
    source          = "convention",
    yaml_path       = NA_character_,
    study_id        = basename(study_dir),
    output_root     = NULL,                          # falls back to --scratch
    expr_path       = file.path(study_dir, "expression.rds"),
    act_path        = {
      p <- file.path(study_dir, "activity.rds")
      if (file.exists(p)) p else NULL
    },
    net_path        = find_networks_path(study_dir),
    meta = list(
      studyID    = basename(study_dir),
      studyAbbr  = basename(study_dir),
      shortTitle = basename(study_dir),
      longTitle  = basename(study_dir),
      species    = "unknown",
      coordinate = "UMAP"
    ),
    cell_id_col     = "CellID",
    cell_type_col   = "CellGroup",
    cell_group_col  = "CellGroup",
    coordinate_col  = "UMAP",
    gene_symbol_col = "geneSymbol",
    cluster_palette = "npg"
  )
}

# Build the list of per-study specs from whichever input mode the
# user picked.
specs <- list()
if (!is.null(config_one)) {
  if (!file.exists(config_one)) {
    stop("--config does not exist: ", config_one)
  }
  specs <- list(spec_from_yaml(config_one))
} else if (!is.null(configs_dir)) {
  if (!dir.exists(configs_dir)) {
    stop("--configs-dir does not exist: ", configs_dir)
  }
  yaml_paths <- discover_configs(configs_dir)
  if (length(yaml_paths) == 0L) {
    stop("No *.yaml / *.yml files in --configs-dir ", configs_dir)
  }
  specs <- lapply(yaml_paths, spec_from_yaml)
  log_msg(sprintf("Found %d config(s) under %s",
                  length(yaml_paths), configs_dir))
} else {
  if (!dir.exists(studies_root)) {
    stop("--studies-root does not exist: ", studies_root)
  }
  study_dirs <- discover_real_studies(studies_root)
  if (length(study_dirs) == 0L) {
    stop("No study folders under ", studies_root,
         " contain config.yaml or expression.rds")
  }
  specs <- lapply(study_dirs, spec_from_dir)
  log_msg(sprintf("Found %d study folder(s) under %s",
                  length(study_dirs), studies_root))
}

# --only filters by study_id across all modes.
if (length(only_set) > 0L) {
  ids <- vapply(specs, `[[`, "", "study_id")
  specs <- specs[ids %in% only_set]
  if (length(specs) == 0L) {
    stop("--only=", only_csv, " matched no studies.")
  }
}

# ---- One-study benchmark ---------------------------------------------------

bench_real_one <- function(spec, scratch, cli_output_root = NULL) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase required; install via BiocManager::install('Biobase')")
  }
  study_id <- spec$study_id
  # Resolve the per-study output dir: CLI flag > YAML cfg$output > scratch.
  out_root <- if (!is.null(cli_output_root) && nzchar(cli_output_root)) {
    cli_output_root
  } else if (!is.null(spec$output_root) && nzchar(spec$output_root)) {
    spec$output_root
  } else {
    scratch
  }
  sub_root <- file.path(out_root, study_id)
  dir.create(sub_root, recursive = TRUE, showWarnings = FALSE)

  log_msg(sprintf("[%s] spec source: %s%s", study_id, spec$source,
                  if (spec$source == "yaml")
                    sprintf(" (%s)", spec$yaml_path) else ""))
  log_msg(sprintf("[%s] output dir : %s", study_id, sub_root))

  # Capture input file sizes up front (before any in-memory load).
  .file_bytes <- function(p) {
    if (is.null(p) || !file.exists(p)) return(0)
    as.numeric(file.info(p)$size)
  }
  expr_input_bytes <- .file_bytes(spec$expr_path)
  act_input_bytes  <- .file_bytes(spec$act_path)
  net_input_bytes  <- .file_bytes(spec$net_path)
  total_input_bytes <- expr_input_bytes + act_input_bytes + net_input_bytes

  log_msg(sprintf("[%s] readRDS %s (%s)", study_id,
                  basename(spec$expr_path),
                  format(expr_input_bytes, big.mark = ",")))
  expr_eset <- readRDS(spec$expr_path)

  # Pre-flight: R's Matrix package caps dgCMatrix at 2^31-1 nonzero
  # entries (sparseMatrix() throws "more than 2^31-1 nonzero entries").
  # Check the actual nnz when the eset is already sparse, otherwise
  # fall back to nrow * ncol (the cap a dense-to-sparse conversion
  # would hit).
  .nnz_cap <- .Machine$integer.max  # 2^31 - 1
  .nnz_or_cells <- function(eset) {
    m <- Biobase::exprs(eset)
    if (inherits(m, "sparseMatrix")) {
      list(nnz = length(m@x), dense = FALSE,
           nr = nrow(m), nc = ncol(m))
    } else {
      list(nnz = as.numeric(nrow(m)) * as.numeric(ncol(m)),
           dense = TRUE, nr = nrow(m), nc = ncol(m))
    }
  }
  .check_size <- function(eset, label) {
    sz <- .nnz_or_cells(eset)
    if (sz$nnz > .nnz_cap) {
      stop(sprintf(paste0(
        "matrix too large for dgCMatrix: %s is %d x %d (%s), ",
        "%.2g entries > 2^31-1 nnz cap"),
        label, sz$nr, sz$nc,
        if (sz$dense) "dense" else "sparse",
        sz$nnz),
        call. = FALSE)
    }
  }
  .check_size(expr_eset, "expression matrix")
  act_eset  <- if (!is.null(spec$act_path) && file.exists(spec$act_path)) {
    log_msg(sprintf("[%s] readRDS %s (%s)", study_id,
                    basename(spec$act_path),
                    format(act_input_bytes, big.mark = ",")))
    readRDS(spec$act_path)
  } else NULL
  if (!is.null(act_eset)) .check_size(act_eset, "activity matrix")
  log_msg(sprintf("[%s] networks: %s%s", study_id,
                  if (is.null(spec$net_path)) "(none)" else
                    basename(spec$net_path),
                  if (net_input_bytes > 0)
                    sprintf(" (%s)", format(net_input_bytes, big.mark = ","))
                  else ""))

  log_msg(sprintf("[%s] prepare_study_from_eset ...", study_id))
  invisible(gc(reset = TRUE, full = TRUE))
  t0 <- Sys.time()
  res <- prepare_study_from_eset(
    out_dir         = sub_root,
    expression_eset = expr_eset,
    activity_eset   = act_eset,
    networks_path   = spec$net_path,
    meta            = spec$meta,
    cell_id_col     = spec$cell_id_col,
    cell_type_col   = spec$cell_type_col,
    cell_group_col  = spec$cell_group_col,
    coordinate_col  = spec$coordinate_col,
    gene_symbol_col = spec$gene_symbol_col,
    cluster_palette = spec$cluster_palette,
    emit            = c("graph", "bundle"),
    verbose         = FALSE
  )
  prepare_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  gc_after <- gc(reset = FALSE, full = FALSE)
  prepare_peak_mb <- sum(as.numeric(gc_after[, ncol(gc_after)]))

  bundle <- res$bundle_path
  bundle_bytes <- file.info(bundle)$size

  shard_files <- list.files(res$out_dir, recursive = TRUE,
                             full.names = TRUE,
                             pattern = "\\.csv\\.gz$")
  shard_bytes <- sum(file.info(shard_files)$size, na.rm = TRUE)

  # Total output size: every regular file written under res$out_dir
  # (bundle + shard tree + graph_files/<sid>/* + meta/clusters TSVs).
  all_output <- list.files(res$out_dir, recursive = TRUE,
                            full.names = TRUE, all.files = TRUE,
                            no.. = TRUE)
  all_output <- all_output[!dir.exists(all_output)]
  total_output_bytes <- sum(file.info(all_output)$size, na.rm = TRUE)
  graph_bytes <- max(0, total_output_bytes - bundle_bytes - shard_bytes)

  # Free up the large in-memory esets before measuring cold load
  rm(expr_eset, act_eset)
  invisible(gc(verbose = FALSE))

  t1 <- Sys.time()
  s <- load_study(bundle)
  load_seconds <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  fetch_pool <- sample(s$expression_index,
                        min(25L, length(s$expression_index)))
  fetch_t <- vapply(fetch_pool, function(g) {
    tg <- Sys.time()
    invisible(gene_values(s, g, "Express_normalized"))
    as.numeric(difftime(Sys.time(), tg, units = "secs"))
  }, numeric(1))

  list(
    studyID            = study_id,
    n_cells            = nrow(s$cells),
    n_genes            = length(s$genes),
    n_clusters         = nrow(s$clusters),
    out_dir            = sub_root,
    # Inputs
    expr_input_bytes   = expr_input_bytes,
    act_input_bytes    = act_input_bytes,
    net_input_bytes    = net_input_bytes,
    total_input_bytes  = total_input_bytes,
    # Outputs
    bundle_bytes       = bundle_bytes,
    shard_bytes        = shard_bytes,
    graph_bytes        = graph_bytes,
    total_output_bytes = total_output_bytes,
    # Wall time / memory / fetch
    prepare_seconds    = prepare_seconds,
    prepare_peak_mb    = prepare_peak_mb,
    load_seconds       = load_seconds,
    fetch_median       = stats::median(fetch_t),
    fetch_mean         = mean(fetch_t),
    fetch_max          = max(fetch_t),
    n_fetched          = length(fetch_t),
    net_tf_edges       = if (!is.null(s$network_tf))  nrow(s$network_tf)  else 0L,
    net_sig_edges      = if (!is.null(s$network_sig)) nrow(s$network_sig) else 0L
  )
}

# ---- Run loop --------------------------------------------------------------

log_msg(sprintf("Benchmarking %d studies (scratch=%s)",
                length(specs), scratch_root))

placeholder_row <- function(spec, status, note = "") {
  data.frame(
    studyID            = spec$study_id,
    n_cells            = NA_integer_,
    n_genes            = NA_integer_,
    n_clusters         = NA_integer_,
    out_dir            = NA_character_,
    expr_input_bytes   = NA_real_,
    act_input_bytes    = NA_real_,
    net_input_bytes    = NA_real_,
    total_input_bytes  = NA_real_,
    bundle_bytes       = NA_real_,
    shard_bytes        = NA_real_,
    graph_bytes        = NA_real_,
    total_output_bytes = NA_real_,
    prepare_seconds    = NA_real_,
    prepare_peak_mb    = NA_real_,
    load_seconds       = NA_real_,
    fetch_median       = NA_real_,
    fetch_mean         = NA_real_,
    fetch_max          = NA_real_,
    n_fetched          = NA_integer_,
    net_tf_edges       = NA_integer_,
    net_sig_edges      = NA_integer_,
    status             = status,
    note               = note,
    stringsAsFactors   = FALSE
  )
}

rows <- list()
for (i in seq_along(specs)) {
  spec <- specs[[i]]
  log_msg(sprintf("\n=== [%d/%d] %s ===",
                  i, length(specs), spec$study_id))

  if (identical(spec$source, "skipped-legacy-genes-csv")) {
    log_msg(sprintf("  SKIPPED: %s -- %s",
                    spec$study_id, spec$skip_reason))
    rows[[length(rows) + 1L]] <- placeholder_row(spec, "skipped",
                                                  spec$skip_reason)
    next
  }

  err_obj <- NULL
  out <- tryCatch(bench_real_one(spec, scratch_root,
                                  cli_output_root = output_root),
                  error = function(err) {
                    err_obj <<- err
                    message(sprintf("  ERROR: %s", conditionMessage(err)))
                    NULL
                  })
  if (is.null(out)) {
    msg <- if (!is.null(err_obj)) conditionMessage(err_obj) else ""
    # Classify the 2^31-1 sparse-matrix overflow as its own status so
    # it shows up alongside legacy-skipped rows rather than as a generic
    # failure -- it is a data-size limit, not a bug.
    if (grepl("2\\^31|nnz cap|too large for dgCMatrix",
              msg, ignore.case = TRUE)) {
      rows[[length(rows) + 1L]] <- placeholder_row(spec,
                                                    "skipped-too-large",
                                                    msg)
    } else {
      rows[[length(rows) + 1L]] <- placeholder_row(spec, "error", msg)
    }
  } else {
    row <- as.data.frame(out, stringsAsFactors = FALSE)
    row$status <- "ok"
    row$note   <- ""
    rows[[length(rows) + 1L]] <- row
  }
}

if (length(rows) == 0L) {
  stop("No studies processed; see error messages above.")
}
# Column-align with NA fill so OK / skipped / error rows merge cleanly.
all_cols <- unique(unlist(lapply(rows, colnames), use.names = FALSE))
rows <- lapply(rows, function(r) {
  miss <- setdiff(all_cols, colnames(r))
  for (m in miss) r[[m]] <- NA
  r[, all_cols, drop = FALSE]
})
df <- do.call(rbind, rows)
utils::write.table(df, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
log_msg(sprintf(paste0(
  "\nWrote %d-row TSV to %s (%d ok, %d skipped, ",
  "%d too-large, %d error)"),
  nrow(df), out_tsv,
  sum(df$status == "ok",                na.rm = TRUE),
  sum(df$status == "skipped",           na.rm = TRUE),
  sum(df$status == "skipped-too-large", na.rm = TRUE),
  sum(df$status == "error",             na.rm = TRUE)))

# ---- Console summary -------------------------------------------------------

display <- df
display$in_mb     <- round(display$total_input_bytes  / 1024^2, 1)
display$out_mb    <- round(display$total_output_bytes / 1024^2, 1)
display$bundle_mb <- round(display$bundle_bytes       / 1024^2, 1)
display$shard_mb  <- round(display$shard_bytes        / 1024^2, 1)
display$prep_s    <- round(display$prepare_seconds, 1)
display$peak_mb   <- round(display$prepare_peak_mb, 0)
display$load_s    <- round(display$load_seconds, 3)
display$fetch_ms  <- round(display$fetch_median * 1000, 1)
print(display[, c("studyID", "n_cells", "n_genes", "n_clusters",
                   "in_mb", "bundle_mb", "shard_mb", "out_mb",
                   "prep_s", "peak_mb", "load_s", "fetch_ms")],
      row.names = FALSE)

log_msg(sprintf(
  paste0("\nTotals: input %.1f-%.1f MB, output %.1f-%.1f MB, ",
         "bundle %.1f-%.1f MB, prepare %.1f-%.1f s, peak %.0f-%.0f MB"),
  min(df$total_input_bytes)  / 1024^2, max(df$total_input_bytes)  / 1024^2,
  min(df$total_output_bytes) / 1024^2, max(df$total_output_bytes) / 1024^2,
  min(df$bundle_bytes)       / 1024^2, max(df$bundle_bytes)       / 1024^2,
  min(df$prepare_seconds),            max(df$prepare_seconds),
  min(df$prepare_peak_mb),            max(df$prepare_peak_mb)
))
