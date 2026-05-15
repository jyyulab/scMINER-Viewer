# paper/benchmarks/methods.R
# Benchmarking helpers for the scMINER Viewer Applications Note.
# Sourced by paper/benchmarks/figures.R.

suppressPackageStartupMessages({
  library(Matrix)
  library(scminerViewer)
})

# ---- 1. Synthetic study generator -----------------------------------------
#
# Returns a list with everything prepare_study_data() needs. Density is the
# expected fraction of non-zero entries per gene; ~10% mirrors typical scRNA
# normalised counts after CPM/log transform.

make_synthetic_study <- function(n_cells, n_genes, n_clusters = 4L,
                                  density = 0.1, seed = 1L,
                                  with_activity = TRUE,
                                  with_networks = TRUE) {
  set.seed(seed)
  cellTypes <- paste0("CT_", LETTERS[seq_len(n_clusters)])
  cell_assign <- sample(cellTypes, n_cells, replace = TRUE)
  cells <- data.frame(
    cellID    = sprintf("cell_%06d", seq_len(n_cells)),
    cellType  = cell_assign,
    cellGroup = cell_assign,
    coord1    = stats::rnorm(n_cells),
    coord2    = stats::rnorm(n_cells),
    stringsAsFactors = FALSE
  )
  genes <- sprintf("Gene%05d", seq_len(n_genes))

  expression <- Matrix::rsparsematrix(
    n_genes, n_cells, density = density,
    rand.x = function(n) abs(stats::rnorm(n, mean = 1, sd = 1.2))
  )
  rownames(expression) <- genes
  colnames(expression) <- cells$cellID

  # scMINER's activity matrices use the FULL master gene list — only
  # the TF (or SIG) rows carry non-zero values; everything else is a
  # zero row. prepare_study_data enforces nrow(activity) == length(genes).
  activity_tf <- activity_sig <- NULL
  if (isTRUE(with_activity)) {
    tf_n  <- max(1L, floor(n_genes * 0.05))
    sig_n <- max(1L, floor(n_genes * 0.30))
    full <- function(n_active, scale) {
      x <- expression
      keep_rows <- seq_len(n_active)
      drop_rows <- setdiff(seq_len(n_genes), keep_rows)
      x[drop_rows, ] <- 0
      x[keep_rows, ] <- x[keep_rows, ] * scale
      methods::as(x, "CsparseMatrix")
    }
    activity_tf  <- full(tf_n,  1.5)
    activity_sig <- full(sig_n, 0.5)
    rownames(activity_tf)  <- genes
    rownames(activity_sig) <- genes
  }

  network_tf <- network_sig <- NULL
  if (isTRUE(with_networks)) {
    # Approximately N edges per gene per cluster — keep it bounded.
    n_edges <- min(50000L, max(500L, n_genes * 5L))
    network_tf <- data.frame(
      source   = sample(genes, n_edges, replace = TRUE),
      target   = sample(genes, n_edges, replace = TRUE),
      cellType = sample(cellTypes, n_edges, replace = TRUE),
      mi       = stats::runif(n_edges),
      pearson  = stats::runif(n_edges, -1, 1),
      spearman = stats::runif(n_edges, -1, 1),
      rho      = stats::runif(n_edges, -1, 1),
      pvalue   = stats::runif(n_edges),
      stringsAsFactors = FALSE
    )
    network_sig <- network_tf
  }

  list(
    n_cells = n_cells, n_genes = n_genes, n_clusters = n_clusters,
    density = density,
    cells = cells, genes = genes,
    expression = expression,
    activity_tf = activity_tf, activity_sig = activity_sig,
    network_tf = network_tf, network_sig = network_sig,
    meta = list(
      studyID    = sprintf("syn_%dc_%dg", n_cells, n_genes),
      studyAbbr  = "syn", longTitle = "Synthetic study",
      shortTitle = "Syn", species = "synthetic",
      coordinate = "UMAP"
    )
  )
}

# ---- 2. Write the study to a per-study folder under root ------------------

write_synthetic_study <- function(study, root) {
  prepare_study_data(
    out_dir      = root,
    meta         = study$meta,
    cells        = study$cells,
    clusters     = NULL,                     # auto-fill via fill_clusters
    genes        = study$genes,
    expression   = study$expression,
    activity_tf  = study$activity_tf,
    activity_sig = study$activity_sig,
    network_tf   = study$network_tf,
    network_sig  = study$network_sig,
    emit         = c("graph", "bundle"),
    verbose      = FALSE
  )
}

# ---- 3. Timing helpers (wall-clock seconds) -------------------------------

time_seconds <- function(expr) {
  t0 <- Sys.time()
  force(expr)
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

# Reproducible micro-benchmark — runs `expr` `reps` times, returns the
# vector of per-iteration seconds.
time_repeats <- function(expr_fn, reps = 5L) {
  vapply(seq_len(reps),
         function(i) time_seconds(expr_fn()),
         numeric(1))
}

# ---- 4. Bundle size + load + gene-fetch benchmarks ------------------------

bench_bundle <- function(study, root) {
  # ---- prepare_study_data: wall time + peak memory ----------------------
  # gc(reset = TRUE) zeroes the "max used" counter so we can measure
  # the peak during the prepare call. After the call, the 6th column
  # of gc()'s matrix gives the peak in Mb for both pools (Ncells +
  # Vcells); summing them is the rough peak working-set size.
  invisible(gc(reset = TRUE, full = TRUE))
  prepare_t <- time_seconds(res <- write_synthetic_study(study, root))
  gc_after  <- gc(reset = FALSE, full = FALSE)
  # Final column is "max used (Mb)" regardless of R's exact column
  # naming across versions — Ncells row + Vcells row = total peak Mb.
  prepare_peak_mb <- sum(as.numeric(gc_after[, ncol(gc_after)]))

  bundle <- res$bundle_path
  size   <- file.info(bundle)$size

  # Shard tree total size (for the "lazy" vs "eager" comparison)
  shard_dir <- file.path(res$out_dir, "expression_files")
  shard_size <- sum(file.info(
    list.files(shard_dir, recursive = TRUE, full.names = TRUE))$size,
                    na.rm = TRUE)

  # Cold load (after explicit gc)
  invisible(gc(verbose = FALSE))
  load_t <- time_seconds(s <- load_study(bundle))

  # First-gene fetch — pick 25 random genes from the expression index
  gene_pool <- sample(s$expression_index,
                       min(25L, length(s$expression_index)))
  fetch_t <- vapply(gene_pool, function(g) {
    time_seconds(gene_values(s, g, "Express_normalized"))
  }, numeric(1))

  list(
    n_cells           = study$n_cells,
    n_genes           = study$n_genes,
    bundle_bytes      = size,
    shard_bytes       = shard_size,
    prepare_seconds   = prepare_t,
    prepare_peak_mb   = prepare_peak_mb,
    load_seconds      = load_t,
    fetch_median      = stats::median(fetch_t),
    fetch_mean        = mean(fetch_t),
    fetch_max         = max(fetch_t),
    n_fetched         = length(fetch_t)
  )
}

# ---- 5. Multi-study discover_studies() scaling ----------------------------

bench_discover <- function(n_studies, cells_per = 200L, genes_per = 200L) {
  root <- tempfile("multi_")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  for (i in seq_len(n_studies)) {
    s <- make_synthetic_study(
      n_cells = cells_per, n_genes = genes_per,
      n_clusters = 3L, density = 0.15, seed = i,
      with_activity = FALSE, with_networks = FALSE
    )
    s$meta$studyID <- sprintf("S%04d", i)
    write_synthetic_study(s, root)
  }
  invisible(gc(verbose = FALSE))
  t_discover <- time_seconds(df <- discover_studies(root))
  list(
    n_studies      = n_studies,
    discover_seconds = t_discover,
    n_found        = nrow(df)
  )
}

# ---- 6. Bundle inspection on the real 2327 study --------------------------

bench_real_study <- function(bundle_path  = "data/2327/2327.scminer.h5",
                              shard_dir    = "data/example") {
  if (!file.exists(bundle_path)) return(NULL)
  invisible(gc(reset = TRUE, full = TRUE))
  load_t <- time_seconds(s <- load_study(bundle_path, shard_dir = shard_dir))
  gc_after <- gc(reset = FALSE, full = FALSE)
  load_peak_mb <- sum(as.numeric(gc_after[, ncol(gc_after)]))

  # If the shard tree isn't available, skip the fetch benchmark
  gene_pool <- if (!is.null(s$expression_index))
    sample(s$expression_index, min(50L, length(s$expression_index)))
  else character(0)
  fetch_t <- if (length(gene_pool) > 0L) {
    vapply(gene_pool, function(g) {
      time_seconds(gene_values(s, g, "Express_normalized"))
    }, numeric(1))
  } else numeric(0)

  list(
    bundle_path   = bundle_path,
    bundle_bytes  = file.info(bundle_path)$size,
    n_cells       = nrow(s$cells),
    n_genes       = length(s$genes),
    n_clusters    = nrow(s$clusters),
    load_seconds  = load_t,
    load_peak_mb  = load_peak_mb,
    fetch_median  = if (length(fetch_t)) stats::median(fetch_t) else NA_real_,
    fetch_mean    = if (length(fetch_t)) mean(fetch_t) else NA_real_,
    fetch_max     = if (length(fetch_t)) max(fetch_t) else NA_real_,
    n_fetched     = length(fetch_t),
    net_tf        = if (!is.null(s$network_tf)) nrow(s$network_tf) else 0L,
    net_sig       = if (!is.null(s$network_sig)) nrow(s$network_sig) else 0L
  )
}

# ---- 7. Pretty-print bytes ------------------------------------------------

human_bytes <- function(x) {
  units <- c("B", "KB", "MB", "GB", "TB")
  if (is.na(x) || x <= 0) return(sprintf("%s", x))
  k <- min(length(units) - 1L, floor(log(x, 1024)))
  sprintf("%.1f %s", x / 1024^k, units[k + 1L])
}
