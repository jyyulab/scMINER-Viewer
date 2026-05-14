# Cluster auto-fill: colours from ggsci palettes (default NPG / Nature
# Publishing Group) and label coordinates from per-cluster cell centroids.
#
# These helpers are exposed (a) as a public `fill_clusters()` entry point
# so callers can pre-fill clusters before write_bundle / write_graph, and
# (b) used internally by prepare_study_data so the auto-generated
# study_meta.csv always carries real colors + label positions.

# Built-in NPG palette — used as a fallback when ggsci isn't installed.
.npg_colors_10 <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
)

# Map every YAML-friendly palette name to a ggsci `pal_<fn>(palette)`
# call. Names mirror the headings in the ggsci vignette
# (https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html).
.cluster_palette_specs <- list(
  npg          = list(fn = "pal_npg",          arg = "nrc"),
  aaas         = list(fn = "pal_aaas",         arg = "default"),
  lancet       = list(fn = "pal_lancet",       arg = "lanonc"),
  nejm         = list(fn = "pal_nejm",         arg = "default"),
  jama         = list(fn = "pal_jama",         arg = "default"),
  jco          = list(fn = "pal_jco",          arg = "default"),
  ucscgb       = list(fn = "pal_ucscgb",       arg = "default"),
  d3           = list(fn = "pal_d3",           arg = "category10"),
  locuszoom    = list(fn = "pal_locuszoom",    arg = "default"),
  igv          = list(fn = "pal_igv",          arg = "default"),
  uchicago     = list(fn = "pal_uchicago",     arg = "default"),
  startrek     = list(fn = "pal_startrek",     arg = "uniform"),
  tron         = list(fn = "pal_tron",         arg = "legacy"),
  futurama     = list(fn = "pal_futurama",     arg = "planetexpress"),
  rickandmorty = list(fn = "pal_rickandmorty", arg = "schwifty"),
  simpsons     = list(fn = "pal_simpsons",     arg = "springfield")
)

.cluster_palette_names <- function() names(.cluster_palette_specs)

# Generate `n` colours from a named ggsci palette. Falls back to the
# baked-in NPG palette when ggsci isn't installed.
.default_cluster_colors <- function(n, palette = "npg") {
  if (n <= 0L) return(character(0))
  palette <- tolower(as.character(palette))

  if (requireNamespace("ggsci", quietly = TRUE)) {
    spec <- .cluster_palette_specs[[palette]]
    if (is.null(spec)) {
      warning(sprintf("Unknown cluster palette '%s'; falling back to npg. ",
                       "Valid names: %s"),
              palette, paste(.cluster_palette_names(), collapse = ", "),
              call. = FALSE)
      spec <- .cluster_palette_specs[["npg"]]
    }
    pal_fn <- getExportedValue("ggsci", spec$fn)
    cols <- tryCatch(pal_fn(palette = spec$arg)(n),
                     error = function(e) NULL,
                     warning = function(w) suppressWarnings(pal_fn(palette = spec$arg)(n)))
    if (!is.null(cols) && length(cols) >= 1) {
      # ggsci returns 9-char hex with alpha (#RRGGBBAA); strip alpha for
      # plotly compatibility (plotly's hex parser wants #RRGGBB).
      cols <- substr(cols, 1L, 7L)
      # If the palette has fewer colors than n, ggsci returns NA tails;
      # recycle the non-NA colors to cover.
      if (any(is.na(cols))) {
        keep <- cols[!is.na(cols)]
        if (length(keep) == 0L) keep <- .npg_colors_10
        cols <- rep_len(keep, n)
      }
      return(cols)
    }
    warning(sprintf("ggsci::%s failed; falling back to built-in npg.",
                    spec$fn), call. = FALSE)
  }
  rep_len(.npg_colors_10, n)
}

# Per-cluster average UMAP coordinates from the cells data.frame. Returns
# a list with label_1 / label_2 vectors aligned to cell_types.
.compute_cluster_labels <- function(cells, cell_types) {
  stopifnot(all(c("cellType", "coord1", "coord2") %in% colnames(cells)))
  m1 <- aggregate(coord1 ~ cellType, data = cells,
                  FUN = mean, na.rm = TRUE)
  m2 <- aggregate(coord2 ~ cellType, data = cells,
                  FUN = mean, na.rm = TRUE)
  label_1 <- m1$coord1[match(cell_types, m1$cellType)]
  label_2 <- m2$coord2[match(cell_types, m2$cellType)]
  # Defensive: if a listed cluster has zero cells in `cells`, drop the
  # NA back to a sane numeric.
  label_1[is.na(label_1)] <- 0
  label_2[is.na(label_2)] <- 0
  list(label_1 = label_1, label_2 = label_2)
}

#' Auto-fill a clusters data.frame from a cells data.frame.
#'
#' Build (or top up) the clusters data.frame consumed by
#' [prepare_study_data()] / [write_bundle()]. Missing columns are
#' computed from the cells:
#'
#' \describe{
#'   \item{`count`}{Number of cells per `cellType` (from `table()`).}
#'   \item{`color`}{Colours from a ggsci palette (default NPG). Set
#'     `palette` to any of `r paste(names(.cluster_palette_specs), collapse = ", ")`.}
#'   \item{`label_1`, `label_2`}{Per-cluster centroid in `(coord1, coord2)`
#'     — `aggregate(coord1 ~ cellType, FUN = mean)` over the cells.}
#' }
#'
#' Existing values are preserved; only NA / empty cells are filled in.
#'
#' @param cells data.frame with `cellID`, `cellType`, `coord1`, `coord2`.
#' @param clusters Existing clusters data.frame, or `NULL` to build from
#'   scratch using `unique(cells$cellType)`.
#' @param palette ggsci palette name. See the ggsci vignette
#'   (<https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html>)
#'   for choices. Defaults to `"npg"` (Nature Publishing Group).
#'
#' @return A data.frame with columns `cellType`, `count`, `color`,
#'   `label_1`, `label_2`.
#' @export
fill_clusters <- function(cells, clusters = NULL, palette = "npg") {
  stopifnot(is.data.frame(cells))
  required <- c("cellType", "coord1", "coord2")
  miss <- setdiff(required, colnames(cells))
  if (length(miss) > 0L) {
    stop("`cells` missing columns: ", paste(miss, collapse = ", "))
  }

  cnt <- as.data.frame(table(cells$cellType), stringsAsFactors = FALSE)
  names(cnt) <- c("cellType", "count")
  cnt$count <- as.integer(cnt$count)

  if (is.null(clusters) || nrow(clusters) == 0L) {
    clusters <- data.frame(
      cellType = cnt$cellType,
      count    = cnt$count,
      stringsAsFactors = FALSE
    )
  } else {
    if (is.null(clusters$count) || all(is.na(clusters$count))) {
      idx <- match(clusters$cellType, cnt$cellType)
      clusters$count <- cnt$count[idx]
      clusters$count[is.na(clusters$count)] <- 0L
    }
  }

  # Fill colors
  if (is.null(clusters$color) ||
      any(is.na(clusters$color) | !nzchar(as.character(clusters$color)))) {
    auto <- .default_cluster_colors(nrow(clusters), palette = palette)
    if (is.null(clusters$color)) {
      clusters$color <- auto
    } else {
      bad <- is.na(clusters$color) | !nzchar(as.character(clusters$color))
      clusters$color[bad] <- auto[bad]
    }
  }

  # Fill label_1 / label_2 (centroids)
  need_labels <- (
    is.null(clusters$label_1) || is.null(clusters$label_2) ||
    any(is.na(clusters$label_1)) || any(is.na(clusters$label_2))
  )
  if (isTRUE(need_labels)) {
    labs <- .compute_cluster_labels(cells, clusters$cellType)
    if (is.null(clusters$label_1)) {
      clusters$label_1 <- labs$label_1
    } else {
      na1 <- is.na(clusters$label_1)
      clusters$label_1[na1] <- labs$label_1[na1]
    }
    if (is.null(clusters$label_2)) {
      clusters$label_2 <- labs$label_2
    } else {
      na2 <- is.na(clusters$label_2)
      clusters$label_2[na2] <- labs$label_2[na2]
    }
  }
  clusters
}
