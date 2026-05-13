synthetic_study <- function(n_cells = 50, n_genes = 30, n_clusters = 3,
                            seed = 1) {
  set.seed(seed)
  cellTypes <- paste0("CT_", LETTERS[seq_len(n_clusters)])
  cell_assign <- sample(cellTypes, n_cells, replace = TRUE)
  cells <- data.frame(
    cellID    = paste0("cell_", seq_len(n_cells)),
    cellType  = cell_assign,
    cellGroup = cell_assign,
    coord1    = rnorm(n_cells),
    coord2    = rnorm(n_cells),
    stringsAsFactors = FALSE
  )
  cnt <- as.data.frame(table(cells$cellType), stringsAsFactors = FALSE)
  clusters <- data.frame(
    cellType = cnt$Var1,
    count    = cnt$Freq,
    color    = grDevices::rainbow(nrow(cnt)),
    label_1  = rnorm(nrow(cnt)),
    label_2  = rnorm(nrow(cnt)),
    stringsAsFactors = FALSE
  )
  genes <- paste0("Gene", sprintf("%03d", seq_len(n_genes)))

  expression <- Matrix::rsparsematrix(
    n_genes, n_cells, density = 0.15,
    rand.x = function(n) abs(stats::rnorm(n))
  )
  rownames(expression) <- genes
  colnames(expression) <- cells$cellID
  activity_tf  <- expression * 1.5
  activity_sig <- expression * 0.5
  rownames(activity_tf)  <- genes
  rownames(activity_sig) <- genes
  colnames(activity_tf)  <- cells$cellID
  colnames(activity_sig) <- cells$cellID

  net <- data.frame(
    source   = sample(genes, 12, replace = TRUE),
    target   = sample(genes, 12, replace = TRUE),
    cellType = sample(cellTypes, 12, replace = TRUE),
    mi       = stats::runif(12),
    pearson  = stats::runif(12, -1, 1),
    spearman = stats::runif(12, -1, 1),
    rho      = stats::runif(12, -1, 1),
    pvalue   = stats::runif(12),
    stringsAsFactors = FALSE
  )

  list(
    meta = list(
      studyID    = "9999",
      studyAbbr  = "Test",
      longTitle  = "A synthetic test study",
      shortTitle = "Test",
      species    = "Mus musculus",
      coordinate = "UMAP"
    ),
    cells        = cells,
    clusters     = clusters,
    genes        = genes,
    expression   = expression,
    activity_tf  = activity_tf,
    activity_sig = activity_sig,
    network_tf   = net,
    network_sig  = net
  )
}
