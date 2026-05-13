"""Shared pytest fixtures for scminer_viewer.

We rely on R + the scminerViewer package to produce a synthetic data
tree (graph layout + shards) plus a bundle. This way the Python reader
is genuinely tested against bytes written by the R writer.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _rscript_available() -> bool:
    return shutil.which("Rscript") is not None


@pytest.fixture(scope="session")
def fixture_bundle(tmp_path_factory) -> Path:
    """Build a synthetic study layout + bundle via the R package.

    Returns the path to the bundle. The bundle's parent directory is the
    shard_dir (so `load_study(bundle)` auto-discovers the shard tree).
    """
    if not _rscript_available():
        pytest.skip("Rscript not on PATH; cannot build fixture bundle.")
    out_dir = tmp_path_factory.mktemp("study")
    out_dir_str = repr(str(out_dir))
    rscript = f"""
    suppressPackageStartupMessages(library(scminerViewer))
    set.seed(7)
    n_cells <- 25; n_genes <- 12
    cellTypes <- c("A", "B", "C")
    cells <- data.frame(
      cellID    = paste0("c", seq_len(n_cells)),
      cellType  = sample(cellTypes, n_cells, replace = TRUE),
      cellGroup = sample(cellTypes, n_cells, replace = TRUE),
      coord1    = rnorm(n_cells),
      coord2    = rnorm(n_cells),
      stringsAsFactors = FALSE
    )
    cnt <- as.data.frame(table(cells$cellType), stringsAsFactors = FALSE)
    clusters <- data.frame(
      cellType = cnt$Var1, count = cnt$Freq,
      color    = c("#aa0000", "#00aa00", "#0000aa")[seq_len(nrow(cnt))],
      label_1  = rnorm(nrow(cnt)),
      label_2  = rnorm(nrow(cnt)),
      stringsAsFactors = FALSE
    )
    genes <- paste0("g", sprintf("%03d", seq_len(n_genes)))
    M <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.25,
                                rand.x = function(n) abs(rnorm(n)))
    rownames(M) <- genes
    colnames(M) <- cells$cellID
    net <- data.frame(
      source   = sample(genes, 8, replace = TRUE),
      target   = sample(genes, 8, replace = TRUE),
      cellType = sample(cellTypes, 8, replace = TRUE),
      mi       = runif(8),
      pearson  = runif(8, -1, 1),
      spearman = runif(8, -1, 1),
      rho      = runif(8, -1, 1),
      pvalue   = runif(8),
      stringsAsFactors = FALSE
    )
    res <- prepare_study_data(
      out_dir = {out_dir_str},
      meta = list(
        studyID    = "9999", studyAbbr = "Fx",
        longTitle  = "Fixture study", shortTitle = "Fx",
        species    = "Mus musculus", coordinate = "UMAP"
      ),
      cells        = cells,
      clusters     = clusters,
      genes        = genes,
      expression   = M,
      activity_tf  = M * 2,
      activity_sig = M * 3,
      network_tf   = net,
      network_sig  = net,
      default_genes = genes[1:2],
      emit          = c("graph", "bundle")
    )
    cat(res$bundle_path)
    """
    res = subprocess.run(
        ["Rscript", "-e", rscript],
        capture_output=True, text=True, cwd=PROJECT_ROOT,
    )
    if res.returncode != 0:
        pytest.skip(f"Could not build fixture bundle via R:\n{res.stderr}")
    bundle_path = Path(res.stdout.strip())
    assert bundle_path.exists(), f"Bundle not produced at {bundle_path}"
    return bundle_path
