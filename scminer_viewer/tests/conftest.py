"""Shared pytest fixtures for scminer_viewer.

We rely on R + the scminerViewer package to produce synthetic .scminer.h5
fixtures (so the Python reader is genuinely tested against bytes written
by the R writer, not by code in this repo).
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _rscript_available() -> bool:
    return shutil.which("Rscript") is not None


@pytest.fixture(scope="session")
def fixture_bundle(tmp_path_factory) -> Path:
    """Build a small .scminer.h5 via the R scminerViewer package."""
    if not _rscript_available():
        pytest.skip("Rscript not on PATH; cannot build fixture bundle.")
    out_dir = tmp_path_factory.mktemp("bundle")
    bundle_path = out_dir / "fixture.scminer.h5"
    bundle_path_str = repr(str(bundle_path))
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
    write_bundle(
      bundle_path  = {bundle_path_str},
      meta         = list(
        studyID    = "9999", studyAbbr = "Fx",
        longTitle  = "Fixture study", shortTitle = "Fx",
        species    = "Mus musculus", coordinate = "UMAP"
      ),
      cells        = cells, clusters = clusters, genes = genes,
      expression   = M, activity_tf  = M * 2, activity_sig = M * 3,
      network_tf   = net, network_sig = net,
      overwrite    = TRUE
    )
    """
    res = subprocess.run(
        ["Rscript", "-e", rscript],
        capture_output=True, text=True, cwd=PROJECT_ROOT,
    )
    if res.returncode != 0:
        pytest.skip(f"Could not build fixture bundle via R:\n{res.stderr}")
    assert bundle_path.exists(), f"Bundle not produced at {bundle_path}"
    return bundle_path
