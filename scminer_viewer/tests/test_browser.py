"""Tests for the Python multi-study browser API."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from scminer_viewer import discover_studies, build_browser


def _rscript_available() -> bool:
    return shutil.which("Rscript") is not None


PROJECT_ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def two_study_root(tmp_path_factory) -> Path:
    """Build a root with two distinct synthetic studies via R."""
    if not _rscript_available():
        pytest.skip("Rscript not on PATH; cannot build fixture studies.")
    root = tmp_path_factory.mktemp("multistudies")
    root_str = repr(str(root))
    rscript = f"""
    suppressPackageStartupMessages(library(scminerViewer))
    make_one <- function(sid, abbr, seed) {{
      set.seed(seed)
      n_cells <- 20; n_genes <- 10
      cellTypes <- c("A","B")
      cells <- data.frame(
        cellID = paste0("c", seq_len(n_cells)),
        cellType = sample(cellTypes, n_cells, replace = TRUE),
        cellGroup = sample(cellTypes, n_cells, replace = TRUE),
        coord1 = rnorm(n_cells), coord2 = rnorm(n_cells),
        stringsAsFactors = FALSE
      )
      cnt <- as.data.frame(table(cells$cellType), stringsAsFactors=FALSE)
      clusters <- data.frame(
        cellType = cnt$Var1, count = cnt$Freq,
        color = c("#aa0000","#0000aa")[seq_len(nrow(cnt))],
        stringsAsFactors = FALSE
      )
      genes <- paste0("g", sprintf("%03d", seq_len(n_genes)))
      M <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3,
                                  rand.x = function(n) abs(rnorm(n)))
      rownames(M) <- genes; colnames(M) <- cells$cellID
      net <- data.frame(
        source = sample(genes, 5, replace=TRUE),
        target = sample(genes, 5, replace=TRUE),
        cellType = sample(cellTypes, 5, replace=TRUE),
        mi = runif(5), pearson = runif(5,-1,1),
        spearman = runif(5,-1,1), rho = runif(5,-1,1),
        pvalue = runif(5), stringsAsFactors=FALSE
      )
      prepare_study_data(
        out_dir = {root_str},
        meta = list(studyID = sid, studyAbbr = abbr,
                    longTitle = paste("Study", sid),
                    shortTitle = abbr,
                    species = "Mus musculus", coordinate = "UMAP"),
        cells = cells, clusters = clusters, genes = genes,
        expression = M, activity_tf = M*2, activity_sig = M*3,
        network_tf = net, network_sig = net,
        default_genes = genes[1:2],
        emit = c("graph", "bundle")
      )
    }}
    make_one("1111", "one", 1)
    make_one("2222", "two", 2)
    """
    res = subprocess.run(
        ["Rscript", "-e", rscript],
        capture_output=True, text=True, cwd=PROJECT_ROOT,
    )
    if res.returncode != 0:
        pytest.skip(f"Could not build fixture studies via R:\n{res.stderr}")
    return root


def test_discover_studies_empty_dir(tmp_path):
    empty = tmp_path / "empty"
    df = discover_studies(empty)
    assert isinstance(df, pd.DataFrame)
    assert df.empty
    assert set(df.columns) >= {
        "studyID", "studyAbbr", "shortTitle", "longTitle",
        "species", "n_cells", "n_genes", "n_clusters",
        "bundle_path", "study_dir",
    }


def test_discover_studies_finds_all_bundles(two_study_root):
    df = discover_studies(two_study_root)
    assert len(df) == 2
    assert set(df["studyID"]) == {"1111", "2222"}
    assert set(df["studyAbbr"]) == {"one", "two"}
    assert (df["n_cells"] > 0).all()
    assert all(Path(p).exists() for p in df["bundle_path"])


def test_discover_studies_sorted_by_id(two_study_root):
    df = discover_studies(two_study_root)
    assert list(df["studyID"]) == sorted(df["studyID"])


def test_build_browser_constructs_app(two_study_root):
    from shiny import App
    app = build_browser(two_study_root)
    assert isinstance(app, App)


def test_build_browser_with_shard_dir(two_study_root):
    from shiny import App
    # In the typical case the shards are already co-located; here we
    # simply verify the parameter is plumbed through.
    app = build_browser(two_study_root, shard_dir=two_study_root)
    assert isinstance(app, App)
