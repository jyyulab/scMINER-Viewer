"""Verify Python load_study reads bundles written by R round-trips correctly."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from scminer_viewer import Study, load_study


def test_load_study_meta(fixture_bundle):
    s = load_study(fixture_bundle)
    assert isinstance(s, Study)
    assert s.meta.studyID == "9999"
    assert s.meta.studyAbbr == "Fx"
    assert s.meta.shortTitle == "Fx"
    assert s.meta.species == "Mus musculus"
    assert s.meta.coordinate == "UMAP"
    assert s.meta.bundleVersion == 1


def test_load_study_cells(fixture_bundle):
    s = load_study(fixture_bundle)
    assert s.n_cells == 25
    assert s.cells.index.name == "cellID"
    assert list(s.cells.columns) == ["cellType", "cellGroup", "coord1", "coord2"]
    assert s.cells["coord1"].dtype == np.float64
    # IDs should round-trip as Python strings (not bytes)
    assert all(isinstance(i, str) for i in s.cells.index[:3])


def test_load_study_clusters(fixture_bundle):
    s = load_study(fixture_bundle)
    assert s.clusters.index.name == "cellType"
    assert "count" in s.clusters.columns
    assert "color" in s.clusters.columns
    assert "label_1" in s.clusters.columns
    assert s.clusters["count"].sum() == s.n_cells
    assert all(c.startswith("#") for c in s.clusters["color"])


def test_load_study_genes_and_matrices(fixture_bundle):
    s = load_study(fixture_bundle)
    assert s.n_genes == 12
    assert all(isinstance(g, str) for g in s.genes)
    assert isinstance(s.expression, csr_matrix)
    assert s.expression.shape == (12, 25)
    assert isinstance(s.activity_tf, csr_matrix)
    assert isinstance(s.activity_sig, csr_matrix)
    # activity_tf = 2 * expression in the fixture
    diff = np.abs(s.activity_tf.toarray() - 2 * s.expression.toarray()).max()
    assert diff < 1e-9
    # activity_sig = 3 * expression
    diff = np.abs(s.activity_sig.toarray() - 3 * s.expression.toarray()).max()
    assert diff < 1e-9


def test_load_study_networks(fixture_bundle):
    s = load_study(fixture_bundle)
    for net in (s.network_tf, s.network_sig):
        assert isinstance(net, pd.DataFrame)
        assert set(["source", "target", "cellType",
                    "mi", "pearson", "spearman", "rho", "pvalue"]).issubset(
            net.columns
        )
        assert net["mi"].dtype == np.float64
        assert len(net) == 8


def test_gene_values_returns_dense_row(fixture_bundle):
    s = load_study(fixture_bundle)
    g = s.genes[0]
    vals = s.gene_values(g)
    assert vals is not None
    assert vals.shape == (s.n_cells,)
    # Match the underlying CSR row
    expected = s.expression.getrow(0).toarray().ravel()
    np.testing.assert_allclose(vals, expected)


def test_gene_values_unknown_gene_returns_none(fixture_bundle):
    s = load_study(fixture_bundle)
    assert s.gene_values("not_a_real_gene") is None


def test_load_study_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_study(tmp_path / "does_not_exist.scminer.h5")
