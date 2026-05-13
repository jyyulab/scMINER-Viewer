"""Verify Python load_study + lazy gene_values against R-written bundles."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scminer_viewer import Study, load_study


def test_load_study_meta(fixture_bundle):
    s = load_study(fixture_bundle)
    assert isinstance(s, Study)
    assert s.meta.studyID == "9999"
    assert s.meta.studyAbbr == "Fx"
    assert s.meta.shortTitle == "Fx"
    assert s.meta.species == "Mus musculus"
    assert s.meta.coordinate == "UMAP"
    assert s.meta.bundleVersion == 2


def test_load_study_cells(fixture_bundle):
    s = load_study(fixture_bundle)
    assert s.n_cells == 25
    assert s.cells.index.name == "cellID"
    assert list(s.cells.columns) == ["cellType", "cellGroup", "coord1", "coord2"]
    assert s.cells["coord1"].dtype == np.float64
    assert all(isinstance(i, str) for i in s.cells.index[:3])


def test_load_study_clusters(fixture_bundle):
    s = load_study(fixture_bundle)
    assert s.clusters.index.name == "cellType"
    assert "count" in s.clusters.columns
    assert "color" in s.clusters.columns
    assert "label_1" in s.clusters.columns
    assert s.clusters["count"].sum() == s.n_cells
    assert all(c.startswith("#") for c in s.clusters["color"])


def test_load_study_indexes(fixture_bundle):
    s = load_study(fixture_bundle)
    assert s.n_genes == 12
    assert all(isinstance(g, str) for g in s.genes)
    # All 12 genes in fixture have shards in every matrix
    assert s.expression_index is not None
    assert s.activity_tf_index is not None
    assert s.activity_sig_index is not None
    assert set(s.expression_index) == set(s.genes)
    assert s.default_genes is not None
    assert len(s.default_genes) == 2


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


def test_gene_values_lazy_read(fixture_bundle):
    s = load_study(fixture_bundle)
    g = s.genes[0]
    vals = s.gene_values(g)
    assert vals is not None
    assert vals.shape == (s.n_cells,)
    # activity_tf = 2 * expression in the fixture
    vals_tf = s.gene_values(g, "Activity_tf")
    assert vals_tf is not None
    # Compare non-zero entries (zeros stay zero after multiply)
    valid = ~np.isnan(vals) & (vals != 0)
    np.testing.assert_allclose(vals_tf[valid], 2 * vals[valid])


def test_gene_values_caches_repeat_calls(fixture_bundle):
    s = load_study(fixture_bundle)
    g = s.genes[0]
    first = s.gene_values(g)
    second = s.gene_values(g)
    assert first is second  # cached by identity


def test_gene_values_unknown_gene_returns_none(fixture_bundle):
    s = load_study(fixture_bundle)
    assert s.gene_values("not_a_real_gene") is None


def test_gene_values_invalid_relationship_returns_none(fixture_bundle):
    s = load_study(fixture_bundle)
    g = s.genes[0]
    assert s.gene_values(g, "NotAValidKind") is None


def test_load_study_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_study(tmp_path / "does_not_exist.scminer.h5")


def test_load_study_explicit_shard_dir(fixture_bundle, tmp_path):
    # Copy the bundle (only) to a new dir; without explicit shard_dir
    # the lazy reads will fail (shards aren't there).
    import shutil as _shutil
    new_bundle = tmp_path / fixture_bundle.name
    _shutil.copyfile(fixture_bundle, new_bundle)
    s_wrong = load_study(new_bundle)  # shard_dir = tmp_path (no shards)
    assert s_wrong.gene_values(s_wrong.genes[0]) is None
    s_right = load_study(new_bundle, shard_dir=fixture_bundle.parent)
    assert s_right.gene_values(s_right.genes[0]) is not None
