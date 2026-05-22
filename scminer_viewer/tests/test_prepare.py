"""Tests for the prepare submodule.

Covers:
* write_bundle round-trip (write then load_study reads back)
* prepare_study_data end-to-end with synthetic structures
* prepare_study_from_anndata end-to-end with a synthetic AnnData
* fill_clusters (palettes + centroids)
* parse_default_genes / validate_default_genes
* read_networks (TSV → DataFrames)
* load_study_config (YAML schema validation)
* read_graph_study round-trip (write graph → read it back)
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from scminer_viewer import load_study
from scminer_viewer.prepare import (
    fill_clusters,
    load_study_config,
    parse_default_genes,
    prepare_study_data,
    prepare_study_from_anndata,
    read_graph_study,
    read_networks,
    validate_default_genes,
    write_bundle,
)


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------


def _synthetic_inputs(n_cells: int = 60, n_genes: int = 40, seed: int = 7):
    rng = np.random.default_rng(seed)
    genes = [f"Gene{i:03d}" for i in range(n_genes)]
    cell_ids = [f"cell{i:03d}" for i in range(n_cells)]
    cell_types = rng.choice(["A", "B", "C"], size=n_cells)
    cells = pd.DataFrame({
        "cellID":    cell_ids,
        "cellType":  cell_types,
        "cellGroup": cell_types,
        "coord1":    rng.normal(size=n_cells),
        "coord2":    rng.normal(size=n_cells),
    })
    expression = rng.random((n_genes, n_cells)).astype(np.float64)
    activity_tf = rng.random((n_genes, n_cells)).astype(np.float64)
    activity_sig = rng.random((n_genes, n_cells)).astype(np.float64)
    network_tf = pd.DataFrame({
        "source":   ["Gene001"] * 3,
        "target":   ["Gene002", "Gene003", "Gene004"],
        "cellType": ["A", "B", "C"],
        "mi":       [0.1, 0.2, 0.3],
        "pearson":  [0.4, 0.5, 0.6],
        "spearman": [0.5, 0.6, 0.7],
        "rho":      [0.6, 0.7, 0.8],
        "pvalue":   [0.01, 0.02, 0.03],
    })
    meta = {
        "studyID":    "99",
        "studyAbbr":  "demo",
        "longTitle":  "Synthetic demo for prepare tests",
        "shortTitle": "Demo",
        "species":    "Mus musculus",
        "coordinate": "UMAP",
    }
    return meta, cells, genes, expression, activity_tf, activity_sig, network_tf


# ---------------------------------------------------------------------------
# write_bundle round-trip
# ---------------------------------------------------------------------------


def test_write_bundle_roundtrip(tmp_path: Path) -> None:
    meta, cells, genes, *_, _ = _synthetic_inputs()
    clusters = fill_clusters(cells)
    bundle_path = tmp_path / "demo.scminer.h5"

    write_bundle(
        bundle_path,
        meta=meta,
        cells=cells,
        clusters=clusters,
        genes=genes,
        expression_genes=genes,
        activity_tf_genes=genes[:10],
        activity_sig_genes=genes[:20],
        default_genes=["Gene001", "Gene005"],
    )
    assert bundle_path.exists()

    study = load_study(bundle_path)
    assert study.meta.studyID == "99"
    assert study.meta.studyAbbr == "demo"
    assert study.n_cells == len(cells)
    assert study.n_genes == len(genes)
    assert len(study.expression_index) == len(genes)
    assert len(study.activity_tf_index) == 10
    assert len(study.activity_sig_index) == 20
    assert list(study.default_genes) == ["Gene001", "Gene005"]


def test_write_bundle_refuses_overwrite(tmp_path: Path) -> None:
    meta, cells, genes, *_ = _synthetic_inputs()
    clusters = fill_clusters(cells)
    bundle = tmp_path / "demo.scminer.h5"
    write_bundle(bundle, meta=meta, cells=cells, clusters=clusters, genes=genes)
    with pytest.raises(FileExistsError):
        write_bundle(bundle, meta=meta, cells=cells, clusters=clusters, genes=genes)
    # Overwrite=True should succeed
    write_bundle(bundle, meta=meta, cells=cells, clusters=clusters,
                 genes=genes, overwrite=True)


def test_write_bundle_validates_network_columns(tmp_path: Path) -> None:
    meta, cells, genes, *_ = _synthetic_inputs()
    clusters = fill_clusters(cells)
    bad_net = pd.DataFrame({"source": ["x"], "target": ["y"]})  # missing cols
    with pytest.raises(ValueError, match="missing columns"):
        write_bundle(tmp_path / "b.h5", meta=meta, cells=cells,
                     clusters=clusters, genes=genes, network_tf=bad_net)


# ---------------------------------------------------------------------------
# fill_clusters
# ---------------------------------------------------------------------------


def test_fill_clusters_from_scratch() -> None:
    rng = np.random.default_rng(0)
    cells = pd.DataFrame({
        "cellID":   [f"c{i}" for i in range(30)],
        "cellType": ["A"] * 10 + ["B"] * 12 + ["C"] * 8,
        "coord1":   rng.normal(size=30),
        "coord2":   rng.normal(size=30),
    })
    clusters = fill_clusters(cells)
    assert set(clusters.columns) >= {"cellType", "count", "color",
                                      "label_1", "label_2"}
    assert sorted(clusters["count"].tolist()) == [8, 10, 12]
    assert all(c.startswith("#") and len(c) == 7 for c in clusters["color"])


def test_fill_clusters_preserves_existing_color() -> None:
    cells = pd.DataFrame({
        "cellID": ["a", "b"], "cellType": ["X", "Y"],
        "coord1": [0.0, 1.0], "coord2": [0.0, 1.0],
    })
    pre = pd.DataFrame({"cellType": ["X", "Y"], "color": ["#FF0000", "#00FF00"]})
    out = fill_clusters(cells, pre)
    assert out["color"].tolist() == ["#FF0000", "#00FF00"]


# ---------------------------------------------------------------------------
# parse_default_genes / validate_default_genes
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("inp,expected", [
    (None, None),
    ([], None),
    (["A", "B", "C"], ["A", "B", "C"]),
    ("A, B, C", ["A", "B", "C"]),
    ("A; B\nC", ["A", "B", "C"]),
    (["A", "A", "B"], ["A", "B"]),
])
def test_parse_default_genes(inp, expected) -> None:
    assert parse_default_genes(inp) == expected


def test_validate_default_genes_drops_unknown(recwarn) -> None:
    out = validate_default_genes(["A", "Bx", "C"], ["A", "B", "C"])
    assert out == ["A", "C"]
    assert any("Bx" in str(w.message) for w in recwarn)


# ---------------------------------------------------------------------------
# read_networks
# ---------------------------------------------------------------------------


def test_read_networks(tmp_path: Path) -> None:
    p = tmp_path / "nets.tsv"
    p.write_text(
        "source\ttarget\tNetworkType\tCellGroup\tmi\tpearson\tspearman\trho\tpvalue\n"
        "TF1\tG1\tTF\tA\t0.10\t0.20\t0.30\t0.40\t0.001\n"
        "TF1\tG2\tTF\tB\t0.15\t0.25\t0.35\t0.45\t0.002\n"
        "SIG1\tG3\tSIG\tA\t0.50\t0.60\t0.70\t0.80\t0.003\n",
        encoding="utf-8",
    )
    nets = read_networks(p)
    assert len(nets["tf"]) == 2
    assert len(nets["sig"]) == 1
    assert list(nets["tf"].columns) == [
        "source", "target", "cellType",
        "mi", "pearson", "spearman", "rho", "pvalue",
    ]
    assert nets["tf"]["mi"].iloc[0] == pytest.approx(0.10)


# ---------------------------------------------------------------------------
# load_study_config
# ---------------------------------------------------------------------------


def test_load_study_config_defaults(tmp_path: Path) -> None:
    cfg_path = tmp_path / "config.yml"
    cfg_path.write_text(
        "output: data\n"
        "study:\n"
        "  ID: '99'\n"
        "  studyAbbr: demo\n"
        "  longTitle: Demo title\n"
        "  shortTitle: Demo\n"
        "input:\n"
        "  expression: ./e.h5ad\n"
        "default_genes: [GeneA, GeneB]\n",
        encoding="utf-8",
    )
    cfg = load_study_config(cfg_path)
    assert cfg["output"] == "data"
    assert cfg["cellID"] == "cellID"
    assert cfg["cellType"] == "cellGroup"
    assert cfg["cellGroup"] == "cellGroup"
    assert cfg["coordinate"] == "UMAP"
    assert cfg["cluster_palette"] == "npg"
    assert cfg["default_genes"] == ["GeneA", "GeneB"]


def test_load_study_config_missing_keys(tmp_path: Path) -> None:
    p = tmp_path / "cfg.yml"
    p.write_text("output: data\nstudy: {}\ninput: {}\n", encoding="utf-8")
    with pytest.raises(ValueError, match="study.ID"):
        load_study_config(p)


# ---------------------------------------------------------------------------
# prepare_study_data end-to-end
# ---------------------------------------------------------------------------


def test_prepare_study_data_end_to_end(tmp_path: Path) -> None:
    meta, cells, genes, expression, activity_tf, activity_sig, network_tf = (
        _synthetic_inputs()
    )

    result = prepare_study_data(
        out_dir=tmp_path,
        meta=meta,
        cells=cells,
        genes=genes,
        expression=expression,
        activity_tf=activity_tf,
        activity_sig=activity_sig,
        network_tf=network_tf,
        default_genes=["Gene002", "Gene003"],
    )

    bundle_path = Path(result["bundle_path"])
    assert bundle_path.exists()
    study = load_study(bundle_path)
    assert study.n_cells == len(cells)
    assert study.n_genes == len(genes)
    assert list(study.default_genes) == ["Gene002", "Gene003"]

    # A few shards must exist with the right shape.
    exp_shard = tmp_path / "99" / "expression_files" / "99" / "g" / "Gene005.csv.gz"
    assert exp_shard.exists()
    with gzip.open(exp_shard, "rt") as fh:
        line = fh.readline().strip()
    vals = np.fromstring(line, sep=",", dtype=np.float64)
    assert vals.shape == (len(cells),)
    np.testing.assert_allclose(vals, expression[5, :], atol=1e-12)

    # Lazy reader should be able to round-trip values through the bundle.
    values = study.gene_values("Gene005", "Express_normalized")
    assert values is not None
    np.testing.assert_allclose(values, expression[5, :], atol=1e-12)


def test_prepare_study_data_graph_only(tmp_path: Path) -> None:
    meta, cells, genes, expression, *_ = _synthetic_inputs()
    result = prepare_study_data(
        out_dir=tmp_path, meta=meta, cells=cells, genes=genes,
        expression=expression, emit=("graph",),
    )
    assert result["bundle_path"] is None
    assert (tmp_path / "99" / "Study" / "99_study.tsv").exists()
    assert (tmp_path / "99" / "Gene" / "99_n_gene.tsv").exists()


# ---------------------------------------------------------------------------
# prepare_study_from_anndata
# ---------------------------------------------------------------------------


def test_prepare_study_from_anndata(tmp_path: Path) -> None:
    anndata = pytest.importorskip("anndata")
    rng = np.random.default_rng(1)
    n_cells, n_genes = 40, 20
    cell_ids = [f"c{i:03d}" for i in range(n_cells)]
    gene_ids = [f"G{i:03d}" for i in range(n_genes)]

    obs = pd.DataFrame({
        "cellGroup": rng.choice(["A", "B"], size=n_cells),
        "UMAP_1":    rng.normal(size=n_cells),
        "UMAP_2":    rng.normal(size=n_cells),
    }, index=cell_ids)
    var = pd.DataFrame({"geneSymbol": gene_ids}, index=gene_ids)
    X = rng.random((n_cells, n_genes)).astype(np.float64)
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    meta = {
        "studyID": "77", "studyAbbr": "syn", "longTitle": "x",
        "shortTitle": "x", "species": "Mus musculus", "coordinate": "UMAP",
    }
    result = prepare_study_from_anndata(
        out_dir=tmp_path, expression_adata=adata, meta=meta,
    )
    study = load_study(result["bundle_path"])
    assert study.n_cells == n_cells
    assert study.n_genes == n_genes

    # Round-trip one gene's values through the lazy shard reader.
    vals = study.gene_values("G005", "Express_normalized")
    assert vals is not None
    # AnnData layout is cells x genes; extractor transposes — so values
    # should match X[:, 5].
    np.testing.assert_allclose(vals, X[:, 5], atol=1e-12)


def test_prepare_study_from_anndata_explicit_coord_columns(
    tmp_path: Path,
) -> None:
    """Spatial-style studies don't follow the <stem>_1 / _2 convention.

    Verify that passing `coord1_col` / `coord2_col` directly picks
    arbitrary obs columns (here ``X`` and ``Y``) instead of looking for
    ``UMAP_1`` / ``UMAP_2``.
    """
    anndata = pytest.importorskip("anndata")
    rng = np.random.default_rng(7)
    n_cells, n_genes = 24, 10
    cell_ids = [f"s{i:03d}" for i in range(n_cells)]
    gene_ids = [f"G{i:03d}" for i in range(n_genes)]

    x_vals = np.linspace(-1.0, 1.0, n_cells)
    y_vals = np.linspace( 5.0, 7.5, n_cells)
    obs = pd.DataFrame({
        "cell_type": rng.choice(["spotA", "spotB"], size=n_cells),
        "X":         x_vals,
        "Y":         y_vals,
    }, index=cell_ids)
    var = pd.DataFrame({"geneSymbol": gene_ids}, index=gene_ids)
    X = rng.random((n_cells, n_genes)).astype(np.float64)
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    meta = {
        "studyID": "spatial1", "studyAbbr": "sp", "longTitle": "x",
        "shortTitle": "x", "species": "Mus musculus", "coordinate": "UMAP",
    }
    result = prepare_study_from_anndata(
        out_dir=tmp_path,
        expression_adata=adata,
        meta=meta,
        cell_type_col="cell_type",
        cell_group_col="cell_type",
        coord1_col="X",
        coord2_col="Y",
    )
    study = load_study(result["bundle_path"])
    np.testing.assert_allclose(study.cells["coord1"].to_numpy(), x_vals)
    np.testing.assert_allclose(study.cells["coord2"].to_numpy(), y_vals)


# ---------------------------------------------------------------------------
# read_graph_study round-trip
# ---------------------------------------------------------------------------


def test_reindex_rows_sparse_correctness() -> None:
    """The optimized CSR-buffer path must match the reference LIL path."""
    import scipy.sparse as sp
    from scminer_viewer.prepare._eset import _reindex_rows

    rng = np.random.default_rng(42)
    n_src, n_cells = 50, 80
    mat = sp.random(n_src, n_cells, density=0.25, format="csr",
                    dtype=np.float64, random_state=42)
    src_names = [f"S{i:02d}" for i in range(n_src)]

    # master with some src genes (mapped to non-sequential target rows),
    # some genes absent from src (must become zero rows), and a few src
    # genes absent from master (must be dropped).
    master = [f"M{i:02d}" for i in range(120)]
    mapping = {0: 5, 1: 17, 2: 100, 3: 42, 4: 0, 5: 60, 6: 119, 7: 80}
    for src_i, tgt_i in mapping.items():
        master[tgt_i] = src_names[src_i]

    out = _reindex_rows(mat, src_names, master)
    assert sp.issparse(out)
    assert out.shape == (120, 80)

    # Spot check: each kept src row lands at its mapped target row.
    for src_i, tgt_i in mapping.items():
        np.testing.assert_array_equal(
            out.getrow(tgt_i).toarray().ravel(),
            mat.getrow(src_i).toarray().ravel(),
        )
    # Spot check: an absent row is all zeros.
    absent_rows = [i for i in range(120) if i not in mapping.values()]
    for r in rng.choice(absent_rows, size=10, replace=False):
        assert out.getrow(r).nnz == 0


def test_reindex_rows_sparse_handles_no_overlap() -> None:
    """No src gene in master → empty CSR of the right shape."""
    import scipy.sparse as sp
    from scminer_viewer.prepare._eset import _reindex_rows
    mat = sp.random(10, 5, density=0.5, format="csr",
                    dtype=np.float64, random_state=0)
    out = _reindex_rows(mat, [f"S{i}" for i in range(10)], ["M0", "M1", "M2"])
    assert out.shape == (3, 5)
    assert out.nnz == 0


def test_read_graph_study_roundtrip(tmp_path: Path) -> None:
    meta, cells, genes, expression, *_ = _synthetic_inputs()
    prepare_study_data(
        out_dir=tmp_path, meta=meta, cells=cells, genes=genes,
        expression=expression, emit=("graph",),
    )
    read = read_graph_study(tmp_path, "99")
    assert read["meta"]["studyID"] == "99"
    assert read["meta"]["studyAbbr"] == "demo"
    assert len(read["cells"]) == len(cells)
    assert len(read["genes"]) == len(genes)
    assert read["expression_genes"] == genes
    # No activity in this fixture
    assert read["activity_tf_genes"] is None
    assert read["activity_sig_genes"] is None
