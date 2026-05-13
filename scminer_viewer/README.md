# scminer-viewer

Shiny-for-Python webui that loads `.scminer.h5` study bundles produced
by the companion R package [`scminerViewer`](../scminerViewer/) and
renders them with the same six-tab layout used at
<https://scminer.stjude.org/study/Tregs> (Cluster / Heatmap / Bubble /
Feature / Violin / Network).

- **Python ≥ 3.9** required. Runtime deps: `shiny`, `shinywidgets`,
  `h5py`, `numpy`, `scipy`, `pandas`, `plotly`.
- Optional: `networkx` (better network-graph layout; falls back to a
  circular layout if absent).

## Install

```sh
# From a checkout, editable, with pytest extras
pip install -e ".[dev]"

# Or from a built wheel
pip install scminer_viewer-0.1.0-py3-none-any.whl
```

## Quick start

```sh
scminer-viewer info path/to/2327.scminer.h5
scminer-viewer run  path/to/2327.scminer.h5 --port 8000
# add --no-browser to skip opening a browser window
```

Programmatic usage:

```python
from scminer_viewer import load_study, build_app, run_app

study = load_study("path/to/2327.scminer.h5")
print(study)
# <Study tex (2327) cells=8464 genes=9861 clusters=3
#  expression=- activity_tf=- activity_sig=- network_tf=yes network_sig=yes>

# Inspect data directly
print(study.cells.head())
print(study.clusters)
expr_row = study.gene_values("Mrpl15", relationship="Express_normalized")
print(expr_row)  # ndarray of length n_cells, or None if the matrix is absent

# Launch (blocking)
run_app("path/to/2327.scminer.h5", port=8000)

# Or build for embedding into a larger Shiny app / tests
app = build_app("path/to/2327.scminer.h5")
```

## Public API

| Symbol | Purpose |
| --- | --- |
| `load_study(bundle_path) -> Study`               | Read a `.scminer.h5` into a `Study` dataclass. |
| `Study`                                           | Dataclass with `meta`, `cells` (pandas, indexed by cellID), `clusters` (pandas, indexed by cellType), `genes` (numpy), `expression`/`activity_tf`/`activity_sig` (scipy CSR or `None`), `network_tf`/`network_sig` (pandas or `None`). |
| `Study.gene_values(gene, relationship)`          | Dense ndarray row for one gene, or `None` if missing. `relationship` is one of `"Express_normalized"`, `"Activity_tf"`, `"Activity_sig"`. |
| `build_app(bundle_path) -> shiny.App`            | Build a Shiny app without launching it. |
| `run_app(bundle_path, host, port, launch_browser)` | Launch + serve the Shiny app via uvicorn (blocking). |
| `plots.cluster_plot(study, ...)`                 | Plotly cluster scatter (coord1, coord2 colored by cellType). |
| `plots.feature_plot(study, gene, ...)`           | Per-cell scatter colored by gene values. |
| `plots.violin_plot(study, gene, ...)`            | Per-cluster violin of gene values. |
| `plots.heatmap_plot(study, genes, ...)`          | Gene × cluster mean heatmap. |
| `plots.bubble_plot(study, genes, ...)`           | Gene × cluster, dot size = pct expressing, color = mean. |
| `plots.network_plot(study, gene, network_type, ...)` | TF/SIG neighbourhood of one gene rendered with `networkx` spring layout. |

CLI:

```text
scminer-viewer info <bundle>            # print a one-line Study summary
scminer-viewer run  <bundle> [options]  # launch the Shiny app
  --host HOST          default 127.0.0.1
  --port PORT          default 8000
  --no-browser         don't open a browser window
```

## Bundle format

The HDF5 bundle layout is the R↔Python contract described in
[`../IMPLEMENTATION.md`](../IMPLEMENTATION.md). Strings are UTF-8;
sparse matrices follow the scipy CSR convention (zero-based indices,
contiguous `data`/`indices`/`indptr`/`shape`). Scalar string datasets
written by R appear as 1-element arrays here — the reader unwraps them
automatically.

## Tests

```sh
pytest tests -q
```

14 tests across `test_data.py` (load_study, Study, gene_values) and
`test_plots.py` (six plot helpers). Each test builds its fixture by
shelling out to `Rscript` → `scminerViewer::write_bundle()`, so the
Python reader is verified against bytes produced by the R writer. The
tests skip cleanly when `Rscript` is unavailable.
