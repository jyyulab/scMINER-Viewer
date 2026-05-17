# scminer-viewer

Python-native package that **both prepares and serves** scMINER study
data:

* **`scminer_viewer.prepare`** — write `.scminer.h5` bundles and the
  on-disk shard tree from an [AnnData](https://anndata.readthedocs.io)
  input (Python equivalent of the R-side `scminerViewer::prepare_study`
  family).
* **`scminer_viewer`** core — Shiny-for-Python webui that renders
  bundles with the same six-tab layout used at
  <https://scminer.stjude.org/study/Tregs> (Cluster / Heatmap / Bubble
  / Feature / Violin / Network), plus a multi-study card-grid browser.

The bundles produced here are byte-compatible with the ones the R
package writes; the R `run_app()` / `run_browser()` and Python
`run` / `browse` are interchangeable readers of the same files.

- **Python ≥ 3.9** required. Runtime deps: `shiny`, `shinywidgets`,
  `h5py`, `numpy`, `scipy`, `pandas`, `plotly`.
- Optional: `networkx` (better network-graph layout; falls back to a
  circular layout if absent).
- Prepare extras: `anndata`, `pyyaml` — installed via
  `pip install scminer-viewer[prepare]`.

## Install

```sh
# Viewer only
pip install -e .

# Viewer + prepare module (anndata + pyyaml)
pip install -e ".[prepare]"

# Everything including test deps
pip install -e ".[dev]"
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

# shard_dir defaults to Path(bundle_path).parent
study = load_study("path/to/data/2327.scminer.h5")
print(study)
# <Study tex (2327) cells=8464 genes=9861 clusters=3
#  expression_index=9861 activity_tf_index=925 activity_sig_index=4708
#  defaults=- network_tf=yes network_sig=yes
#  shard_dir=path/to/data>

# Inspect data directly
print(study.cells.head())
print(study.clusters)

# Lazy: reads path/to/data/expression_files/2327/m/Mrpl15.csv.gz
expr_row = study.gene_values("Mrpl15", relationship="Express_normalized")
print(expr_row)  # ndarray of length n_cells, or None if missing

# Launch (blocking)
run_app("path/to/data/2327.scminer.h5", port=8000)

# Or build for embedding into a larger Shiny app / tests
app = build_app("path/to/data/2327.scminer.h5")
```

## Public API

| Symbol | Purpose |
| --- | --- |
| `load_study(bundle_path, shard_dir=None) -> Study` | Read a `.scminer.h5` into a `Study` dataclass. `shard_dir` defaults to `Path(bundle_path).parent`. |
| `Study`                                           | Dataclass with `meta`, `cells` (pandas, indexed by cellID), `clusters` (pandas, indexed by cellType), `genes` (numpy), `expression_index`/`activity_tf_index`/`activity_sig_index` (numpy or `None`), `default_genes` (numpy or `None` — auto-selected in the Shiny app's gene picker on launch; set this via `default_genes:` in the YAML on the R side, see [`scminerViewer` README](../scminerViewer/README.md#pre-selecting-genes-on-app-launch-default_genes)), `network_tf`/`network_sig` (pandas or `None`), `shard_dir`. |
| `Study.gene_values(gene, relationship)`          | Lazily read one gene's row from the shard tree; ndarray of length `n_cells` aligned to `Study.cells.index`, or `None` if missing. Cached per gene. `relationship` is one of `"Express_normalized"`, `"Activity_tf"`, `"Activity_sig"`. |
| `Study.has_gene(gene, relationship)`             | True if `gene` is in the corresponding bundle index. |
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
scminer-viewer info     <bundle>          # print a one-line Study summary
scminer-viewer run      <bundle> [...]    # launch the Shiny app
scminer-viewer browse   <root_dir> [...]  # multi-study card-grid browser
scminer-viewer list     <root_dir>        # enumerate every bundle under root
scminer-viewer prepare  <config.yaml>     # build a study from a YAML config
  --emit graph,bundle    subset of {graph,bundle} (default: both)
  --quiet                suppress per-shard progress messages
```

## Preparing a study (`scminer_viewer.prepare`)

The prepare submodule writes the same on-disk artifacts as the R
package's `prepare_study()`: a per-study folder containing the
`.scminer.h5` bundle plus the gzipped per-gene shard tree plus the
graph-import layout (Cell/, Gene/, Network_*/, study_meta/,
study_gene_*/).

```python
import anndata
import scminer_viewer

# 1) Fast path — already have an AnnData (or anndata-compatible) object
adata = anndata.read_h5ad("path/to/expression.h5ad")

scminer_viewer.prepare_study_from_anndata(
    out_dir="data",
    expression_adata=adata,
    meta={
        "studyID":    "99",
        "studyAbbr":  "demo",
        "longTitle":  "My demo study",
        "shortTitle": "Demo",
        "species":    "Mus musculus",
        "coordinate": "UMAP",
    },
)
# → writes data/99/99.scminer.h5 + shard tree + graph layout
scminer_viewer.run_app("data/99/99.scminer.h5")

# 2) YAML-driven (same schema as the R config, but `input.expression`
#    points at an .h5ad instead of an .rds)
scminer_viewer.prepare_study("config-99.yaml")

# 3) Low-level — bring already-extracted structures yourself
scminer_viewer.prepare_study_data(
    out_dir="data", meta=meta, cells=cells_df, genes=gene_list,
    expression=sparse_matrix,  # genes x cells
    default_genes=["GeneA", "GeneB"],
)
```

### Converting an scMINER `.rds` ExpressionSet to `.h5ad` first

The R-side `prepare_study()` reads Biobase ExpressionSet RDS files.
The Python prepare consumes [AnnData](https://anndata.readthedocs.io)
instead — the de-facto Python single-cell standard. If you have RDS
inputs (the scMINER pipeline's default output), convert them once on
the R side and the resulting `.h5ad` works with `prepare_study()`:

```r
# In R, with the sceasy + anndata packages installed
library(sceasy)
library(reticulate)
use_python(Sys.which("python3"))   # or a venv with anndata installed

eset <- readRDS("path/to/expression.rds")
sceasy::convertFormat(
  eset,
  from = "Biobase",
  to   = "anndata",
  outFile = "path/to/expression.h5ad",
)
```

Activity ExpressionSets convert the same way; the activity AnnData's
`var_names` must end in `_TF` / `.TF` or `_SIG` / `.SIG` so
`extract_activity()` can split them.

### YAML config schema

The Python `load_study_config()` accepts the same schema as the R
loader, with `input.expression` (and optional `input.activity`)
pointing at `.h5ad` files instead of `.rds`:

```yaml
output: data
study:
  ID:         "99"
  studyAbbr:  demo
  longTitle:  My demo study
  shortTitle: Demo
species:     Mus musculus
coordinate:  UMAP

cellID:      cellID
cellType:    cellGroup
cellGroup:   cellGroup
geneSymbol:  geneSymbol

input:
  expression: ./expression.h5ad
  activity:   ./activity.h5ad        # optional
  networks:   ./networks.txt          # optional, plain TSV

cluster_palette: npg
default_genes: [Cd8a, Pdcd1, Tox]    # auto-selected on app launch
```

### Prepare public API

| Symbol | Purpose |
| --- | --- |
| `prepare_study(config_path)` | YAML-driven; reads the referenced `.h5ad` / `.tsv` files and writes the full output. |
| `prepare_study_from_anndata(out_dir, expression_adata, ...)` | Accepts one or two `anndata.AnnData` objects. |
| `prepare_study_data(out_dir, meta, cells, genes, ...)` | Lowest-level form — already-extracted structures. |
| `write_bundle(bundle_path, meta, cells, clusters, genes, ...)` | HDF5 writer alone (no shard tree). |
| `read_graph_study(data_dir, study_id)` | Round-trip the on-disk graph layout into a dict the orchestrators consume. |
| `fill_clusters(cells, palette="npg")` | Per-cluster counts + ggsci-ported palette colors + centroid labels. |
| `load_study_config(path)` | Parse + validate the YAML; fill in defaults. |
| `parse_default_genes(x)` | Normalise a YAML value (list / string / `None`) into a deduped list. |
| `validate_default_genes(defaults, master)` | Filter `defaults` to those in the bundle's master gene list. |
| `read_networks(path)` | TSV → `{"tf": df_or_None, "sig": df_or_None}`. |
| `extract_cells/genes/expression/activity(adata, ...)` | AnnData extractors, for granular control. |

### Multi-study workflow

This Python package is **single-study** in v1 — pass the path to one
`.scminer.h5` bundle and the app serves that study. The companion R
package adds a multi-study card-grid landing page via
`scminerViewer::run_browser(root_dir)`; the same bundles work for both,
so you can use the R browser as the index and link out to either
package's single-study viewer.

Per-study bundles are co-located with their shard trees under a
`<root>/<studyID>/` subfolder:

```
data/
├── 2327/
│   ├── 2327.scminer.h5
│   ├── Cell/, Gene/, Network_*/, study_meta/, study_gene_*/
│   ├── expression_files/2327/
│   └── activity_files/2327/
├── 9999/
│   └── …
```

To launch a specific study from the multi-study root:

```sh
scminer-viewer run data/2327/2327.scminer.h5 --port 8000
```

`load_study(bundle_path)` defaults `shard_dir = Path(bundle_path).parent`,
so no extra config is needed when bundles are co-located with shards.
Pass `shard_dir=` (or `--shard-dir` on a future CLI release) when the
shards live elsewhere.

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

17 tests across `test_data.py` (load_study, Study, lazy gene_values,
caching, explicit shard_dir) and `test_plots.py` (six plot helpers).
Each test builds its fixture by shelling out to `Rscript` →
`scminerViewer::prepare_study_data()` (which emits both the graph
layout and the bundle), so the Python lazy reader is verified against
bytes produced by the R writer. The tests skip cleanly when `Rscript`
is unavailable.
