# scMINER Viewer — Architecture & Status

Two sibling packages for viewing scMINER studies offline, both file-based
and standalone (no Java services, no graph DB, no MySQL):

| Package | Language | Version | Path | Role |
| --- | --- | --- | --- | --- |
| `scminerViewer`  | R      | 0.1.0 | [`scminerViewer/`](scminerViewer/)   | Prepares scMINER studies (graph-import layout + `.scminer.h5` bundle); serves a single-study Shiny app and a multi-study browser. |
| `scminer_viewer` | Python | 0.2.0 | [`scminer_viewer/`](scminer_viewer/) | Loads the same HDF5 bundle and serves a Shiny-for-Python single-study viewer and multi-study browser. |

Deliverable: [`data/2327/2327.scminer.h5`](data/2327/) — ~77 MB (lazy-mode, v2 bundle).

📖 **Full guide**: [`book/`](book/) — bookdown with overview, install,
tutorial, YAML schema, bundle/shard deep-dive, R + Python API
reference, Shiny app walkthrough, troubleshooting. Built output lives
in [`docs/`](docs/) and is GitHub-Pages-ready.

---

## Locked design decisions

| Question | Choice | Rationale |
| --- | --- | --- |
| Backend           | Standalone, file-based              | No Java, no graph DB, no MySQL; simplest for sharing studies offline. |
| Storage           | Lazy HDF5 bundle + on-disk shard tree | Bundle holds only metadata + per-matrix gene indexes (~80 MB even for studies whose underlying matrices are GB-scale). Matrix values stay on disk as gzipped per-gene shards and are read on demand. |
| Layout            | One subfolder per study under a shared root | `<root>/<studyID>/<studyID>.scminer.h5` + co-located graph layout + shard tree. Multiple studies live side-by-side; `run_browser(root)` discovers all of them. |
| Python framework  | Shiny for Python                    | Same reactive model as R Shiny — minimises divergence between the two packages. |
| v1 feature scope  | Six core plot tabs + cluster table + downsampling + multi-study browser | Cluster / Heatmap / Bubble / Feature / Violin / Network. Cluster table with checkbox selection and live cell counts. Sampling % auto-applies above 65K cells. Default-genes auto-load. 3-level nested tabs (Gene → CellType → Expression/TF/SIG) for the per-gene panels. |

---

## The bundle format (R ↔ Python contract)

A small `.scminer.h5` file (`bundleVersion = 2`). Strings are
variable-length UTF-8. The bundle stores **metadata + per-matrix gene
indexes only** — expression / activity values stay on disk as gzipped
per-gene shards and are read lazily by `gene_values()`.

```
/meta/
  studyID, studyAbbr, longTitle, shortTitle, species, coordinate  (string)
  bundleVersion                                                   (int, = 2)

/cells/
  cellID    (N,)   string
  cellType  (N,)   string
  cellGroup (N,)   string
  coord1    (N,)   float64
  coord2    (N,)   float64

/clusters/
  cellType (K,)   string
  count    (K,)   int32
  color    (K,)   string  (hex)
  label_1  (K,)   float64    optional — x for cluster label placement
  label_2  (K,)   float64    optional — y for cluster label placement

/genes/
  symbol   (G,)   string         master picker list

/index/                          optional, "which genes are available"
  expression   (G_e,) string     present iff expression shards exist
  activity_tf  (G_t,) string
  activity_sig (G_s,) string

/defaults/                       optional, "load these on app startup"
  genes        (D,)   string

/network_tf/                     optional
  source, target, cellType                       string
  mi, pearson, spearman, rho, pvalue             float64
/network_sig/                    optional, same columns
```

The bundle is normally co-located with the shard tree
(`<shard_dir>/<studyID>.scminer.h5` and
`<shard_dir>/{expression_files,activity_files}/<studyID>/...`). The
reader auto-discovers `shard_dir = dirname(bundle_path)` unless
overridden.

**Compatibility rules** — both packages MUST:

- Treat unknown groups as future-bundle-versions and ignore them.
- Tolerate any optional group being absent.
- Read scalar string datasets that the R writer emits as 1-element
  arrays (the Python reader unwraps with `arr.item()` / `arr[0]`).
- Tolerate a missing parent group when checking optional paths
  (`/index/expression` may be absent because the whole `/index` group
  is absent).

---

## Package layouts

```
scminerViewer/                              scminer_viewer/
├── DESCRIPTION                             ├── pyproject.toml
├── LICENSE                                 ├── README.md
├── NAMESPACE  (roxygen-generated)          ├── src/scminer_viewer/
├── README.md                               │   ├── __init__.py
├── R/                                      │   ├── __main__.py
│   ├── utils.R                shared       │   ├── cli.py
│   ├── write_bundle.R         HDF5 writer  │   ├── data.py
│   ├── read_bundle.R          load_study   │   │     Study + load_study + gene_values
│   ├── graph_read.R           read_graph   │   ├── app.py
│   ├── write_graph.R          internal     │   │     build_app + run_app
│   ├── prepare_study.R        3 entries +  │   ├── browser.py
│   │                           staged      │   │     discover_studies + build/run_browser
│   │                           helpers     │   └── plots/
│   ├── plots.R                7 plot       │       ├── _common.py
│   │                           helpers     │       ├── cluster.py
│   ├── app_ui.R               Shiny UI     │       ├── feature.py
│   ├── app_server.R           Shiny server │       ├── violin.py
│   ├── run_app.R              run/build_app│       ├── heatmap.py
│   └── run_browser.R          multi-study  │       ├── bubble.py
├── inst/                                   │       └── network.py
│   ├── extdata/example_config.yml          └── tests/
│   └── scripts/                                ├── conftest.py     fixtures via Rscript
│       ├── build_2327_bundle.R                 ├── test_data.py
│       └── install.R          safe_install     ├── test_plots.py
├── tests/testthat/                             └── test_browser.py
│   ├── helper-fixtures.R
│   ├── test-bundle-roundtrip.R
│   ├── test-graph-read.R
│   ├── test-shard-reader.R
│   ├── test-prepare-study.R
│   ├── test-staged-helpers.R
│   ├── test-browser.R
│   └── test-app.R
└── man/  (roxygen-generated *.Rd files)

book/                                       docs/
├── _bookdown.yml, _output.yml              └── rendered gitbook html, ready for
├── index.Rmd + 9 chapter Rmds                 GitHub Pages from /docs on main
├── README-deploy.md
├── render.R
├── style.css, preamble.tex
```

---

## Public API

### `scminerViewer` (R)

#### Orchestrators (one-shot prep)

| Function | Purpose |
| --- | --- |
| `prepare_study(config_path, emit, verbose)`                    | YAML-driven entry. Requires `yaml` + `Biobase`. |
| `prepare_study_from_eset(out_dir, expression_eset, ...)`       | Accepts a Biobase `ExpressionSet`. Splits activity rows into TF/SIG by `_TF`/`_SIG` row suffix. Requires `Biobase`. |
| `prepare_study_data(out_dir, meta, cells, clusters, genes, expression, ..., default_genes, emit, verbose)` | Lowest-level orchestrator — plain R structures only. Writes to `<out_dir>/<studyID>/`. |

#### Staged helpers (granular control / debug)

Call these directly when you want to inspect an intermediate before
committing to disk; `prepare_study()` is composed of them.

| Function | Purpose |
| --- | --- |
| `load_study_config(config_path)`                  | Parse + validate YAML; fill defaults. Pure parsing — does not touch any RDS / TSV. |
| `extract_cells(eset, ...)`                        | `pData(eset)` → cells data.frame. |
| `extract_genes(eset, gene_symbol_col)`            | `fData(eset)` → character vector of gene symbols. |
| `extract_expression(eset, genes)`                 | `exprs(eset)` → sparse `Matrix` (genes × cells). |
| `extract_activity(activity_eset, master_genes)`   | Splits rows by `_TF`/`_SIG` suffix; returns `list(tf, sig)` reindexed to `master_genes`. |
| `read_networks(path)`                             | Parse a scMINER networks TSV; returns `list(tf, sig)` data.frames. |
| `fill_clusters(cells, clusters = NULL, palette = "npg")` | Auto-populate a clusters data.frame: counts, colours (ggsci palette — default Nature/NPG), and label centroids (`mean(coord1)` / `mean(coord2)` per cluster). Existing non-NA values preserved. |

#### Bundle + study handles

| Function | Purpose |
| --- | --- |
| `write_bundle(bundle_path, meta, cells, clusters, genes, expression_genes, activity_tf_genes, activity_sig_genes, default_genes, network_tf, network_sig, overwrite)` | Write the `.scminer.h5` bundle. Indexes + metadata only. |
| `load_study(bundle_path, shard_dir = NULL)`       | Read a `.scminer.h5` into an S3 `scminer_study` list. `shard_dir` defaults to `dirname(bundle_path)`. |
| `gene_values(study, gene, relationship)`          | Lazily read one gene's row from the shard tree. Cached per gene. |
| `read_graph_study(data_dir, study_id)`            | Reconstruct study inputs from the on-disk graph layout. Per-matrix indexes from the manifest CSVs; matrix values not loaded. Auto-detects both flat (`<data_dir>/Study/...`) and wrapped (`<data_dir>/<sid>/Study/...`) layouts. |

#### Shiny entry points

| Function | Purpose |
| --- | --- |
| `run_app(bundle_path, shard_dir = NULL, ...)`     | Launch the single-study Shiny app. |
| `build_app(bundle_path, shard_dir = NULL)`        | Build a single-study `shiny.appobj` without launching. |
| `discover_studies(root_dir)`                      | Return one row per `<studyID>/<studyID>.scminer.h5` bundle found. |
| `run_browser(root_dir, shard_dir = NULL, ...)`    | Launch the **multi-study** browser — card-grid index → drill-in. |
| `build_browser(root_dir, shard_dir = NULL)`       | Build the multi-study `shiny.appobj` without launching. |

#### Dev helpers

| Helper | Purpose |
| --- | --- |
| `inst/scripts/install.R` → `safe_install(pkg, lib)` | Unloads namespace, GCs, wipes the install dir, then runs `R CMD INSTALL` in a `callr` subprocess. Sidesteps the *"lazy-load database is corrupt / interrupted promise evaluation"* trap that fires when the package is reinstalled while another R session has it loaded. |
| `inst/scripts/build_2327_bundle.R`                 | One-shot script to migrate the on-disk 2327 study into the `<data_dir>/<sid>/<sid>.scminer.h5` form. |

### `scminer_viewer` (Python)

| Symbol | Purpose |
| --- | --- |
| `load_study(bundle_path, shard_dir=None) -> Study`        | Read a `.scminer.h5` into a `Study` dataclass. |
| `Study`                                                    | Dataclass with `meta`, `cells` (pandas, indexed by cellID), `clusters` (pandas, indexed by cellType), `genes` (numpy), `expression_index` / `activity_tf_index` / `activity_sig_index` (numpy or None), `default_genes` (numpy or None), `network_tf` / `network_sig` (pandas or None), `shard_dir`. |
| `Study.gene_values(gene, relationship)`                    | Lazily read one gene's row from the shard tree; numpy ndarray aligned to `Study.cells.index`, or `None` if absent. Cached per gene. |
| `Study.has_gene(gene, relationship)`                       | True if `gene` is in the corresponding bundle index. |
| `build_app(bundle_path, shard_dir=None) -> shiny.App`      | Build a single-study Shiny app. |
| `run_app(bundle_path, shard_dir=None, host, port, launch_browser)` | Build + serve the single-study app via uvicorn (blocking). |
| `discover_studies(root_dir) -> pandas.DataFrame`           | One row per `<studyID>/<studyID>.scminer.h5` bundle found. |
| `build_browser(root_dir, shard_dir=None) -> shiny.App`     | Build the multi-study browser app. |
| `run_browser(root_dir, shard_dir=None, host, port, launch_browser)` | Serve the multi-study browser. |
| `plots.{cluster,feature,violin,heatmap,bubble,network}_plot(study, ..., cell_mask)` | Plotly figure factories. Every helper takes an optional `cell_mask` (boolean ndarray) that intersects with the cluster filter — used by the Sampling % control. |

**CLI** (installed via `pip install -e .`):

```text
scminer-viewer info  <bundle> [--shard-dir DIR]
scminer-viewer list  <root_dir>
scminer-viewer run   <bundle>  [--shard-dir DIR] [--host H] [--port N] [--no-browser]
scminer-viewer browse <root_dir> [--shard-dir DIR] [--host H] [--port N] [--no-browser]
```

---

## Input layouts the R package reads

Two distinct trees, both rooted under `data_dir` (with `shard_dir`
defaulting to the same path).

### 1. Graph-import layout

```
<data_dir>/<studyID>/                                 (wrapped form)
├── Study/<studyID>_study.tsv                         (4 tab-sep fields)
├── Cell/<studyID>_n_cell.tsv                         (7 columns: cellID, cellID, cellType, cellGroup, coord1, coord2, coordinateName)
├── Gene/<studyID>_n_gene.tsv                         (one gene per line)
├── Network_TF_Activity/<studyID>_TF.tsv              (10 columns; reader keeps cols 1,2,5..10)
├── Network_SIG_Activity/<studyID>_SIG.tsv
├── study_meta/<studyID>_study_meta.csv               (cluster colors + label coords)
├── Study_Contains_Cell/, Study_Contains_Gene/        (relationship tsvs — not consumed)
├── study_gene_expression/<studyID>_expression.csv    (manifest)
├── study_gene_tf/<studyID>_activity_tf.csv
└── study_gene_sig/<studyID>_activity_sig.csv
```

The reader also accepts the **flat** form (`<data_dir>/Study/...`) and
the legacy shared `study_meta/study_meta.csv` form as fallbacks. The
`Cell` TSV duplicates `cellID` in columns 1 and 2 — an artifact of
the original `sapply` over a data.frame row. The reader tolerates
both 7-column and the older 6-column form; the writer preserves the
7-column form for byte-compatibility with the existing Java backend.

### 2. Sharded matrix tree

```
<shard_dir>/<studyID>/
├── expression_files/<studyID>/
│   ├── meta.csv                                      (one line, comma-separated cellIDs)
│   └── <letter>/<gene>.csv.gz                        (gzipped CSV, one row of N values)
└── activity_files/<studyID>/
    ├── meta.csv                                      (shared by both TF and SIG)
    ├── TF/<letter>/<gene>.csv.gz
    └── SIG/<letter>/<gene>.csv.gz
```

`<letter>` is `tolower(substr(gene, 1, 1))` if alphabetic, else `nm`.
Paths are derived from the gene name + matrix kind — the manifests'
`FileHeader` and `File` columns are **ignored**. The activity
`meta.csv` may list a subset of the cells (activity is computed only
for cells passing some upstream filter); the reader matches each
shard's columns to `cells$cellID` via a permutation index so missing
cells become zero columns automatically.

---

## How to run

```sh
# ----- R side ---------------------------------------------------------------
R CMD build scminerViewer
R CMD INSTALL scminerViewer_0.1.0.tar.gz

# Run the R test suite
Rscript -e 'testthat::test_dir("scminerViewer/tests/testthat")'

# Build / refresh the 2327 bundle from data/ (writes data/2327/2327.scminer.h5)
Rscript scminerViewer/inst/scripts/build_2327_bundle.R

# Launch the single-study Shiny app
Rscript -e 'scminerViewer::run_app("data/2327/2327.scminer.h5", port=8000)'

# OR launch the multi-study browser (card-grid index of every
# <studyID>/<studyID>.scminer.h5 found under the root)
Rscript -e 'scminerViewer::run_browser("data", port=8000)'

# Build the bookdown (output → docs/)
Rscript book/render.R

# ----- Python side ----------------------------------------------------------
python3 -m venv .venv && source .venv/bin/activate
pip install -e "scminer_viewer[dev]"

# Run the Python test suite (builds fixtures via Rscript)
pytest scminer_viewer/tests -q

# Inspect or launch
scminer-viewer info   data/2327/2327.scminer.h5
scminer-viewer list   data
scminer-viewer run    data/2327/2327.scminer.h5 --port 8000
scminer-viewer browse data                       --port 8000
```

The bundle is small (~80 MB for 2327) because the expression /
activity matrix values are never embedded; the apps fetch each gene's
`.csv.gz` from the shard tree on demand. If the
`expression_files/<studyID>/` or `activity_files/<studyID>/` trees are
missing, the relevant index is still in the bundle but `gene_values()`
returns `NULL` for any gene in that index and the apps show
"no data" placeholders.

---

## Documentation

| Doc | What it covers |
| --- | --- |
| [`README.md`](README.md)                              | Project front page — high-level overview, lazy-mode TL;DR, repo layout, three-command quick-start, status. |
| [`scminerViewer/README.md`](scminerViewer/README.md)  | R package — full API surface (orchestrators / staged helpers / bundle + Shiny), config YAML reference, multi-study tutorial, install + RStudio + safe-install, troubleshooting. |
| [`scminer_viewer/README.md`](scminer_viewer/README.md) | Python package — install, `Study` dataclass, CLI (`info` / `list` / `run` / `browse`), R-vs-Python differences. |
| [`book/`](book/)                                       | bookdown source: 9 chapters covering overview, install, tutorial, YAML schema, bundle/shards, R API, Python API, Shiny app walkthrough, troubleshooting. Render with `Rscript book/render.R` to `docs/`. |
| [`docs/`](docs/)                                       | Rendered gitbook, ready to serve from GitHub Pages (`Settings → Pages → main → /docs`). See [`book/README-deploy.md`](book/README-deploy.md) for deploy options. |

---

## Status

| Component | Status | Tests |
| --- | --- | --- |
| `write_bundle` / `load_study` (v2 bundle round-trip)     | done | `test-bundle-roundtrip.R` |
| Lazy `gene_values()` over the shard tree (R)             | done | `test-shard-reader.R` |
| `read_graph_study` (manifest-derived indexes; flat + wrapped layouts) | done | `test-graph-read.R` |
| `prepare_study*` (graph writer + bundle co-emit; `<out>/<sid>/` wrap) | done | `test-prepare-study.R` |
| Staged helpers (`load_study_config`, `extract_*`, `read_networks`) | done | `test-staged-helpers.R` |
| Cluster auto-fill (`fill_clusters` — ggsci colours + centroid labels) | done | `test-fill-clusters.R` |
| R Shiny app — 6 plot tabs + cluster table + downsampling + 3-level nested per-gene panels | done | `test-app.R` |
| Multi-study browser (`run_browser` / `discover_studies`) — R | done | `test-browser.R` |
| Python `load_study` + `Study.gene_values` (lazy)         | done | `test_data.py` |
| Python Shiny app (UI + 6 plot tabs, sampling, 3-level tabs, live cluster counts) | done | `test_plots.py` |
| Multi-study browser (`discover_studies` / `run_browser`) — Python | done | `test_browser.py` |
| `safe_install()` helper (sidesteps lazy-load corruption) | done | manual smoke test |
| Bookdown (9 chapters, gitbook + pdf output)              | done | `bookdown::render_book` round-trip |
| `2327.scminer.h5` deliverable                            | done | 77 MB at `data/2327/2327.scminer.h5`; 9861 expression / 925 TF / 4708 SIG indexed; values fetched lazily from shards |

**Test totals**: R **194** passing in `scminerViewer/tests/testthat/`;
Python **22** passing in `scminer_viewer/tests/`. Python fixtures are
built by calling `Rscript` → `prepare_study_data` (which writes both
the shard tree and the bundle), so the Python lazy reader is verified
against the same on-disk layout the R apps consume.

### Remaining / nice-to-have

- `R CMD check --as-cran` clean run (mostly lint-only diagnostics
  currently; no functional issues).
- sdist + wheel publishing for the Python package (`uv build` /
  `twine check`).
- Deeper cross-language parity tests: identical top-3 expressed genes
  per cluster, identical network neighbour counts.
- GitHub Pages auto-deploy: today the rendered `docs/` is committed to
  `main` and served from there; users who'd rather not commit
  artefacts can switch to a manual push to `gh-pages` (see
  `book/README-deploy.md`). A workflow-based auto-deploy was
  considered and explicitly removed.
