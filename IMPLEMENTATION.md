# scMINER Viewer — Architecture & Status

Two sibling packages for viewing scMINER studies offline, both file-based
and standalone (no Java services, no graph DB, no MySQL):

| Package | Language | Path | Role |
| --- | --- | --- | --- |
| `scminerViewer`  | R      | [`scminerViewer/`](scminerViewer/)   | Refactor of `scMINER-portal-datapre-R/`. Emits the existing graph-import layout *and* an HDF5 study bundle. Ships a Shiny app reading the bundle. |
| `scminer_viewer` | Python | [`scminer_viewer/`](scminer_viewer/) | Loads the same HDF5 bundle and serves a Shiny-for-Python webui mirroring the Vue layout at <https://scminer.stjude.org/study/Tregs>. |

Deliverable: [`data/2327/2327.scminer.h5`](data/2327/) — ~77 MB (lazy-mode, v2 bundle).

---

## Locked design decisions

| Question | Choice | Rationale |
| --- | --- | --- |
| Backend           | Standalone, file-based            | No Java services, no graph DB, no MySQL; simplest for sharing studies offline. |
| Storage           | HDF5 bundle alongside graph layout | Single file the apps load in one shot; preserves compatibility with the existing import pipeline. |
| Python framework  | Shiny for Python                  | Same reactive model as R Shiny — minimises divergence between the two packages. |
| v1 feature scope  | Core plots + cluster table        | Six tabs (Cluster / Heatmap / Bubble / Feature / Violin / Network), cluster table with row-selection visibility, gene selectize input. No auth, no downsampling slider, no drag-reorder. |

---

## The bundle format (R ↔ Python contract)

A small `.scminer.h5` file (bundle version 2). Strings are
variable-length UTF-8. The bundle stores **metadata + per-matrix gene
indexes only** — expression / activity values stay on disk as gzipped
per-gene shards and are read lazily by `gene_values()` when an app
actually needs them.

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

The bundle is normally co-located with the shard tree:
`<shard_dir>/<studyID>.scminer.h5` and
`<shard_dir>/{expression_files,activity_files}/<studyID>/...`. The
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
├── NAMESPACE                               ├── src/scminer_viewer/
├── README.md                               │   ├── __init__.py
├── R/                                      │   ├── __main__.py
│   ├── utils.R           shared helpers    │   ├── cli.py
│   ├── write_bundle.R    HDF5 writer       │   ├── data.py            Study + load_study
│   ├── read_bundle.R     load_study        │   ├── app.py             build_app / run_app
│   ├── graph_read.R      read_graph_study  │   └── plots/
│   ├── write_graph.R     internal helpers  │       ├── _common.py
│   ├── prepare_study.R   3 public entries  │       ├── cluster.py
│   ├── plots.R           6 plot helpers    │       ├── feature.py
│   ├── app_ui.R          Shiny UI          │       ├── violin.py
│   ├── app_server.R      Shiny server      │       ├── heatmap.py
│   └── run_app.R         run_app/build_app │       ├── bubble.py
├── inst/scripts/                           │       └── network.py
│   └── build_2327_bundle.R                 └── tests/
└── tests/testthat/                             ├── conftest.py        builds fixtures via Rscript
    ├── helper-fixtures.R                       ├── test_data.py
    ├── test-bundle-roundtrip.R                 └── test_plots.py
    ├── test-graph-read.R
    ├── test-shard-reader.R
    ├── test-prepare-study.R
    └── test-app.R
```

---

## Public API

### `scminerViewer` (R)

| Function | Purpose |
| --- | --- |
| `prepare_study(config_path, emit, verbose)`                   | YAML-driven entry; drop-in for the original `main.R`. Requires `yaml` + `Biobase`. |
| `prepare_study_from_eset(out_dir, expression_eset, ...)`      | Accepts a Biobase `ExpressionSet`. Requires `Biobase`. Activity matrix is split into TF/SIG by `_TF`/`_SIG` row suffix. |
| `prepare_study_data(out_dir, meta, cells, clusters, genes, expression, activity_tf, activity_sig, network_tf, network_sig, emit, verbose)` | Lowest-level entry; plain R structures. Writes graph layout, HDF5 bundle, or both per `emit`. |
| `write_bundle(bundle_path, meta, cells, clusters, genes, expression_genes, activity_tf_genes, activity_sig_genes, default_genes, network_tf, network_sig, overwrite)` | Write the `.scminer.h5` bundle (indexes only — no matrix values). |
| `load_study(bundle_path, shard_dir = NULL)`                   | Read a `.scminer.h5` into an S3 `scminer_study` list. `shard_dir` defaults to `dirname(bundle_path)`. |
| `gene_values(study, gene, relationship)`                      | Lazily read one gene's row from the on-disk shard tree, aligned to `study$cells$cellID`. Cached per gene. |
| `read_graph_study(data_dir, study_id)`                        | Reconstruct bundle inputs (cells, clusters, gene indexes from manifests, networks) from the on-disk graph layout. Does **not** read shard values. |
| `run_app(bundle_path, shard_dir = NULL, host, port, launch_browser, ...)` | Launch the single-study Shiny app for a bundle. |
| `build_app(bundle_path, shard_dir = NULL)`                    | Build a single-study `shiny.appobj` without launching (for tests / embedding). |
| `discover_studies(root_dir)`                                  | Return a data.frame describing every `<studyID>/<studyID>.scminer.h5` bundle found under `root_dir`. |
| `run_browser(root_dir, shard_dir = NULL, host, port, launch_browser, ...)` | Launch the **multi-study** browser — card-grid landing page; click → `?study=<id>` opens the standard viewer with a "← Back" link. |
| `build_browser(root_dir, shard_dir = NULL)`                   | Build the multi-study `shiny.appobj` without launching. |

### `scminer_viewer` (Python)

| Symbol | Purpose |
| --- | --- |
| `load_study(bundle_path, shard_dir=None) -> Study` | Read a `.scminer.h5` into a `Study` dataclass. `shard_dir` defaults to `Path(bundle_path).parent`. |
| `Study`                                          | `meta`, `cells` (pandas, indexed by cellID), `clusters` (pandas, indexed by cellType), `genes` (numpy), `expression_index`/`activity_tf_index`/`activity_sig_index` (numpy or `None`), `default_genes` (numpy or `None`), `network_tf`/`network_sig` (pandas or `None`), `shard_dir`. |
| `Study.gene_values(gene, relationship)`         | Lazily read one gene's values from the shard tree; numpy ndarray aligned to `Study.cells.index`, or `None` if absent. Cached per gene. |
| `Study.has_gene(gene, relationship)`            | True if `gene` is in the corresponding bundle index. |
| `build_app(bundle_path) -> shiny.App`           | Build a Shiny app without launching it. |
| `run_app(bundle_path, host, port, launch_browser)` | Launch the Shiny app (blocking, via uvicorn). |
| `plots.{cluster,feature,violin,heatmap,bubble,network}_plot(study, ...)` | Plotly figure factories used by the app. |
| CLI: `scminer-viewer info <bundle>`             | Print a one-line summary of a bundle. |
| CLI: `scminer-viewer run <bundle> --port N`     | Build + serve the Shiny app. |

---

## Input layouts the R package reads

Two distinct trees, both rooted under `data_dir` (with `shard_dir`
defaulting to the same path):

### 1. graph-import layout — produced by `scMINER-portal-datapre-R/`

```
<data_dir>/
├── Study/<studyID>_study.tsv                       (4 tab-sep fields: id, abbr, longTitle, shortTitle)
├── Cell/<studyID>_n_cell.tsv                       (7 columns: cellID, cellID, cellType, cellGroup, coord1, coord2, coordinateName)
├── Gene/<studyID>_n_gene.tsv                       (one gene per line)
├── Network_TF_Activity/<studyID>_TF.tsv            (10 columns; reader keeps cols 1,2,5..10)
├── Network_SIG_Activity/<studyID>_SIG.tsv          (same shape)
├── study_meta/<studyID>_study_meta.csv             (StudyID, StudyAbbr, CellType, CellGroup, Color, Label_1, Label_2, NetworkCellType, [Species])
├── Study_Contains_Cell/, Study_Contains_Gene/      (relationship tsvs — not consumed by the reader)
├── study_gene_expression/<studyID>_expression.csv  (manifest: GeneSymbol, Species, StudyID, StudyAbbr, Type, FileHeader, File)
├── study_gene_tf/<studyID>_activity_tf.csv
└── study_gene_sig/<studyID>_activity_sig.csv
```

The reader also accepts the legacy shared `study_meta/study_meta.csv` form
as a fallback (filtered by `StudyID`).

The `Cell` TSV duplicates `cellID` in columns 1 and 2 — an artifact of
the original `sapply` over a data.frame row. The reader tolerates this
(both 7-column and the older 6-column form). The writer preserves the
7-column form to stay byte-compatible with the existing Java backend.

The `study_meta.csv` `Species` column is optional; when absent, the
reader falls back to the `Species` column of the first matrix manifest
it can read.

### 2. Sharded matrix tree — referenced by the manifests

```
<shard_dir>/
├── expression_files/
│   └── <studyID>/
│       ├── meta.csv                        (one line, comma-separated cellIDs in shard column order)
│       └── <letter>/<gene>.csv.gz          (gzipped CSV, one row of N comma-separated values)
└── activity_files/
    └── <studyID>/
        ├── meta.csv                        (shared by both TF and SIG)
        ├── TF/<letter>/<gene>.csv.gz
        └── SIG/<letter>/<gene>.csv.gz
```

`<letter>` is `tolower(substr(gene, 1, 1))` if alphabetic, else `nm`.
The manifests' `FileHeader` and `File` columns are **ignored** by the
reader — paths are derived from the gene name + matrix kind. The
manifest is used only to enumerate which genes have shards.

The expression `meta.csv` records the cell-ID column order for every
expression shard. The activity `meta.csv` is shared between the TF and
SIG kinds, and may list a subset of the cells (activity is computed
only for cells passing some upstream filter). The reader matches each
shard's columns to `cells$cellID` via a permutation index, so column
reordering / partial cell coverage is handled automatically (cells
missing from the shard header become zero columns).

---

## How to run

```sh
# ----- R side ---------------------------------------------------------------
R CMD build scminerViewer
R CMD INSTALL scminerViewer_0.1.0.tar.gz

# Run the R test suite
Rscript -e 'testthat::test_dir("scminerViewer/tests/testthat")'

# Build / refresh the 2327 bundle from data/ (writes to data/2327.scminer.h5
# next to the shard tree)
Rscript scminerViewer/inst/scripts/build_2327_bundle.R

# Launch the R Shiny app
Rscript -e 'scminerViewer::run_app("data/2327/2327.scminer.h5", port=8000)'

# OR launch the multi-study browser (card-grid index of every
# <studyID>/<studyID>.scminer.h5 found under the root)
Rscript -e 'scminerViewer::run_browser("data", port=8000)'

# ----- Python side ----------------------------------------------------------
python3 -m venv .venv && source .venv/bin/activate
pip install -e "scminer_viewer[dev]"

# Run the Python test suite (builds fixtures via Rscript)
pytest scminer_viewer/tests -q

# Inspect or launch
scminer-viewer info data/2327/2327.scminer.h5
scminer-viewer run  data/2327/2327.scminer.h5 --port 8000
```

The bundle is small (~80 MB for 2327 — mostly networks) because the
expression / activity matrix values are never embedded; the apps fetch
each gene's `.csv.gz` from the shard tree on demand. If the
`expression_files/<studyID>/` or `activity_files/<studyID>/` trees are
missing, the relevant index will still be embedded but `gene_values()`
returns `NULL` for any gene in that index (the app gracefully shows
"no data" placeholders).

---

## Status

| Component | Status | Tests |
| --- | --- | --- |
| `write_bundle` / `load_study` (v2 bundle round-trip)     | done | `test-bundle-roundtrip.R` |
| Lazy `gene_values()` over the shard tree (R)             | done | `test-shard-reader.R` |
| `read_graph_study` (manifest-derived indexes)            | done | `test-graph-read.R` |
| `prepare_study*` (graph writer + bundle co-emit)         | done | `test-prepare-study.R` |
| R Shiny app (UI + 6 plot tabs, lazy plots)               | done | `test-app.R` (construction + helpers) |
| Python `load_study` + `Study.gene_values` (lazy)         | done | `test_data.py` |
| Python Shiny app (UI + 6 plot tabs)                      | done | `test_plots.py` (helpers) |
| `2327.scminer.h5` deliverable                            | done | 77 MB at `data/2327/2327.scminer.h5`; 9861 expression / 925 TF / 4708 SIG indexed; matrix values fetched lazily from shards |
| Multi-study browser (`run_browser` / `discover_studies`) | done | `test-browser.R` |

**Test totals**: R 169 passing in `scminerViewer/tests/testthat/`;
Python 17 passing in `scminer_viewer/tests/`. Python fixtures are
built by calling `Rscript` → `prepare_study_data` (which writes both
the shard tree and the bundle), so the Python lazy reader is verified
against the same on-disk layout the R apps consume.

### Remaining / nice-to-have

- `R CMD check --as-cran` clean run (mostly lint-only diagnostics
  currently; no functional issues).
- roxygen2-generated man pages (DESCRIPTION + NAMESPACE today).
- sdist + wheel publishing for the Python package (`uv build` /
  `twine check`).
- Deeper cross-language parity tests: identical top-3 expressed genes
  per cluster, identical network neighbour counts.

---

## Out of scope for v1

- Private studies, authentication, user accounts.
- Downsampling slider, drag-to-reorder gene tags, per-cluster color
  picker, dot-size live-redraw on feature plots.
- Multi-nest-study tabs (one bundle = one study).
- Submitting your own data (`submitfile.vue`, `visualizeyourdata.vue`).
- Writing data back to a graph DB or any other DB.
