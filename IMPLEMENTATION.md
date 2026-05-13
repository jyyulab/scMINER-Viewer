# scMINER Viewer — Architecture & Status

Two sibling packages for viewing scMINER studies offline, both file-based
and standalone (no Java services, no graph DB, no MySQL):

| Package | Language | Path | Role |
| --- | --- | --- | --- |
| `scminerViewer`  | R      | [`scminerViewer/`](scminerViewer/)   | Refactor of `scMINER-portal-datapre-R/`. Emits the existing graph-import layout *and* an HDF5 study bundle. Ships a Shiny app reading the bundle. |
| `scminer_viewer` | Python | [`scminer_viewer/`](scminer_viewer/) | Loads the same HDF5 bundle and serves a Shiny-for-Python webui mirroring the Vue layout at <https://scminer.stjude.org/study/Tregs>. |

Deliverable: [`data-bundles/2327.scminer.h5`](data-bundles/) — 76 MB.

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

A single `.scminer.h5` file. Strings are variable-length UTF-8; sparse
matrices use the scipy CSR convention (zero-based indices) so
`scipy.sparse.csr_matrix((data, indices, indptr), shape=shape)`
reconstructs them directly. Optional groups are absent (not empty) when
the underlying data isn't supplied.

```
/meta/
  studyID, studyAbbr, longTitle, shortTitle, species, coordinate  (string)
  bundleVersion                                                   (int, =1)

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
  symbol   (G,)   string

/expression/                  optional, sparse CSR (genes × cells)
  data    (nnz,) float64
  indices (nnz,) int32       zero-based column indices
  indptr  (G+1,) int32
  shape   (2,)   int32       [G, N]

/activity_tf/                 optional, same shape as /expression
/activity_sig/                optional, same shape as /expression

/network_tf/                  optional
  source, target, cellType                       string
  mi, pearson, spearman, rho, pvalue             float64
/network_sig/                 optional, same columns
```

**Compatibility rules** — both packages MUST:

- Treat unknown groups as future-bundle-versions and ignore them.
- Tolerate any optional group being absent.
- Read scalar string datasets that the R writer emits as 1-element
  arrays (the Python reader unwraps with `arr.item()` / `arr[0]`).

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
| `write_bundle(bundle_path, meta, cells, clusters, genes, ..., overwrite)` | Write a `.scminer.h5` bundle directly. |
| `load_study(bundle_path)`                                     | Read a `.scminer.h5` into an S3 `scminer_study` list. |
| `read_graph_study(data_dir, study_id, shard_dir, load_expression, load_activity_tf, load_activity_sig, verbose)` | Reconstruct bundle inputs from the on-disk graph layout, optionally including per-gene shards. |
| `run_app(bundle_path, host, port, launch_browser, ...)`       | Launch the Shiny app for a bundle. |
| `build_app(bundle_path)`                                      | Build a `shiny.appobj` without launching (for tests / embedding). |

### `scminer_viewer` (Python)

| Symbol | Purpose |
| --- | --- |
| `load_study(bundle_path) -> Study`              | Read a `.scminer.h5` into a `Study` dataclass. |
| `Study`                                          | `meta`, `cells` (pandas, indexed by cellID), `clusters` (pandas, indexed by cellType), `genes` (numpy), `expression`/`activity_tf`/`activity_sig` (scipy CSR or `None`), `network_tf`/`network_sig` (pandas or `None`). |
| `Study.gene_values(gene, relationship)`         | Dense ndarray row for one gene, or `None` if missing. |
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
├── Study/<studyID>_study.tsv               (4 tab-sep fields: id, abbr, longTitle, shortTitle)
├── Cell/<studyID>_n_cell.tsv               (7 columns: cellID, cellID, cellType, cellGroup, coord1, coord2, coordinateName)
├── Gene/<studyID>_n_gene.tsv               (one gene per line)
├── Network_TF_Activity/<studyID>_TF.tsv    (10 columns; reader keeps cols 1,2,5..10)
├── Network_SIG_Activity/<studyID>_SIG.tsv  (same shape)
├── study_meta/study_meta.csv               (StudyID, StudyAbbr, CellType, CellGroup, Color, Label_1, Label_2, NetworkCellType, [Species])
├── Study_Contains_Cell/, Study_Contains_Gene/   (relationship tsvs — not consumed by the reader)
├── study_gene_expression/<studyID>_expression.csv   (manifest: GeneSymbol, Species, StudyID, StudyAbbr, Type, FileHeader, File)
├── study_gene_tf/<studyID>_activity_tf.csv
└── study_gene_sig/<studyID>_activity_sig.csv
```

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
│   ├── meta.csv                            (one line, comma-separated cellIDs in shard column order)
│   └── <letter>/<gene>.csv.gz              (gzipped CSV, one row of N comma-separated values)
└── activity_files/
    ├── TF/
    │   ├── meta.csv
    │   └── <letter>/<gene>.csv.gz
    └── SIG/
        ├── meta.csv
        └── <letter>/<gene>.csv.gz
```

`<letter>` is `tolower(substr(gene, 1, 1))` if alphabetic, else `nm`.
The manifests' `FileHeader` and `File` columns are **ignored** by the
reader — paths are derived from the gene name + matrix kind. The
manifest is used only to enumerate which genes have shards. The
`meta.csv` for each matrix type carries the cell-ID column order for
every shard beneath it; the reader matches each shard's columns to
`cells$cellID` via a permutation index, so column reordering / partial
cell coverage is handled automatically (missing cells become zeros with
a warning).

---

## How to run

```sh
# ----- R side ---------------------------------------------------------------
R CMD build scminerViewer
R CMD INSTALL scminerViewer_0.1.0.tar.gz

# Run the R test suite
Rscript -e 'testthat::test_dir("scminerViewer/tests/testthat")'

# Build / refresh the 2327 bundle from data/
Rscript scminerViewer/inst/scripts/build_2327_bundle.R \
  data data-bundles/2327.scminer.h5 2327

# Launch the R Shiny app
Rscript -e 'scminerViewer::run_app("data-bundles/2327.scminer.h5", port=8000)'

# ----- Python side ----------------------------------------------------------
python3 -m venv .venv && source .venv/bin/activate
pip install -e "scminer_viewer[dev]"

# Run the Python test suite (builds fixtures via Rscript)
pytest scminer_viewer/tests -q

# Inspect or launch
scminer-viewer info data-bundles/2327.scminer.h5
scminer-viewer run  data-bundles/2327.scminer.h5 --port 8000
```

The bundle skips matrices it can't find on disk and emits a warning
naming the missing file; once the shard tree
(`<shard_dir>/expression_files/`,
`<shard_dir>/activity_files/{TF,SIG}/`) is populated, re-running the
build script will populate `expression`, `activity_tf`, `activity_sig`
with no code changes. Both apps will then automatically light up the
Heatmap / Bubble / Feature / Violin tabs.

---

## Status

| Component | Status | Tests |
| --- | --- | --- |
| `write_bundle` / `load_study` (HDF5 round-trip)          | done | `test-bundle-roundtrip.R` |
| `read_graph_study` (data/ TSVs + sharded matrices)       | done | `test-graph-read.R`, `test-shard-reader.R` |
| `prepare_study*` (graph writer refactor + bundle co-emit)| done | `test-prepare-study.R` |
| R Shiny app (UI + 6 plot tabs)                           | done | `test-app.R` (construction + helpers) |
| Python `load_study` + `Study` dataclass                  | done | `test_data.py` |
| Python Shiny app (UI + 6 plot tabs)                      | done | `test_plots.py` (helpers) |
| `2327.scminer.h5` deliverable                            | done | 8464 cells × 9861 genes; expression/activity slots awaiting shard upload |

**Test totals**: R 119 passing in `scminerViewer/tests/testthat/`;
Python 14 passing in `scminer_viewer/tests/`. Python fixtures are built
by calling `Rscript` → `write_bundle`, so the Python reader is verified
against bytes produced by the R writer (the strongest source of
cross-language parity).

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
