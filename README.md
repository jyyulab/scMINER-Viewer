# scMINER Viewer

Standalone, offline viewer for [scMINER](https://github.com/jyyulab/scMINER)
single-cell studies. Two sibling packages share a small HDF5 bundle
format and the on-disk shard tree produced by the existing scMINER
data-prep pipeline:

| Package | Lang. | What it does |
| --- | --- | --- |
| [`scminerViewer/`](scminerViewer/)   | R      | Prepares scMINER studies (writes the graph-import layout + the `.scminer.h5` bundle) and serves a Shiny app that mirrors the layout at <https://scminer.stjude.org/study/Tregs>. |
| [`scminer_viewer/`](scminer_viewer/) | Python | Reads the same `.scminer.h5` bundle and serves a Shiny-for-Python webui with the same six-tab layout. |

Both apps are file-based: no Java services, no graph DB, no SQL. Drop
a study directory on disk and either package can render it.

## Lazy by design

A study consists of two co-located parts:

1. **Bundle** — `data/<studyID>.scminer.h5`. ~80 MB for the 8464-cell
   2327 study. Holds study metadata, cell + cluster info, the master
   gene list, per-matrix gene **indexes** (which genes exist in
   expression / activity_tf / activity_sig), optional default genes,
   and the networks. **No matrix values.**
2. **Shard tree** — gzipped per-gene CSVs under
   `data/expression_files/<studyID>/<letter>/<gene>.csv.gz` and
   `data/activity_files/<studyID>/{TF,SIG}/<letter>/<gene>.csv.gz`.
   This is the bulk of the data.

At startup the apps load only the bundle. When the user selects a gene
(or for default genes), the relevant `<gene>.csv.gz` is read once and
cached. Heatmap / Bubble / Feature / Violin tabs all use this lazy
accessor; the Cluster Plot and the Clusters table need no shard reads
at all.

## Repository layout

```
scMINER-Viewer/
├── IMPLEMENTATION.md              # full architecture & status doc
├── README.md                      # you are here
├── data/                          # the source study tree
│   ├── 2327.scminer.h5            # 77 MB bundle (lazy-mode, v2)
│   ├── Cell/, Gene/, Network_*/   # graph-import metadata TSVs
│   ├── study_meta/<studyID>_study_meta.csv
│   ├── study_gene_{expression,tf,sig}/<studyID>_*.csv  (manifests)
│   ├── expression_files/<studyID>/{meta.csv, <letter>/<gene>.csv.gz}
│   └── activity_files/<studyID>/{meta.csv, TF/<letter>/..., SIG/<letter>/...}
├── scminerViewer/                 # R package
│   ├── R/, tests/, inst/extdata/example_config.yml
│   └── inst/scripts/build_2327_bundle.R
└── scminer_viewer/                # Python package
    └── src/scminer_viewer/{data,app,cli}.py + plots/
```

## Quick start

### R side

**Install** — either from the shell:

```sh
R CMD build scminerViewer
R CMD INSTALL scminerViewer_0.1.0.tar.gz
```

…or from R / RStudio:

```r
devtools::install("scminerViewer")     # run from the project root
```

In RStudio, you can also **File → Open Project…** the `scminerViewer/`
folder and use **Build → Install Package** (`Cmd/Ctrl + Shift + B`).

**Use it**:

```sh
# Build the 2327 bundle from the data/ tree
Rscript scminerViewer/inst/scripts/build_2327_bundle.R

# Launch the Shiny app
Rscript -e 'scminerViewer::run_app("data/2327.scminer.h5", port = 8000)'
```

See [`scminerViewer/README.md`](scminerViewer/README.md) for the YAML
config schema, the full exported API, and how to prepare a fresh study
from a Biobase `ExpressionSet`.

### Python side

```sh
python3 -m venv .venv && source .venv/bin/activate
pip install -e "scminer_viewer[dev]"

scminer-viewer info data/2327.scminer.h5
scminer-viewer run  data/2327.scminer.h5 --port 8000
```

See [`scminer_viewer/README.md`](scminer_viewer/README.md) for the
`Study` dataclass + programmatic API.

## Bundle format (R ↔ Python contract)

Documented in full in [`IMPLEMENTATION.md`](IMPLEMENTATION.md). Highlights:

- Single HDF5 file, `bundleVersion = 2`.
- Strings are UTF-8.
- Three optional gene-index datasets under `/index/` enumerate which
  genes have shards in each matrix; the apps gate gene selection on
  these.
- An optional `/defaults/genes` dataset lists genes the app should
  auto-load on startup (mirrors the Vue portal's `preGenes`).
- The bundle and the shard tree are normally co-located:
  `load_study(bundle_path)` defaults `shard_dir = dirname(bundle_path)`.

## Status

- **R**: 122/122 tests pass (`scminerViewer/tests/testthat/`).
- **Python**: 17/17 tests pass (`scminer_viewer/tests/`). Fixtures are
  built by shelling out to `Rscript`, so the Python reader is verified
  against bytes produced by the R writer.
- **2327 bundle**: 77 MB at `data/2327.scminer.h5`. 8464 cells × 9861
  genes × 3 clusters; 9861 expression / 925 TF / 4708 SIG genes
  indexed; 743K total network edges; values fetched lazily from the
  shard tree.

## Out of scope for v1

- Private studies, authentication, user accounts.
- Downsampling slider, drag-to-reorder gene tags, per-cluster color
  picker, dot-size live-redraw on feature plots.
- Multi-nest-study tabs (one bundle = one study).
- Submitting your own data through the app
  (`submitfile.vue`, `visualizeyourdata.vue` in the original portal).
- Writing data back to a graph DB or any other DB.
