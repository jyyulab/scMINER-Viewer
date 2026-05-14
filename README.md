# scMINER Viewer

Standalone, offline viewer for [scMINER](https://github.com/jyyulab/scMINER)
single-cell studies. Two sibling packages share a small HDF5 bundle
format and the on-disk shard tree produced by the existing scMINER
data-prep pipeline:

| Package | Lang. | What it does |
| --- | --- | --- |
| [`scminerViewer/`](scminerViewer/)   | R      | Prepares scMINER studies (writes the graph-import layout + the `.scminer.h5` bundle) and serves a Shiny app + multi-study browser. |
| [`scminer_viewer/`](scminer_viewer/) | Python | Reads the same `.scminer.h5` bundle; serves a Shiny-for-Python single-study viewer and the same multi-study browser. |

📖 **Full guide**: [`book/`](book/) — a complete bookdown with
overview, installation, tutorial, YAML reference, bundle/shard
deep-dive, R + Python API reference, Shiny app walkthrough, and
troubleshooting.

| To view | How |
| --- | --- |
| Online (GitHub Pages)  | Browse `https://<user>.github.io/<repo>/` once you've enabled Pages — see [`book/README-deploy.md`](book/README-deploy.md) for the one-time setup. |
| Locally               | `Rscript book/render.R` then open `docs/index.html`. |
| Source                | Read the `book/*.Rmd` chapters directly on GitHub (they render as markdown — readable but unlinked). |

Both apps are file-based: no Java services, no graph DB, no SQL. Drop
a study directory on disk and either package can render it.

## Lazy by design + multi-study by default

Each study lives in **its own subfolder** under a shared root so one
root can host many studies side-by-side:

```
data/                              # multi-study root → pass to run_browser()
├── 2327/
│   ├── 2327.scminer.h5            # ~80 MB bundle (metadata + gene indexes)
│   ├── Cell/, Gene/, Network_*/, study_meta/, study_gene_*/
│   ├── expression_files/2327/<letter>/<gene>.csv.gz   # sharded values
│   └── activity_files/2327/{meta.csv, TF/, SIG/}
├── 9999/
│   └── …
```

Inside each study folder:

1. **Bundle** — `<studyID>/<studyID>.scminer.h5`. Holds study metadata,
   cell + cluster info, master gene list, per-matrix gene **indexes**
   (which genes exist in expression / activity_tf / activity_sig),
   optional default genes, networks. **No matrix values.**
2. **Shard tree** — gzipped per-gene CSVs. The bulk of the data.

At startup the apps load only the bundle. When the user selects a gene
(or for default genes), the relevant `<gene>.csv.gz` is read once and
cached. Heatmap / Bubble / Feature / Violin tabs all use this lazy
accessor; the Cluster Plot and the Clusters table need no shard reads
at all. `run_browser(root_dir)` discovers every study under the root
and presents them as a card grid; clicking a card drills into that
study's viewer (back link to return).

## Repository layout

```
scMINER-Viewer/
├── IMPLEMENTATION.md                  # architecture & status doc
├── README.md                          # you are here
├── data/                              # multi-study root
│   └── 2327/                          # per-study subfolder
│       ├── 2327.scminer.h5
│       ├── Cell/, Gene/, Network_*/
│       ├── study_meta/2327_study_meta.csv
│       ├── study_gene_{expression,tf,sig}/2327_*.csv
│       ├── expression_files/2327/{meta.csv, <letter>/<gene>.csv.gz}
│       └── activity_files/2327/{meta.csv, TF/<letter>/..., SIG/<letter>/...}
├── scminerViewer/                     # R package
│   ├── R/, tests/, inst/extdata/example_config.yml
│   └── inst/scripts/build_2327_bundle.R
└── scminer_viewer/                    # Python package
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

**Use it** — three things you can do:

```sh
# (A) Build the 2327 bundle from an on-disk graph layout (no RDS needed)
Rscript scminerViewer/inst/scripts/build_2327_bundle.R
# → writes data/2327/2327.scminer.h5

# (B) Full prepare_study pipeline driven by a YAML config (needs the
#     scMINER ExpressionSet RDS / networks TSV; see data/example/2327.yml)
Rscript -e 'scminerViewer::prepare_study("data/example/2327.yml")'

# (C1) Multi-study browser — one URL per study, card-grid index
Rscript -e 'scminerViewer::run_browser("data", port = 8000)'

# (C2) Same, but shards still live in the legacy data/example/ tree
Rscript -e 'scminerViewer::run_browser("data", shard_dir = "data/example", port = 8000)'

# (D) Single-study app
Rscript -e 'scminerViewer::run_app("data/2327/2327.scminer.h5", port = 8000)'
```

Full multi-study walkthrough (prep → browse → add another) is in
[`scminerViewer/README.md → Tutorial`](scminerViewer/README.md#tutorial-from-zero-to-a-multi-study-browser).

A concrete, runnable YAML for the 2327 study lives at
[`data/example/2327.yml`](data/example/2327.yml) — copy it as a starting
point for new studies and edit the `input.*` paths to point at your
ExpressionSet RDS files. The package also ships a fully-annotated
template at
`system.file("extdata", "example_config.yml", package = "scminerViewer")`.

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

- **R**: 169/169 tests pass (`scminerViewer/tests/testthat/`).
- **Python**: 17/17 tests pass (`scminer_viewer/tests/`). Fixtures are
  built by shelling out to `Rscript`, so the Python reader is verified
  against bytes produced by the R writer.
- **2327 bundle**: 77 MB at `data/2327/2327.scminer.h5`. 8464 cells ×
  9861 genes × 3 clusters; 9861 expression / 925 TF / 4708 SIG genes
  indexed; 743K total network edges; values fetched lazily from the
  shard tree.

## Out of scope for v1

- Private studies, authentication, user accounts.
- Drag-to-reorder gene tags, per-cluster color picker, dot-size
  live-redraw on feature plots.
- Submitting your own data through the app
  (`submitfile.vue`, `visualizeyourdata.vue` in the original portal).
- Writing data back to a graph DB or any other DB.
