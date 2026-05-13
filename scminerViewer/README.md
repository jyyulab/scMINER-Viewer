# scminerViewer

R package for preparing scMINER study data and viewing it via Shiny.
Emits both the graph-import layout used by the scMINER portal **and** a
self-contained HDF5 study bundle that the Shiny app and the companion
[`scminer_viewer`](../scminer_viewer/) Python package both consume.

- **R 4.x** required. Hard deps: `hdf5r`, `Matrix`, `methods`,
  `data.table`, `R.utils`.
- Optional: `Biobase` (for ExpressionSet input), `yaml` (for the YAML
  config entry point), `shiny`, `bslib`, `plotly`, `DT`, `visNetwork`,
  `htmlwidgets` (for the Shiny app).

## Install

```sh
R CMD build scminerViewer
R CMD INSTALL scminerViewer_0.1.0.tar.gz
```

## Quick start

```r
library(scminerViewer)

# 1) YAML-driven: drop-in replacement for scMINER-portal-datapre-R/main.R
prepare_study("path/to/config.yml")

# 2) From a Biobase ExpressionSet (no YAML)
prepare_study_from_eset(
  out_dir         = "out",
  expression_eset = readRDS("expression.rds"),
  activity_eset   = readRDS("activity.rds"),
  networks_path   = "networks.tsv",
  meta = list(studyID="2327", studyAbbr="tex",
              longTitle="...", shortTitle="Tex",
              species="Mus musculus", coordinate="UMAP"),
  cell_type_col   = "cellGroup",
  emit            = c("graph", "bundle")
)

# 3) From plain R structures (no Biobase, no YAML)
prepare_study_data(
  out_dir  = "out",
  meta     = list(studyID="2327", studyAbbr="tex",
                  longTitle="...", shortTitle="Tex",
                  species="Mus musculus", coordinate="UMAP"),
  cells    = cells_df,       # cellID, cellType, cellGroup, coord1, coord2
  clusters = clusters_df,    # cellType, count (optional), color, label_1, label_2
  genes    = gene_vec,       # character()
  expression   = expr_sparse,   # G x N Matrix
  activity_tf  = tf_sparse,
  activity_sig = sig_sparse,
  network_tf   = tf_net_df,
  network_sig  = sig_net_df,
  emit = c("graph", "bundle")
)

# 4) Migrate an existing on-disk graph-layout study to a bundle
study <- read_graph_study(
  data_dir          = "path/to/data",
  study_id          = "2327",
  shard_dir         = "path/to/data",   # defaults to data_dir
  load_expression   = TRUE,
  load_activity_tf  = TRUE,
  load_activity_sig = TRUE,
  verbose           = TRUE
)
write_bundle("2327.scminer.h5",
             meta = study$meta, cells = study$cells,
             clusters = study$clusters, genes = study$genes,
             expression  = study$expression,
             activity_tf = study$activity_tf,
             activity_sig = study$activity_sig,
             network_tf = study$network_tf,
             network_sig = study$network_sig)

# 5) Serve the Shiny app
run_app("2327.scminer.h5", port = 8000)
```

## Exported API

| Function | Purpose |
| --- | --- |
| `prepare_study(config_path, emit, verbose)`                    | YAML-driven entry, drop-in for the original `main.R`. Requires `yaml` + `Biobase`. |
| `prepare_study_from_eset(out_dir, expression_eset, ...)`       | Accepts a Biobase ExpressionSet. Splits activity rows into TF/SIG by `_TF`/`_SIG` row suffix. Requires `Biobase`. |
| `prepare_study_data(out_dir, meta, cells, clusters, genes, expression, ..., emit, verbose)` | Lowest-level entry, plain R structures only. |
| `write_bundle(bundle_path, meta, cells, clusters, genes, ..., overwrite)` | Write a `.scminer.h5` bundle directly. |
| `load_study(bundle_path)`                                      | Read a `.scminer.h5` into an S3 `scminer_study` list (with `print` method). |
| `read_graph_study(data_dir, study_id, shard_dir, load_expression, load_activity_tf, load_activity_sig, verbose)` | Reconstruct study inputs from the on-disk graph layout + optional shard tree. |
| `run_app(bundle_path, host, port, launch_browser, ...)`        | Launch the Shiny app for a bundle. |
| `build_app(bundle_path)`                                       | Build a `shiny.appobj` without launching it (for tests / embedding). |

## Shard tree layout

`read_graph_study(..., load_expression = TRUE)` reads matrices from a
sharded tree alongside the graph-layout `Cell/`, `Gene/` etc.:

```
<shard_dir>/
├── expression_files/
│   ├── meta.csv                  # one line, comma-separated cellIDs
│   └── <letter>/<gene>.csv.gz    # gzipped 1-row CSV of N values
└── activity_files/
    ├── TF/{meta.csv, <letter>/<gene>.csv.gz}
    └── SIG/{meta.csv, <letter>/<gene>.csv.gz}
```

`<letter>` is `tolower(substr(gene,1,1))` if alphabetic, else `nm`. The
reader derives paths from the gene name and ignores the manifests'
`File` / `FileHeader` columns. Cell-ID column ordering in each shard is
read once from the top-level `meta.csv`; a permutation index aligns
shard columns to `cells$cellID`.

## Bundle format

Documented in [`../IMPLEMENTATION.md`](../IMPLEMENTATION.md). One HDF5
file per study; sparse matrices follow the scipy CSR convention
(zero-based indices) so the Python reader reconstructs them with
`scipy.sparse.csr_matrix((data, indices, indptr), shape=shape)` directly.

## Tests

```sh
Rscript -e 'testthat::test_dir("scminerViewer/tests/testthat")'
```

119 tests across six files (`test-bundle-roundtrip.R`,
`test-graph-read.R`, `test-shard-reader.R`, `test-prepare-study.R`,
`test-app.R`, plus shared `helper-fixtures.R`). Set
`SCMINER_DATA_DIR=/path/to/data` to enable the integration tests that
read real graph-layout outputs.

## Reproducing the 2327 bundle

```sh
Rscript scminerViewer/inst/scripts/build_2327_bundle.R \
  data data-bundles/2327.scminer.h5 2327 [shard_dir]
```

Positional args: `data_dir`, `out_path`, `study_id`, `shard_dir`
(defaults to `data_dir`). The script tolerates missing matrices,
warns with the path of the file it can't find, and produces the
bundle with whatever it could load.
