# scminerViewer

R package for preparing scMINER study data and viewing it via Shiny.
Writes the graph-import layout used by the scMINER portal **and** a
small `.scminer.h5` bundle file that the R Shiny app and the companion
[`scminer_viewer`](../scminer_viewer/) Python webui both consume.

**Lazy by design**: the bundle holds only study metadata + per-matrix
gene indexes; the actual expression and activity values stay on disk as
gzipped per-gene shards and are read on demand by `gene_values()` when
a user picks a gene. So a typical study bundle is tens of MB even when
the underlying matrices are hundreds of MB.

- **R 4.x** required. Hard deps: `hdf5r`, `Matrix`, `methods`,
  `data.table`, `R.utils`.
- Optional: `Biobase` (for ExpressionSet input), `yaml` (for the YAML
  config entry point), `shiny`, `bslib`, `plotly`, `DT`, `visNetwork`,
  `htmlwidgets` (for the Shiny app).

## Install

### Command line

```sh
R CMD build scminerViewer
R CMD INSTALL scminerViewer_0.1.0.tar.gz
```

### From R or RStudio

```r
# devtools (or pak / remotes work too) — install from the local source dir
install.packages("devtools")           # if not already installed
devtools::install("scminerViewer")     # run from the project root
```

Or, in **RStudio**:

1. **File → Open Project…** and pick `scminerViewer/` (the package
   folder is its own RStudio project).
2. **Build → Install Package** (or `Cmd/Ctrl + Shift + B`). RStudio
   runs `R CMD INSTALL` on the package and reloads it.
3. To install with all suggested deps (Shiny app, ExpressionSet input,
   YAML configs), tick **Build → Configure Build Tools… → Install and
   Restart** options, or run once:
   ```r
   install.packages(c("yaml", "shiny", "bslib", "plotly", "DT",
                       "visNetwork", "htmlwidgets"))
   # Bioconductor:
   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
   BiocManager::install("Biobase")
   ```

Verify the install:

```r
library(scminerViewer)
?prepare_study     # or any other exported function
```

## How it works (one minute)

A study on disk consists of **two co-located parts**:

1. The **graph layout** — `Cell/`, `Gene/`, `Network_*/`, `study_meta/`,
   manifest CSVs, plus the sharded per-gene `.csv.gz` tree under
   `expression_files/<studyID>/` and `activity_files/<studyID>/`.
2. The **bundle file** — a single `.scminer.h5` written next to the
   shard tree, e.g. `data/2327.scminer.h5`. It carries study metadata,
   cell + cluster info, the master gene list, gene indexes (which genes
   have shards in each matrix), optional default genes, and networks.

`load_study(bundle_path)` reads the bundle; `gene_values(study, gene,
relationship)` then reads exactly one `<gene>.csv.gz` per call from
the shard tree (the parent dir of the bundle by default, overridable
with `shard_dir = ...`).

## Quick start

```r
library(scminerViewer)

# Prepare from a scMINER ExpressionSet via a YAML config (see below for
# the schema). Writes the graph layout + a .scminer.h5 bundle.
prepare_study("path/to/config.yml")

# Or migrate an existing on-disk graph-layout study to a bundle without
# touching the shard tree (only manifests + metadata are read).
study <- read_graph_study(data_dir = "path/to/data", study_id = "2327")
write_bundle(
  bundle_path        = "path/to/data/2327.scminer.h5",
  meta               = study$meta,
  cells              = study$cells,
  clusters           = study$clusters,
  genes              = study$genes,
  expression_genes   = study$expression_genes,
  activity_tf_genes  = study$activity_tf_genes,
  activity_sig_genes = study$activity_sig_genes,
  network_tf         = study$network_tf,
  network_sig        = study$network_sig
)

# Inspect a bundle + lazily read one gene's row aligned to the cells
study <- load_study("path/to/data/2327.scminer.h5")
print(study)
vals <- gene_values(study, "Mrpl15", "Express_normalized")

# Launch the Shiny app on the bundle
run_app("path/to/data/2327.scminer.h5", port = 8000)
```

For programmatic pipelines that already have matrices in memory, use
`prepare_study_from_eset()` (Biobase ExpressionSet input) or
`prepare_study_data()` (plain R structures).

## Exported API

| Function | Purpose |
| --- | --- |
| `prepare_study(config_path, emit, verbose)`                    | YAML-driven entry. Requires `yaml` + `Biobase`. |
| `prepare_study_from_eset(out_dir, expression_eset, ...)`       | Accepts a Biobase ExpressionSet. Splits activity rows into TF/SIG by `_TF`/`_SIG` row suffix. Requires `Biobase`. |
| `prepare_study_data(out_dir, meta, cells, clusters, genes, expression, ..., default_genes, emit, verbose)` | Lowest-level entry — plain R structures only. |
| `write_bundle(bundle_path, meta, cells, clusters, genes, expression_genes, activity_tf_genes, activity_sig_genes, default_genes, network_tf, network_sig, overwrite)` | Write the `.scminer.h5` bundle. Indexes + metadata only; values stay on disk. |
| `read_graph_study(data_dir, study_id)`                         | Reconstruct study inputs from the on-disk graph layout. Per-matrix gene **indexes** come from the manifest CSVs; shard values are never loaded eagerly. |
| `load_study(bundle_path, shard_dir = NULL)`                    | Read a `.scminer.h5` into an S3 `scminer_study` list. `shard_dir` defaults to `dirname(bundle_path)`. |
| `gene_values(study, gene, relationship)`                       | Lazily read one gene's row from the shard tree (`Express_normalized`, `Activity_tf`, or `Activity_sig`). Aligned to `study$cells$cellID`; cached per gene. |
| `run_app(bundle_path, host, port, launch_browser, ...)`        | Launch the Shiny app for a bundle. |
| `build_app(bundle_path)`                                       | Build a `shiny.appobj` without launching it (for tests / embedding). |

## Config YAML

`prepare_study(config_path)` reads a YAML file pointing at a scMINER
ExpressionSet RDS (and optional activity + networks files) and writes
both the graph layout and the bundle. A fully annotated example ships
with the package; copy it as a starting point:

```r
file.copy(
  system.file("extdata", "example_config.yml", package = "scminerViewer"),
  "config.yml"
)
```

Then run:

```sh
Rscript -e 'scminerViewer::prepare_study("config.yml")'
```

Pass `emit = "bundle"` (or `"graph"`) to skip emitting the other format.

### Schema

| YAML key | Required | Default | Notes |
| --- | --- | --- | --- |
| `output`            | yes  | —             | Directory for graph layout + bundle. The bundle lands at `<output>/<studyID>.scminer.h5`. |
| `study.ID`          | yes  | —             | Numeric or string. |
| `study.studyAbbr`   | yes  | —             | Short identifier (used in URLs / labels). |
| `study.longTitle`   | yes  | —             | Full study description. |
| `study.shortTitle`  | yes  | —             | One-line label. |
| `species`           | no   | `""`          | e.g. `"Mus musculus"`. |
| `coordinate`        | no   | `"UMAP"`      | Both the meta field and the pData column prefix (`<coordinate>_1`, `<coordinate>_2`). |
| `cellID`            | no   | `"cellID"`    | pData column name carrying the cellID. Rownames are also assigned to this column if missing. |
| `cellType`          | no   | `"cellGroup"` | pData column name carrying the cell-type label. |
| `cellGroup`         | no   | `cellType`    | pData column name carrying the cell-group label. |
| `geneSymbol`        | no   | `"geneSymbol"`| fData column name carrying the gene symbol. |
| `input.expression`  | yes  | —             | Path to the expression `ExpressionSet` RDS. |
| `input.activity`    | no   | none          | Path to the activity `ExpressionSet` RDS. Rows ending in `_TF` / `.TF` go to `activity_tf`; `_SIG` / `.SIG` go to `activity_sig`. |
| `input.networks`    | no   | none          | Path to a TSV with columns `source, target, NetworkType, CellGroup, mi, pearson, spearman, rho, pvalue`. |

## Shard tree layout

The bundle's `expression_index` / `activity_tf_index` / `activity_sig_index`
list which genes have on-disk shards. `gene_values()` reads one shard
per call from this tree:

```
<shard_dir>/
├── expression_files/<studyID>/
│   ├── meta.csv                  # one line, comma-separated cellIDs
│   └── <letter>/<gene>.csv.gz    # gzipped 1-row CSV of N values
└── activity_files/<studyID>/
    ├── meta.csv                  # shared by TF and SIG (may be a cell subset)
    ├── TF/<letter>/<gene>.csv.gz
    └── SIG/<letter>/<gene>.csv.gz
```

`<letter>` is `tolower(substr(gene, 1, 1))` if alphabetic, else `nm`.
Paths are derived from the gene name + matrix kind — the manifests'
`File` / `FileHeader` columns are ignored. Cell-ID column ordering is
read once from the per-matrix `meta.csv`; a permutation index aligns
shard columns to `cells$cellID`, so missing cells (activity has a
subset of the cells expression has) become zero columns automatically.

The full bundle schema (HDF5 group layout, dtypes, optional vs required
groups) lives in [`../IMPLEMENTATION.md`](../IMPLEMENTATION.md).

## Tests

```sh
Rscript -e 'testthat::test_dir("scminerViewer/tests/testthat")'
```

122 tests across `test-bundle-roundtrip.R`, `test-graph-read.R`,
`test-shard-reader.R`, `test-prepare-study.R`, and `test-app.R` (plus
shared `helper-fixtures.R`). Set `SCMINER_DATA_DIR=/path/to/data` to
enable the integration tests that read real graph-layout outputs.

## Reproducing the 2327 bundle

```sh
Rscript scminerViewer/inst/scripts/build_2327_bundle.R
# or with explicit args:
Rscript scminerViewer/inst/scripts/build_2327_bundle.R [data_dir] [out_path] [study_id]
```

Defaults: `data_dir = data`, `study_id = 2327`,
`out_path = <data_dir>/<study_id>.scminer.h5` — by default the bundle
lands next to the shard tree so `load_study()` auto-discovers shards
without an explicit `shard_dir`. The script tolerates missing matrices
and warns on the path of any file it can't find.
