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

> If `?prepare_study` returns "No documentation found", the `man/`
> directory wasn't generated. From the package source, run once:
>
> ```r
> roxygen2::roxygenise("scminerViewer")   # regenerates NAMESPACE + man/*.Rd
> ```
>
> then reinstall.

## How it works (one minute)

`prepare_study()` writes everything for one study into
`<out_dir>/<studyID>/`, so a single root can hold many studies:

```
<root>/
├── <studyID-A>/
│   ├── <studyID-A>.scminer.h5    bundle (metadata + gene indexes)
│   ├── Cell/, Gene/, Network_*/, study_meta/, study_gene_*/
│   ├── expression_files/<studyID-A>/<letter>/<gene>.csv.gz
│   └── activity_files/<studyID-A>/{meta.csv, TF/, SIG/}
├── <studyID-B>/
│   └── …
```

- `load_study(bundle_path)` reads the bundle (no matrix values);
- `gene_values(study, gene, relationship)` reads exactly one
  `<gene>.csv.gz` per call from the shard tree (`dirname(bundle_path)`
  by default; override with `shard_dir = ...`).
- `run_app(bundle_path)` serves the **single-study** viewer.
- `run_browser(root_dir)` serves the **multi-study** viewer: a
  card-grid index of every study at `<root>/<sid>/<sid>.scminer.h5`;
  clicking a card opens `?study=<sid>` with the single-study layout
  and a "← Back to studies" link.

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
study <- load_study("path/to/data/2327/2327.scminer.h5")
print(study)
vals <- gene_values(study, "Mrpl15", "Express_normalized")

# Launch the single-study Shiny app
run_app("path/to/data/2327/2327.scminer.h5", port = 8000)

# OR launch the multi-study browser (card-grid index of every
# <studyID>/<studyID>.scminer.h5 found under the root)
run_browser("path/to/data", port = 8000)
```

For programmatic pipelines that already have matrices in memory, use
`prepare_study_from_eset()` (Biobase ExpressionSet input) or
`prepare_study_data()` (plain R structures).

## Tutorial: from zero to a multi-study browser

### A. Prepare one study from a YAML config

1. Copy the annotated template that ships with the package:

   ```r
   file.copy(
     system.file("extdata", "example_config.yml", package = "scminerViewer"),
     "config.yml"
   )
   ```

2. Edit `config.yml` — fill in `output:` (the **root** containing all
   your studies, e.g. `"data"`), the `study.*` block (ID/abbr/titles),
   and the `input.*` paths to your scMINER `ExpressionSet` RDS files.

3. Run:

   ```r
   library(scminerViewer)
   prepare_study("config.yml")
   ```

   This writes `data/<studyID>/<studyID>.scminer.h5` plus the full
   graph layout and the shard tree inside `data/<studyID>/`.

### B. Or migrate an existing on-disk study (no RDS needed)

If you already have the graph-layout TSVs + shard tree on disk (e.g.
under `data/example/`), build the bundle directly without invoking the
Biobase pipeline:

```r
study <- read_graph_study(data_dir = "data/example", study_id = "2327")
dir.create("data/2327", recursive = TRUE, showWarnings = FALSE)
write_bundle(
  bundle_path        = "data/2327/2327.scminer.h5",
  meta               = study$meta,
  cells              = study$cells,
  clusters           = study$clusters,
  genes              = study$genes,
  expression_genes   = study$expression_genes,
  activity_tf_genes  = study$activity_tf_genes,
  activity_sig_genes = study$activity_sig_genes,
  network_tf         = study$network_tf,
  network_sig        = study$network_sig,
  overwrite          = TRUE
)
```

If you also want the shard tree co-located with the new bundle, move
it once with `mv` / `file.rename()`; otherwise pass `shard_dir =` to
`run_app()` / `run_browser()` (see C).

### C. Browse one or many studies

```r
# Single study — bundle and shards co-located
run_app("data/2327/2327.scminer.h5", port = 8000)

# Single study — shards live somewhere else
run_app("data/2327/2327.scminer.h5",
        shard_dir = "data/example", port = 8000)

# Multi-study card-grid index of every <studyID>/<studyID>.scminer.h5
# found under the root; click a card to drill into that study's viewer
run_browser("data", port = 8000)

# Multi-study with external shard root (single shard_dir is shared by
# every study under root_dir — useful while migrating)
run_browser("data", shard_dir = "data/example", port = 8000)
```

`discover_studies("data")` returns a data.frame of every bundle the
browser would find; useful for scripting or sanity-checking the layout
before launching the app.

### D. Adding a second study

Just point `prepare_study()` (or another migration) at the same root
with a new `study.ID`. Each call creates a sibling `<sid>/` folder:

```
data/
├── 2327/                      # first study
│   └── 2327.scminer.h5
└── 9999/                      # second study
    └── 9999.scminer.h5
```

`run_browser("data")` now shows two cards. No re-build of existing
studies is required — `discover_studies()` re-scans on each launch.

## Exported API

### Orchestrators (do everything in one call)

| Function | Purpose |
| --- | --- |
| `prepare_study(config_path, emit, verbose)`                    | YAML-driven entry. Requires `yaml` + `Biobase`. |
| `prepare_study_from_eset(out_dir, expression_eset, ...)`       | Accepts a Biobase ExpressionSet. Splits activity rows into TF/SIG by `_TF`/`_SIG` row suffix. Requires `Biobase`. |
| `prepare_study_data(out_dir, meta, cells, clusters, genes, expression, ..., default_genes, emit, verbose)` | Lowest-level orchestrator — plain R structures only. |

### Staged helpers (granular control / debugging)

Call these directly when you want to inspect an intermediate before
committing to disk. `prepare_study()` and `prepare_study_from_eset()`
are thin wrappers over them.

| Function | Purpose |
| --- | --- |
| `load_study_config(config_path)`                               | Parse + validate the YAML; fill in defaults. Pure parsing — does not touch the RDS / TSV files referenced in `input.*`. |
| `extract_cells(eset, cell_id_col, cell_type_col, cell_group_col, coordinate_col)` | `pData(eset)` → cells data.frame. Errors loudly on missing columns. |
| `extract_genes(eset, gene_symbol_col)`                         | `fData(eset)` → character vector of gene symbols. |
| `extract_expression(eset, genes = NULL)`                       | `exprs(eset)` → sparse `Matrix` (genes × cells). Sets rownames from `genes` if provided. |
| `extract_activity(activity_eset, master_genes)`                | Splits rows by `_TF`/`_SIG` suffix; returns `list(tf, sig)` reindexed to `master_genes`. Either element may be `NULL`. |
| `read_networks(path)`                                          | Parse a scMINER networks TSV; returns `list(tf, sig)` data.frames. |

### Bundle + Shiny

| Function | Purpose |
| --- | --- |
| `write_bundle(bundle_path, meta, cells, clusters, genes, expression_genes, activity_tf_genes, activity_sig_genes, default_genes, network_tf, network_sig, overwrite)` | Write the `.scminer.h5` bundle. Indexes + metadata only; values stay on disk. |
| `read_graph_study(data_dir, study_id)`                         | Reconstruct study inputs from the on-disk graph layout. Per-matrix gene **indexes** come from the manifest CSVs; shard values are never loaded eagerly. |
| `load_study(bundle_path, shard_dir = NULL)`                    | Read a `.scminer.h5` into an S3 `scminer_study` list. `shard_dir` defaults to `dirname(bundle_path)`. |
| `gene_values(study, gene, relationship)`                       | Lazily read one gene's row from the shard tree (`Express_normalized`, `Activity_tf`, or `Activity_sig`). Aligned to `study$cells$cellID`; cached per gene. |
| `run_app(bundle_path, host, port, launch_browser, ...)`        | Launch the single-study Shiny app for a bundle. |
| `build_app(bundle_path)`                                       | Build a single-study `shiny.appobj` without launching it. |
| `discover_studies(root_dir)`                                   | Return one row per `<studyID>/<studyID>.scminer.h5` bundle found under `root_dir` (data.frame with meta + n_cells / n_genes / paths). |
| `run_browser(root_dir, shard_dir = NULL, host, port, launch_browser, ...)` | Launch the **multi-study** browser. Card-grid index page; click → `?study=<id>` opens that study's viewer with a "← Back" link. `shard_dir` overrides where shards live (each bundle defaults to `dirname(bundle_path)`). |
| `build_browser(root_dir, shard_dir = NULL)`                    | Build a multi-study `shiny.appobj` without launching it. |

### Debug / granular workflow

The orchestrators do `parse → load → extract → write` in one call. To
inspect or override any step, run them manually:

```r
library(scminerViewer)

# 1) Parse the config without touching the RDS files
cfg <- load_study_config("config.yml")
str(cfg)            # see exactly what defaults got applied

# 2) Load the inputs (slow — gigabyte-sized RDS files)
expr_eset <- readRDS(cfg$input$expression)
act_eset  <- if (!is.null(cfg$input$activity)) readRDS(cfg$input$activity)

# 3) Extract each piece independently — every helper returns a plain R
#    structure you can print, summary(), or modify before continuing
cells <- extract_cells(
  expr_eset,
  cell_id_col    = cfg$cellID,
  cell_type_col  = cfg$cellType,
  cell_group_col = cfg$cellGroup,
  coordinate_col = cfg$coordinate
)
genes <- extract_genes(expr_eset, gene_symbol_col = cfg$geneSymbol)
expr  <- extract_expression(expr_eset, genes = genes)
act   <- if (!is.null(act_eset)) extract_activity(act_eset, genes)
        else list(tf = NULL, sig = NULL)
nets  <- if (!is.null(cfg$input$networks))
          read_networks(cfg$input$networks)
        else list(tf = NULL, sig = NULL)

# 4) Inspect / tweak before writing
head(cells)
dim(expr)
table(sapply(nets, NROW))

# 5) Hand the parts to prepare_study_data — it writes the graph layout,
#    the bundle, or both
prepare_study_data(
  out_dir      = cfg$output,
  meta         = list(
    studyID    = cfg$study$ID,    studyAbbr  = cfg$study$studyAbbr,
    longTitle  = cfg$study$longTitle, shortTitle = cfg$study$shortTitle,
    species    = cfg$species,     coordinate = cfg$coordinate
  ),
  cells = cells, genes = genes,
  expression = expr,
  activity_tf = act$tf, activity_sig = act$sig,
  network_tf  = nets$tf, network_sig = nets$sig,
  emit = c("graph", "bundle"),
  verbose = TRUE
)
```

Use `emit = "bundle"` (or `"graph"`) to write only one of the two.

## Config YAML

`prepare_study(config_path)` reads a YAML file pointing at a scMINER
ExpressionSet RDS (and optional activity + networks files) and writes
both the graph layout and the bundle. Two examples are available:

- **Annotated template** shipped inside the package — every field
  documented inline. Copy it with `file.copy(system.file(...))` below.
- **Concrete runnable config** at
  [`data/example/2327.yml`](../data/example/2327.yml) keyed to the
  actual 2327 (Tex) study in this repo. Adjust the `input.*` paths to
  point at your scMINER RDS / TSV files and it runs end-to-end.

Copy the annotated template:

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

169 tests across `test-bundle-roundtrip.R`, `test-graph-read.R`,
`test-shard-reader.R`, `test-prepare-study.R`, `test-staged-helpers.R`,
`test-browser.R`, and `test-app.R` (plus shared `helper-fixtures.R`).
Set `SCMINER_DATA_DIR=/path/to/data` to enable the integration tests
that read real graph-layout outputs.

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
