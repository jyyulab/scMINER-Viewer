# scMINER Viewer — Applications Note

Short paper modelled after Khatamian *et al.* (*Bioinformatics*
35:2165, 2019). Source manuscript, reproducible benchmarks, generated
figures + tables, and the LSF-driven HPC pipeline for the live scMINER
Portal studies.

## Directory layout

```
paper/
├── README.md                       (this file)
├── scminer-viewer.md               manuscript draft (Markdown source)
├── Manuscript-V1-HZ.docx           Word-formatted submission draft
├── Figures.pptx                    figure-assembly working file
├── code/                           all reproducibility code
│   ├── benchmarks/                   figure + table generators
│   │   ├── figure_architecture.R       Figure 1
│   │   ├── methods.R                   bench primitives (sourced by figures.R)
│   │   ├── figures.R                   synthetic-sweep → Figure 2 panels + TSVs
│   │   ├── figure_portal.R             real-portal Figure 3 panels
│   │   ├── figure_portal_compare.R     Supplemental: expr-only + full↔expr-only
│   │   ├── tables.R                    Table 1 + Table 2
│   │   ├── tables_compare.R            Supplemental tables (S + delta)
│   │   └── run_figure1.sh              one-command wrapper for Figure 2
│   ├── portal/                       29-study HPC pipeline
│   │   ├── portal_studies.R            worker: prepare → load → fetch, per YAML
│   │   ├── portal_studies.sh           local wrapper
│   │   ├── portal_studies.bsub         LSF script (single + array body)
│   │   ├── portal_studies_hpc.sh       array driver + dependent merge
│   │   ├── portal_studies_single.sh    sequential bsub driver
│   │   ├── portal_studies_compare.sh   29 studies × 2 modes (full / expr-only)
│   │   ├── portal_compare.R            join per-mode TSVs into compare table
│   │   ├── portal_merge.R              concatenate per-study TSVs
│   │   ├── sparseify_eset.R / .sh      one-time dense→dgCMatrix converter
│   │   └── example.yaml                annotated YAML template
│   └── configs/                      29 YAML configs, one per study
│       ├── 2202.yaml, 2203.yaml, ..., 2342.yaml
│       └── (paths point at HPC data roots)
├── metrics/                        generated TSVs (output of code/)
│   ├── benchmarks/                   synthetic-sweep results
│   │   ├── bundle_scaling.tsv          per-replicate (140 rows)
│   │   ├── bundle_scaling_summary.tsv  per-config mean / sd / SE
│   │   ├── discover_scaling.tsv
│   │   └── discover_scaling_summary.tsv
│   └── portal_studies/               per-study real-data results
│       ├── portal_studies_<id>_full.tsv               (29 studies)
│       ├── portal_studies_<id>_expression-only.tsv    (29 studies)
│       ├── portal_studies_compare.tsv                 joined wide table
│       └── portal_studies_summary.tsv                 status counts
├── figures/                        rendered PDF + PNG, one pair per panel
│   ├── architecture.{pdf,png}            Figure 1
│   ├── figure2.{pdf,png}                 Figure 2 (combined)
│   ├── figure2_{A..F}_*.{pdf,png}        Figure 2 panels
│   ├── figure3.{pdf,png}                 Figure 3 (combined)
│   ├── figure3_{A..F}_*.{pdf,png}        Figure 3 panels
│   └── compare/                          Supplemental + paired-comparison artefacts
│       ├── figureS_expr_only.{pdf,png}      (combined)
│       ├── figureS_expr_only_{A..D}.{pdf,png}
│       ├── compare.{pdf,png}                (combined)
│       └── compare_{A..D}_*.{pdf,png}
└── tables/                         rendered tables (one .md + one .tsv each)
    ├── figure2_scaling.{md,tsv}          Table 2 — synthetic sweep
    ├── figure3_portal_studies.{md,tsv}   Table 1 — 13 TF/sig-eligible portal studies
    ├── tableS_expr_only.{md,tsv}         Supplemental — expression-only baseline
    └── compare_delta.{md,tsv}            Supplemental — paired full vs expression-only delta
```

## Contents at a glance

| Path | What it is |
| --- | --- |
| [`scminer-viewer.md`](scminer-viewer.md)                                                | Manuscript draft (~2 pages, Markdown). |
| [`Manuscript-V1-HZ.docx`](Manuscript-V1-HZ.docx)                                        | Word-formatted draft for review / submission. |
| [`Figures.pptx`](Figures.pptx)                                                          | PowerPoint workbook used for final figure assembly. |
| [`code/benchmarks/figure_architecture.R`](code/benchmarks/figure_architecture.R)        | Figure 1 (data-flow diagram). No data deps; ~ 1 s. |
| [`code/benchmarks/methods.R`](code/benchmarks/methods.R)                                | Synthetic study generator + bench primitives. Sourced by `figures.R`. |
| [`code/benchmarks/figures.R`](code/benchmarks/figures.R)                                | 7 × 4 synthetic sweep × 5 reps + 2327 row + discover sweep. Writes 6 Figure 2 panels + combined `figure2.{pdf,png}` + scaling TSVs. |
| [`code/benchmarks/figure_portal.R`](code/benchmarks/figure_portal.R)                    | Renders 6 Figure 3 panels + combined `figure3.{pdf,png}` for the 13 TF/sig-eligible studies. |
| [`code/benchmarks/figure_portal_compare.R`](code/benchmarks/figure_portal_compare.R)    | Supplemental: expression-only panels + paired full↔expr-only comparison. |
| [`code/benchmarks/tables.R`](code/benchmarks/tables.R)                                  | Generates Table 1 (`figure3_portal_studies`) + Table 2 (`figure2_scaling`). |
| [`code/benchmarks/tables_compare.R`](code/benchmarks/tables_compare.R)                  | Generates Supplemental Table (`tableS_expr_only`) + paired delta (`compare_delta`). |
| [`code/benchmarks/run_figure1.sh`](code/benchmarks/run_figure1.sh)                      | One-command Figure 2 driver: runs `figures.R` then `tables.R`. Local by default; pass `--hpc --mem ... --wall ...` to bsub. |
| [`code/configs/`](code/configs/)                                                        | 29 YAML configs (`<studyID>.yaml`) pointing at HPC data. |
| [`code/portal/portal_studies.R`](code/portal/portal_studies.R)                          | HPC worker: per-YAML `prepare_study_from_eset → load_study → gene_values`, writes per-study metrics. |
| [`code/portal/portal_studies_hpc.sh`](code/portal/portal_studies_hpc.sh)                | One-command HPC driver — submits a parallel job array + dependent merge. |
| [`code/portal/portal_studies_single.sh`](code/portal/portal_studies_single.sh)          | Sequential single-bsub driver (configurable mem / cores / wall). |
| [`code/portal/portal_studies_compare.sh`](code/portal/portal_studies_compare.sh)        | With-vs-without-TFsig comparison (29 studies × 2 modes). |
| [`code/portal/portal_compare.R`](code/portal/portal_compare.R)                          | Joins per-mode TSVs into one wide compare table with `delta_*` columns. |
| [`code/portal/portal_merge.R`](code/portal/portal_merge.R)                              | Concatenates per-study TSVs. |
| [`code/portal/sparseify_eset.{R,sh}`](code/portal/)                                     | One-time dense→`dgCMatrix` converter for OOM-prone source rds files. |
| [`code/portal/example.yaml`](code/portal/example.yaml)                                  | Annotated YAML schema. |
| [`metrics/benchmarks/`](metrics/benchmarks/)                                            | Synthetic-sweep TSVs (replicate-level + summary). |
| [`metrics/portal_studies/`](metrics/portal_studies/)                                    | Per-study × per-mode TSVs + joined `portal_studies_compare.tsv`. |
| [`figures/`](figures/)                                                                  | Flat tree of rendered figure PDFs + PNGs (one pair per panel + one combined per figure). |
| [`tables/`](tables/)                                                                    | Rendered Markdown + TSV tables. |

## Quick reference: regenerate figures and tables

Run from the **project root** with `scminerViewer` installed. The
portal figure / Table 1 path expects the per-mode TSVs under
`paper/metrics/portal_studies/` (or the path the scripts are pointed
at via `--metrics-dir`).

```sh
# --- Figure 1 (architecture diagram) -------------------------------
Rscript paper/code/benchmarks/figure_architecture.R
# Writes: paper/figures/architecture.{pdf,png}                (~1 s)

# --- Figure 2 (synthetic-sweep, 28 configs × 5 reps) ----------------
# Local (~60–90 min on a laptop):
./paper/code/benchmarks/run_figure1.sh
# Or via LSF on HPC (recommended for the full sweep):
./paper/code/benchmarks/run_figure1.sh --hpc --mem 64000 --wall 6:00
# Or the R steps directly (skip the wrapper):
Rscript paper/code/benchmarks/figures.R
Rscript paper/code/benchmarks/tables.R
# Writes: paper/figures/figure2_{A..F}_*.{pdf,png} + figure2.{pdf,png}
#         paper/metrics/benchmarks/{bundle,discover}_scaling{,_summary}.tsv
#         paper/tables/figure2_scaling.{md,tsv}

# --- Figure 3 + Table 1 (13 TF/sig-eligible portal studies) ---------
Rscript paper/code/benchmarks/figure_portal.R
Rscript paper/code/benchmarks/tables.R
# Writes: paper/figures/figure3_{A..F}_*.{pdf,png} + figure3.{pdf,png}
#         paper/tables/figure3_portal_studies.{md,tsv}

# --- Supplemental (expression-only + paired delta) ------------------
Rscript paper/code/benchmarks/figure_portal_compare.R
Rscript paper/code/benchmarks/tables_compare.R
# Writes: Supplemental panels + figureS_expr_only.{pdf,png}
#         paper/tables/tableS_expr_only.{md,tsv}
#         paper/tables/compare_delta.{md,tsv}
```

**Pointing the portal figure / table at a different snapshot** (e.g.
two HPC runs you want to compare):

```sh
Rscript paper/code/benchmarks/figure_portal.R \
    --metrics-dir paper/metrics/05162026 \
    --figures-dir paper/figures/05162026
```

## Synthetic-sweep benchmark (Figure 2)

The wrapper runs `figures.R` → `tables.R` and tees output to
`paper/logs/figure1_<timestamp>.log`:

```sh
# Local (laptop / HPC interactive node)
./paper/code/benchmarks/run_figure1.sh

# HPC (single bsub job; recommended for the full 28-config sweep
# since 10K × 10K configs may need ~16–32 GB):
./paper/code/benchmarks/run_figure1.sh --hpc --mem 64000 --wall 6:00

# Inspect the bsub command without submitting:
./paper/code/benchmarks/run_figure1.sh --hpc --dry-run

# R steps directly:
Rscript paper/code/benchmarks/figures.R
Rscript paper/code/benchmarks/tables.R
```

Grid: `n_cells ∈ {500, 1000, 2000, 4000, 6000, 8000, 10000}` ×
`n_genes ∈ {2000, 5000, 8000, 10000}` (28 configs). Each
configuration runs `N_REPS = 5` times with independent seeds; mean ±
SE drives Figure 2 error bars and Table 2 cells. Total wall time:
~ 60–90 min on a laptop (140 runs), ~ 30 s for the multi-study
`discover_studies()` sweep (5 × 6 = 30 runs), < 1 s for the real-2327
row. Drop `N_REPS` to `1` in `figures.R` for an iteration cycle.

`run_figure1.sh` flags (HPC mode only):

| Flag | Default | Meaning |
| --- | --- | --- |
| `--hpc` | — | Submit via `bsub` (otherwise: run in-shell). |
| `--mem <MB>` | `32000` | Memory request. |
| `--cores <n>` | `4` | Cores. |
| `--wall <hh:mm>` | `4:00` | Wall-clock limit. |
| `--queue <name>` | `standard` | LSF queue. |
| `--project <id>` | `scminer` | LSF charge code. |
| `--dry-run` | — | Print the bsub command without submitting. |

Each panel is rendered as its own `{pdf,png}` pair:

| File stem | What it shows |
| --- | --- |
| `figure2_A_size_vs_cells`     | Bundle + shard tree size vs cells, faceted by gene count |
| `figure2_B_load_latency`      | `load_study()` cold-start latency |
| `figure2_C_fetch_latency`     | `gene_values()` median fetch latency (bar to max) |
| `figure2_D_discover_scaling`  | `discover_studies()` scaling on a multi-study root |
| `figure2_E_prepare_time`      | `prepare_study_data()` wall time |
| `figure2_F_peak_memory`       | `prepare_study_data()` peak resident memory |

Edit `SCALING_GRID` / `DISCOVER_GRID` in `figures.R` to widen / shrink
the sweep.

## Portal benchmark (29 real studies on HPC)

`paper/code/configs/` holds **29 YAML configs**, one per study,
named `<study.ID>.yaml`. Each config tells `portal_studies.R` exactly
where its `expression.rds`, `activity.rds`, and `networks.txt` live.

### Config schema (`paper/code/configs/2327.yaml`)

```yaml
study:
    ID:         "2327"
    studyAbbr:  tex
    longTitle:  Single-cell RNA sequencing of progenitor, intermediate and terminally exhausted CD8+ T cells from tumors and chronic viral infection
    shortTitle: Tex
species:    Mus musculus
coordinate: UMAP

# Column names inside the ExpressionSets
cellType:    CellGroup
cellGroup:   CellGroup
geneSymbol:  GeneSymbol

input:
    expression: /research_jude/.../Studies/Tex/expression.rds
    activity:   /research_jude/.../Studies/Tex/activity.rds
    networks:   /research_jude/.../Studies/Tex/networks.txt
output:        /research_jude/.../scMINER_Portal/scMINERViewerMetrics
```

* Required keys: `study.ID`, `input.expression`.
* Optional with sensible defaults: `cellID` (`cellID`), `cellType`
  (`cellGroup`), `cellGroup` (= `cellType`), `geneSymbol`
  (`geneSymbol`), `coordinate` (`UMAP`), `cluster_palette` (`npg`).
* Path resolution: paths inside `input.*` may be absolute (most are)
  or relative — relative paths resolve against `input_root` if set,
  else the YAML's own parent dir.
* `output:` is the directory where each study's bundle, shard tree,
  and graph files land. `portal_studies.R` resolves it in this order:
    1. `--output-root <dir>` CLI flag (single root for every study)
    2. The YAML's `output:` value (per-study)
    3. `--scratch <dir>` (default: ephemeral `tempfile()`)

  With the unified `output: /research/.../scMINERViewerMetrics` set
  in every config, a study with `study.ID: 2327` ends up at:
  ```
  /research/.../scMINERViewerMetrics/2327/2327.scminer.h5
  /research/.../scMINERViewerMetrics/2327/expression_files/...
  /research/.../scMINERViewerMetrics/2327/activity_files/...
  /research/.../scMINERViewerMetrics/2327/graph_files/...
  ```

### One-command HPC run

```sh
./paper/code/portal/portal_studies_hpc.sh --configs-dir paper/code/configs
```

The driver enumerates every `*.yaml` in `--configs-dir`, writes the
manifest to `paper/hpc/manifest_<timestamp>.txt` (gitignored), submits
a parallel LSF job array (one task per YAML), and a dependent merge
job that concatenates per-study TSVs.

Driver flags:

| Flag | Default | Meaning |
| --- | --- | --- |
| `--configs-dir <dir>` | — | Folder of YAML configs. |
| `--studies-root <dir>` | — | Alternative: folder of per-study subfolders. |
| `--queue "<q1 q2 ...>"` | `"standard priority"` | Space-separated LSF queue list. Quote when listing >1. |
| `--mem <MB>` | `32000` | Memory per task. |
| `--cores <n>` | `4` | Cores per task. |
| `--wall <hh:mm>` | `6:00` | Wall-clock limit per task. |
| `--project <id>` | `scminer` | LSF charge code. |
| `--dry-run` | — | Print the bsub commands without submitting. |

Monitor & collect:

```sh
bjobs -A                          # array summary
bjobs -J scminer_portal_*         # individual tasks
tail -f paper/logs/portal_studies_*_*.out
ls   paper/metrics/portal_studies/portal_studies_*.tsv    # per-study TSVs
```

**Aggregate + render real-data figures** once the per-study TSVs
exist:

```sh
Rscript paper/code/benchmarks/figure_portal.R
Rscript paper/code/benchmarks/tables.R
```

Each panel is its own file so it can be placed independently:

| File stem | What it shows |
| --- | --- |
| `figure3_A_size_vs_cells` | Bundle + shard tree size vs cell count (log-log) |
| `figure3_B_prepare_time`  | `prepare_study()` wall time vs `n_cells × n_genes` |
| `figure3_C_peak_memory`   | Peak memory vs total input size |
| `figure3_D_load_latency`  | `load_study()` cold-start latency vs bundle size |
| `figure3_E_fetch_latency` | `gene_values()` median fetch latency vs gene count |
| `figure3_F_size_ratio`    | Output : input compression ratio per study |

`figure_portal.R` reads per-mode `*_full.tsv` files (filtered to
`net_tf_edges + net_sig_edges > 0`, so Figure 3 / Table 1 only contain
studies whose YAML actually supplies a `networks.txt`). Studies
dropped by this filter still appear in the Supplemental
expression-only artefacts produced by `figure_portal_compare.R` /
`tables_compare.R`.

### With-vs-without-TFsig comparison (29 studies × 2 modes)

[`portal_studies_compare.sh`](code/portal/portal_studies_compare.sh)
runs the benchmark twice per study so you can quantify what TF/sig
data costs in `prepare_study()` time, peak memory, cold/warm
`load_study()` latency, and on-disk bundle size:

* **`expression-only`** — `activity.rds` and `networks.txt` are
  skipped; `prepare_study_from_eset()` receives `activity_eset = NULL`
  and `networks_path = NULL`.
* **`full`** — historical default: expression + activity + networks.

The two modes share each YAML's `output:` directory, so concurrent
execution would race on the same `<output>/<studyID>/` path. The
driver offers two execution modes:

* **Sequential (default)** — submits `expression-only` first, then
  `full` with `-w "ended(expression-only)"`. Safe on a single shared
  output folder.
* **`--parallel`** — both 29-task arrays fire concurrently. Requires
  disjoint output roots, threaded through to `portal_studies.R`'s
  `--output-root` flag.

```sh
# Sequential (default) -- one shared output root, safest:
./paper/code/portal/portal_studies_compare.sh \
    --configs-dir paper/code/configs \
    --mem 64000 --wall 8:00 --queue "standard priority"

# Parallel -- two arrays run concurrently, ~half the wall time:
./paper/code/portal/portal_studies_compare.sh \
    --configs-dir paper/code/configs \
    --parallel \
    --expr-only-root /research_jude/.../scMINERViewerMetrics/expression-only \
    --full-root      /research_jude/.../scMINERViewerMetrics/full \
    --mem 64000 --wall 8:00 --queue "standard priority"

# Re-run a subset:
./paper/code/portal/portal_studies_compare.sh \
    --configs-dir paper/code/configs --only 2317,2327
```

Driver flags:

| Flag | Default | Meaning |
| --- | --- | --- |
| `--configs-dir <dir>` | `paper/code/configs` | Folder of YAML configs. |
| `--parallel` | — | Fire both arrays concurrently. Requires `--expr-only-root` + `--full-root`. |
| `--expr-only-root <dir>` | — | Output root for the `expression-only` mode. |
| `--full-root <dir>` | — | Output root for the `full` mode. |
| `--queue "<q1 q2 ...>"` | `"standard priority"` | LSF queue list. |
| `--mem <MB>` | `64000` | Memory per task. |
| `--cores <n>` | `4` | Cores per task. |
| `--wall <hh:mm>` | `8:00` | Wall-clock limit per task. |
| `--project <id>` | `scminer` | LSF charge code. |
| `--only <csv>` | (all) | Comma-separated study IDs to restrict to. |
| `--dry-run` | — | Print the bsub commands without submitting. |

Per-mode TSVs land at:

```
paper/metrics/portal_studies/portal_studies_<studyID>_full.tsv
paper/metrics/portal_studies/portal_studies_<studyID>_expression-only.tsv
```

The dependent compare job produces the joined wide table at
`paper/metrics/portal_studies/portal_studies_compare.tsv`, with these
metric triples (`*_full`, `*_expr_only`, `delta_*` where the delta is
`full − expr_only`, so positive = TF/sig adds cost):

| Metric | Columns |
| --- | --- |
| Wall time | `prepare_seconds_full`, `prepare_seconds_expr_only`, `delta_prepare_seconds` |
| Peak memory (MB) | `prepare_peak_mb_full`, `prepare_peak_mb_expr_only`, `delta_prepare_peak_mb` |
| Cold load | `load_seconds_full`, `load_seconds_expr_only`, `delta_load_seconds` |
| Warm load | `load_seconds_warm_full`, `load_seconds_warm_expr_only`, `delta_load_seconds_warm` |
| Bundle size | `bundle_bytes_full`, `bundle_bytes_expr_only`, `delta_bundle_bytes` |
| Total output | `total_output_bytes_full`, `total_output_bytes_expr_only`, `delta_total_output_bytes` |

Identity / sanity columns: `studyID`, `n_cells_full`, `n_genes_full`,
`n_clusters_full`, `net_tf_edges_full`, `net_sig_edges_full`,
`status_full`, `status_expr_only`, `note_full`, `note_expr_only`.

Render the figures and tables that consume this compare TSV:

```sh
Rscript paper/code/benchmarks/figure_portal_compare.R
Rscript paper/code/benchmarks/tables_compare.R
```

Re-run just the aggregator (no LSF):

```sh
Rscript paper/code/portal/portal_compare.R \
    --expr-only-glob "paper/metrics/portal_studies/portal_studies_*_expression-only.tsv" \
    --full-glob      "paper/metrics/portal_studies/portal_studies_*_full.tsv" \
    --out            "paper/metrics/portal_studies/portal_studies_compare.tsv"
```

Monitor:

```sh
bjobs -A
bjobs -J 'scminer_cmp_*_expression-only_*'   # mode 1
bjobs -J 'scminer_cmp_*_full_*'              # mode 2
bjobs -J 'scminer_cmp_*_compare'             # join
tail -f paper/logs/portal_studies_expression-only_*_*.out
tail -f paper/logs/portal_studies_full_*_*.out
```

Notes:

* The dependency is `ended(...)` (not `done(...)`), so a per-study
  failure in one mode does **not** block the other 28. Failed
  per-mode tasks leave no TSV; the compare aggregator handles missing
  rows with `NA` deltas.
* For 2317 / Covid650k and 2333 / ATRT, prime the dense
  `expression.rds` with [`sparseify_eset.sh`](code/portal/sparseify_eset.sh)
  and update the YAML before submitting — see *Priming large studies*
  below.

### Tables

`paper/code/benchmarks/tables.R` regenerates the two main manuscript
tables from the TSVs:

```sh
Rscript paper/code/benchmarks/tables.R
```

Outputs land under `paper/tables/`:

| Table | Source | Output |
| --- | --- | --- |
| Table 1 — per-study metrics for the 13 TF/sig-eligible portal studies (backs Figure 3) | `paper/metrics/portal_studies/portal_studies_*_full.tsv` (filtered to `net_tf_edges + net_sig_edges > 0`) | `paper/tables/figure3_portal_studies.{md,tsv}` |
| Table 2 — synthetic-sweep grid (backs Figure 2) | `paper/metrics/benchmarks/bundle_scaling.tsv` | `paper/tables/figure2_scaling.{md,tsv}` |
| Supplemental — expression-only baseline | `paper/metrics/portal_studies/portal_studies_compare.tsv` (status_expr_only == `ok`) | `paper/tables/tableS_expr_only.{md,tsv}` — run `Rscript paper/code/benchmarks/tables_compare.R` |
| Paired delta — full vs expression-only per metric | same compare TSV | `paper/tables/compare_delta.{md,tsv}` — same script |

The Markdown files ship a leading caption block + a GFM table. The
TSV versions are clean for any downstream post-processing.

### Local / sequential modes (development & reruns)

```sh
# Local interactive — one or many studies on the current machine:
./paper/code/portal/portal_studies.sh --configs-dir paper/code/configs
./paper/code/portal/portal_studies.sh --config       paper/code/configs/2327.yaml
./paper/code/portal/portal_studies.sh --configs-dir  paper/code/configs --only 2327,2326

# Single bsub (sequential walk, no array) -- mem / wall / cores knobs:
./paper/code/portal/portal_studies_single.sh --configs-dir paper/code/configs
./paper/code/portal/portal_studies_single.sh --configs-dir paper/code/configs --mem 64000 --wall 12:00
./paper/code/portal/portal_studies_single.sh --configs-dir paper/code/configs --only 2327,2326

# Single bsub (minimal, env-only -- defaults to mem=32000):
CONFIGS_DIR=$(pwd)/paper/code/configs bsub < paper/code/portal/portal_studies.bsub
ONLY=2327,2326 CONFIGS_DIR=$(pwd)/paper/code/configs bsub < paper/code/portal/portal_studies.bsub

# Single bsub with mem override (bsub flag overrides #BSUB directive):
CONFIGS_DIR=$(pwd)/paper/code/configs bsub -R "rusage[mem=64000]" -M 64000 -W 12:00 \
    < paper/code/portal/portal_studies.bsub
```

`portal_studies.R` flags:

| Flag | Default | Meaning |
| --- | --- | --- |
| `--configs-dir <dir>` | — | Process every YAML in the directory. |
| `--config <yaml>`     | — | Process exactly one YAML (used by array tasks). |
| `--studies-root <dir>` | — | Walk per-study subfolders (legacy mode). |
| `--only a,b,c`        | (all) | CSV of study IDs to restrict to. |
| `--mode <name>`       | `full` | `full` reads expression + activity + networks. `expression-only` skips activity and networks; used by `portal_studies_compare.sh` to isolate the cost of TF/sig data. Each mode writes into its own per-study subdir (`<out>/<studyID>/<mode>/`). |
| `--out <tsv>`         | `paper/metrics/portal_studies.tsv` | Output TSV path. (Older scripts may still target this legacy flat path.) |
| `--output-root <dir>` | — | Force every study's bundle / shards / graph files into `<dir>/<studyID>/`, overriding the YAML's `output:`. |
| `--scratch <dir>`     | tempfile | Fallback when neither `--output-root` nor the YAML's `output:` is set. |
| `--quiet`             | — | Suppress progress messages. |

### Metrics columns (per-study TSV)

| Group | Columns |
| --- | --- |
| Identification | `studyID`, `mode`, `n_cells`, `n_genes`, `n_clusters`, `out_dir` |
| **Inputs** | `expr_input_bytes`, `act_input_bytes`, `net_input_bytes`, `total_input_bytes` |
| **Outputs** | `bundle_bytes`, `shard_bytes`, `graph_bytes`, `total_output_bytes` |
| Wall time & memory | `prepare_seconds`, `prepare_peak_mb`, `load_seconds`, `load_seconds_warm` |
| Gene fetch | `fetch_median`, `fetch_mean`, `fetch_max`, `n_fetched` |
| Networks | `net_tf_edges`, `net_sig_edges` |
| Status | `status` (`ok` / `skipped` / `error-too-large` / `error`), `note` |

Status values:

| `status` | When | TSV columns |
| --- | --- | --- |
| `ok` | Benchmark completed | All metric columns filled |
| `skipped` | YAML uses the legacy `input.genes` (raw matrix + `genes.csv`) layout that `prepare_study_from_eset()` can't consume | Metrics NA; `note` carries reason |
| `error-too-large` | Run-time failure attributable to R's `dgCMatrix` 2^31 − 1 nnz cap (typically inside `.reindex_rows()` for very large activity matrices) | Metrics NA; `note` carries the offending message |
| `error` | Any other failure | Metrics NA; `note` carries the error message |

The bsub task always exits 0 regardless of status (skip / cap / error
are all "expected outcomes"), so the merge job's `ended(array)`
dependency always fires.

### Resource sizing

Empirical numbers from the 16 May 2026 run over all 28 OK studies
(see `paper/metrics/portal_studies/` and
`paper/tables/figure3_portal_studies.md` for per-study rows):

| Cell range | Representative studies | `prepare_seconds` | Peak Mb | `--mem` recommendation |
| --- | --- | --- | --- | --- |
| ≤ 1 k cells (ground-truth set) | Buettner / Chung / Yan / Zeisel / Klein / Kolod / Goolam / Pollen / Usoskin | 5–20 min | 350–500 MB | `--mem 8000` |
| 1 k – 10 k cells | Tex (2327), Tregs (2326), PBMC14k (2325), Hepatoblastoma (2318), iCCA (2319), GSE155446 (2332), Glutamatergic (2331) | 10–60 min | 0.5–3.5 GB | `--mem 16000` |
| 10 k – 50 k cells | Bcell_Dev (2203, 54 k cells × 33 k genes) | 1–2 h | ~ 36 GB | `--mem 64000` |
| 50 k – 100 k cells | study 2339 (50 k), study 2342 (57 k), HMC76k (76 k), Covid97k (97 k), study 2338 (96 k) | 1–17 h | 17–49 GB | `--mem 64000` (Covid97k) → `--mem 96000` (2338) |
| 100 k+ cells | study 2341 (105 k), **ATRT / 2333 (138 k cells × 18 k genes)** | 2.5–16 h | 46–49 GB | `--mem 96000 --wall 24:00` |
| 650 k cells | **Covid650k / 2317 (647 k cells × 25 k genes)** — completes after sparseification (see *Priming large studies* below) | ~15 h | 121 GB | `--mem 130000 --wall 24:00` |

Heuristic for a new study: peak memory ≈ 1–3× the on-disk
`expression.rds` size. If the rds is stored densely (e.g. ATRT's
`dgeMatrix`), the multiplier can climb to 10×; run
`paper/code/portal/sparseify_eset.sh --verify-only --in …` to check
the matrix class before requesting memory.

**Total wall-time for a full 28-task array**: roughly equal to the
single longest study (~17 h on study 2341, ~15 h on Covid650k), since
the array runs in parallel across hosts and queues.

### Priming large studies (dense-backed `ExpressionSet`s)

A handful of portal studies ship `expression.rds` files whose
`exprs()` slot is stored as a **dense** matrix (`dgeMatrix` or base
`matrix`) — Covid650k (2317) and ATRT (2333) are the canonical
examples. For 650 k × 30 k or 138 k × 18 k data, the dense in-memory
footprint is 100–250 GB, so `readRDS()` alone can blow past the
queue ceiling and trigger `TERM_MEMLIMIT` before the benchmark ever
runs. Symptom in the LSF tail:

```
TERM_MEMLIMIT: job killed after reaching LSF memory usage limit.
Exited with exit code 143.
```

**One-time fix** — convert the source rds to a `dgCMatrix`-backed eset
on a high-mem node, point the YAML at the new file, then submit the
benchmark with normal memory. Everything downstream (`portal_studies.R`,
`prepare_study_from_eset`, the shard writer) works unchanged because
the row-streaming pipeline reads both classes uniformly.

**Step 1 — verify-only probe** (cheap, no write; reports class +
dims + nnz + density):

```sh
./paper/code/portal/sparseify_eset.sh \
    --in /research_jude/.../Studies/Covid650k/expression.rds \
    --verify-only --mem 400000 --wall 1:00
```

If the output reports `class: dgCMatrix`, the rds is already sparse —
skip to Step 3. If it reports `class: dgeMatrix` or `matrix`, proceed
to Step 2.

**Step 2 — submit the conversion job**:

```sh
./paper/code/portal/sparseify_eset.sh \
    --in /research_jude/.../Studies/Covid650k/expression.rds \
    --mem 500000 --wall 8:00 --queue priority
```

This loads the rds once on a high-mem node, runs
`exprs(eset) <- as(m, "CsparseMatrix")`, and writes a new
`expression.sparse.rds` next to the source (overridable with `--out`).
Logs land in `paper/logs/sparseify_<tag>_<timestamp>.{out,err}` and
include load / convert / save wall time plus a size-delta report
(typically 5–10× smaller on disk).

Repeat for the activity matrix if it's also dense:

```sh
./paper/code/portal/sparseify_eset.sh \
    --in /research_jude/.../Studies/Covid650k/activity.rds \
    --mem 200000 --wall 4:00
```

**Step 3 — point the YAML at the converted files**:

```yaml
# paper/code/configs/2317.yaml
input:
    expression: /research_jude/.../Studies/Covid650k/expression.sparse.rds
    activity:   /research_jude/.../Studies/Covid650k/activity.sparse.rds
    networks:   /research_jude/.../Studies/Covid650k/networks.txt
```

**Step 4 — re-run the benchmark with ordinary memory**:

```sh
bash paper/code/portal/portal_studies_single.sh \
    --configs-dir paper/code/configs \
    --queue priority --mem 64000 --cores 1 --wall 4:00 --only 2317
```

Flags accepted by `sparseify_eset.sh`:

| Flag | Default | Meaning |
| --- | --- | --- |
| `--in <path>` | (required) | Source `.rds` (Biobase ExpressionSet). |
| `--out <path>` | `<stem>.sparse.rds` next to source | Destination path. |
| `--force` | — | Overwrite existing destination. |
| `--verify-only` | — | Load + report only; do not write. |
| `--mem <MB>` | `400000` | Memory request. Bump to `500000`+ for 650 k-cell studies. |
| `--cores <n>` | `1` | Cores. The conversion is single-threaded. |
| `--wall <hh:mm>` | `6:00` | Wall-clock limit. |
| `--queue <name>` | `standard` | LSF queue. |
| `--project <id>` | `scminer` | Charge code. |
| `--dry-run` | — | Print the bsub command without submitting. |

### Known limitations

* **Legacy raw-matrix layout (none currently in `configs/`).** YAMLs
  using the historical raw `.rds` matrix + separate `.genes.csv`
  layout are not Biobase `ExpressionSet`s. `portal_studies.R` detects
  `input.genes` in the YAML and emits a `skipped` row (the bsub task
  still exits 0). The 29 active configs no longer include any.

* **R `dgCMatrix` 2^31 nnz cap.** `Matrix::sparseMatrix()` refuses
  more than 2^31 − 1 non-zero entries. `portal_studies.R` reports
  matrix nnz at load time but only warns; the run is allowed to
  proceed since scminerViewer's shard writer (`.write_graph_shards`)
  streams genes one at a time and never materializes a full sparse
  copy. The cap can still bite for *activity* matrices — they pass
  through `.reindex_rows()` which builds a master-shape sparse
  matrix — surfacing as `status = "error-too-large"` in the TSV.

* **Dense-backed source `ExpressionSet`s OOM during `readRDS()`.**
  Studies whose `exprs()` is stored as a `dgeMatrix` / base `matrix`
  rather than `dgCMatrix` can need 150–300 GB to deserialize, even
  when the on-disk rds is only a few GB. Symptom: `TERM_MEMLIMIT`,
  exit 143. Fix: see *Priming large studies* above — one-time
  `sparseify_eset.sh` job, then point the YAML at the new
  `.sparse.rds`.
