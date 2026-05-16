# scMINER Viewer — Applications Note

Short paper modelled after Khatamian *et al.* (*Bioinformatics*
35:2165, 2019). Source manuscript, reproducible benchmarks, and
LSF-driven HPC pipeline for the live scMINER Portal studies.

## Contents

| File | What it is |
| --- | --- |
| [`scminer-viewer.md`](scminer-viewer.md)         | Manuscript draft (~2 pages, Markdown). Pandoc-renderable to PDF / DOCX. |
| [`benchmarks/`](benchmarks/)                     | Synthetic-sweep code (figure 1) + architecture renderer (figure 2). |
| &nbsp;&nbsp;[`methods.R`](benchmarks/methods.R)                | Synthetic study generator + bench primitives. Sourced by `figures.R`. |
| &nbsp;&nbsp;[`figures.R`](benchmarks/figures.R)                | 7x4 synthetic-sweep benchmark (cells × genes) + real 2327 row + discover scaling. Writes 6 standalone figure 1 panels (`figure1_{A..F}_*.{pdf,png}`) plus the scaling TSVs. |
| &nbsp;&nbsp;[`figure_architecture.R`](benchmarks/figure_architecture.R) | Renders the architecture diagram (figure 2). No data deps; ~ 1 s. |
| &nbsp;&nbsp;[`figure_portal.R`](benchmarks/figure_portal.R)             | Renders 6 standalone real-data figures (A–F) from the aggregated portal metrics. Each panel is its own `{pdf,png}` pair so it can be placed independently. |
| &nbsp;&nbsp;[`tables.R`](benchmarks/tables.R)                           | Generates pandoc-renderable tables (table 1 / portal studies, table 2 / figure 1 sweep) from the benchmark TSVs. Writes both `.md` and `.tsv` under `tables/`. |
| &nbsp;&nbsp;[`run_figure1.sh`](benchmarks/run_figure1.sh)               | One-command driver: runs `figures.R` then `tables.R`. Local by default; pass `--hpc --mem ... --wall ...` to submit as a single LSF bsub job. |
| [`configs/`](configs/)                           | 29 YAML configs (one per study, named `<study.ID>.yaml`) pointing at HPC data. |
| [`portal/`](portal/)                             | 29-study HPC pipeline: drivers, R workers, bsub task body, helpers. |
| &nbsp;&nbsp;[`portal_studies.R`](portal/portal_studies.R)               | Loads each YAML, runs `prepare_study_from_eset → load_study → gene_values`, writes per-study metrics. |
| &nbsp;&nbsp;[`portal_studies.sh`](portal/portal_studies.sh)             | Local / interactive wrapper around `portal_studies.R`. |
| &nbsp;&nbsp;[`portal_studies.bsub`](portal/portal_studies.bsub)         | LSF script: single-job mode and job-array task body. |
| &nbsp;&nbsp;[`portal_studies_hpc.sh`](portal/portal_studies_hpc.sh)     | One-command HPC driver — submits the array (one task per study) + dependent merge job. |
| &nbsp;&nbsp;[`portal_studies_single.sh`](portal/portal_studies_single.sh) | One-command HPC driver — submits a single sequential bsub with configurable mem/cores/wall. |
| &nbsp;&nbsp;[`portal_studies_compare.sh`](portal/portal_studies_compare.sh) | One-command HPC driver — runs the with-vs-without-TFsig comparison (29 studies × 2 modes, back-to-back arrays + compare job). |
| &nbsp;&nbsp;[`portal_compare.R`](portal/portal_compare.R)               | Joins the per-mode TSVs into one wide comparison table with `delta_*` columns. |
| &nbsp;&nbsp;[`portal_merge.R`](portal/portal_merge.R)                   | Concatenates per-study TSVs into `metrics/portal_studies.tsv`. |
| &nbsp;&nbsp;[`sparseify_eset.R`](portal/sparseify_eset.R)               | One-time converter: rewrites a dense-backed `ExpressionSet` rds into a `dgCMatrix`-backed one. |
| &nbsp;&nbsp;[`sparseify_eset.sh`](portal/sparseify_eset.sh)             | LSF bsub wrapper around `sparseify_eset.R` with high-mem defaults. |
| &nbsp;&nbsp;[`example.yaml`](portal/example.yaml)                       | Reference YAML schema (annotated template). |
| [`figures/`](figures/)                           | Generated standalone panels: `figure1_{A..F}_*.{pdf,png}` (synthetic, 6 panels), `architecture.{pdf,png}` (figure 2), `figure_portal_{A..F}_*.{pdf,png}` (real-data, 6 panels). |
| [`metrics/`](metrics/)                           | Generated TSVs: `bundle_scaling.tsv` (per-rep), `bundle_scaling_summary.tsv` (mean / sd / SE per config), `discover_scaling.tsv` + `_summary.tsv`, `real_study.tsv`, `portal_studies.tsv` (merged), `portal_studies_<id>.tsv` (per-study), `portal_studies_summary.tsv` (status counts). |
| [`tables/`](tables/)                             | Generated tables (one TSV + one Markdown per table). `table 1 = portal_studies`, `table 2 = figure1_scaling`. Re-render with `Rscript paper/benchmarks/tables.R`. |
| [`logs/`](logs/), [`hpc/`](hpc/)                 | Runtime artefacts (manifests, bsub stdout/stderr). |

## Synthetic-sweep benchmark (figure 1)

The wrapper runs the full pipeline — `figures.R` → `tables.R` — and
tees output to `paper/logs/figure1_<timestamp>.log`:

```sh
# Local (laptop / HPC interactive node)
./paper/benchmarks/run_figure1.sh

# HPC (single bsub job; recommended for the full 28-config sweep
# since 10K x 10K configs may need ~16-32 GB):
./paper/benchmarks/run_figure1.sh --hpc --mem 64000 --wall 6:00

# Inspect the bsub command without submitting:
./paper/benchmarks/run_figure1.sh --hpc --dry-run

# Skip the wrapper and run the R steps directly:
Rscript paper/benchmarks/figures.R
Rscript paper/benchmarks/tables.R
```

The sweep is a 7 × 4 grid: `n_cells ∈ {500, 1000, 2000, 4000, 6000,
8000, 10000}`, `n_genes ∈ {2000, 5000, 8000, 10000}` (28 configs).
Each configuration is benchmarked **`N_REPS = 5`** times with
independent random seeds; mean ± standard error across replicates
drives both the figure 1 error bars and the table 2 cells. Total
wall time: ~ 60–90 min on a laptop (140 runs), ~ 30 s for the
multi-study discover sweep (also 5×6 = 30 runs), < 1 s for the real
2327 row. Drop `N_REPS` to `1` in `figures.R` for a single-run
iteration cycle.

The replicate-level TSVs at `paper/metrics/{bundle,discover}_scaling.tsv`
keep every individual run; the aggregated summaries (one row per
config with `_mean / _sd / _se / n_reps` columns) live at
`paper/metrics/{bundle,discover}_scaling_summary.tsv`.

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

Each panel is rendered as its own `{pdf,png}` pair (no patchwork grid):

| File stem | What it shows |
| --- | --- |
| `figure1_A_size_vs_cells`     | Bundle + shard tree size vs cells, faceted by gene count |
| `figure1_B_load_latency`      | `load_study()` cold-start latency |
| `figure1_C_fetch_latency`     | `gene_values()` median fetch latency (bar to max) |
| `figure1_D_discover_scaling`  | `discover_studies()` scaling on a multi-study root |
| `figure1_E_prepare_time`      | `prepare_study_data()` wall time |
| `figure1_F_peak_memory`       | `prepare_study_data()` peak resident memory |

Edit `SCALING_GRID` / `DISCOVER_GRID` in `figures.R` to widen or
shrink the sweep. Outputs land in `paper/figures/` and
`paper/metrics/{bundle_scaling,discover_scaling,real_study}.tsv`.

## Portal benchmark (29 real studies on HPC)

The `paper/configs/` folder holds **29 YAML configs**, one per study,
named `<study.ID>.yaml`. Each config tells `portal_studies.R` exactly
where its `expression.rds`, `activity.rds`, and `networks.txt` live
on HPC.

### Config schema (`paper/configs/2327.yaml`)

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
* Path resolution: paths inside `input.*` may be absolute (most
  configs are) or relative — relative paths resolve against
  `input_root` if set, else the YAML's own parent dir.
* `output:` is the directory where each study's bundle, gene-shard
  tree, and graph files land. `portal_studies.R` resolves this in
  this order (highest wins):
    1. `--output-root <dir>` CLI flag (single root for every study)
    2. The YAML's `output:` value (per-study, default on HPC)
    3. `--scratch <dir>` (default: an ephemeral `tempfile()`)

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
./paper/portal/portal_studies_hpc.sh --configs-dir paper/configs
```

The driver enumerates every `*.yaml` in `--configs-dir`, writes the
manifest to `paper/hpc/manifest_<timestamp>.txt`, submits a parallel
LSF job array (one task per YAML), and a dependent merge job that
concatenates per-study TSVs into `paper/metrics/portal_studies.tsv`
when the array finishes.

Driver flags:

| Flag | Default | Meaning |
| --- | --- | --- |
| `--configs-dir <dir>` | — | Folder of YAML configs. |
| `--studies-root <dir>` | — | Alternative: folder of per-study subfolders. |
| `--queue "<q1 q2 ...>"` | `"standard priority"` | Space-separated LSF queue list. LSF dispatches each task to the first queue with capacity, so a list of 2–3 queues reduces start-up latency when one is congested. Quote when listing >1. |
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
ls   paper/metrics/portal_studies_*.tsv    # per-study TSVs (one per array task)
cat  paper/metrics/portal_studies.tsv      # final merged table
```

**Aggregate + render real-data figures** once the per-study
TSVs are on disk:

```sh
# 1. Concatenate all per-study TSVs into paper/metrics/portal_studies.tsv:
Rscript paper/portal/portal_merge.R paper/metrics/portal_studies.tsv \
    paper/metrics/portal_studies_*.tsv

# 2. Render 6 standalone panels (one {pdf,png} each)
#    + the status-count summary TSV:
Rscript paper/benchmarks/figure_portal.R
```

Each panel is its own file so it can be placed independently in the
manuscript or supplementary slides — no patchwork grid:

| File stem | What it shows |
| --- | --- |
| `figure_portal_A_size_vs_cells` | Bundle + shard tree size vs cell count (log-log) |
| `figure_portal_B_prepare_time`  | `prepare_study()` wall time vs `n_cells × n_genes` |
| `figure_portal_C_peak_memory`   | Peak memory vs total input size |
| `figure_portal_D_load_latency`  | `load_study()` cold-start latency vs bundle size |
| `figure_portal_E_fetch_latency` | `gene_values()` median fetch latency vs gene count |
| `figure_portal_F_size_ratio`    | Output : input compression ratio per study (split by activity presence) |

If the merged TSV is absent, `figure_portal.R` falls back to globbing
the per-study TSVs directly. The status-count summary lands at
`paper/metrics/portal_studies_summary.tsv` so the manuscript can cite
how many runs landed in each status bucket.

### With-vs-without-TFsig comparison (29 studies × 2 modes)

[`portal_studies_compare.sh`](portal/portal_studies_compare.sh) runs
the benchmark twice per study so you can quantify what TF/sig data
costs in `prepare_study()` time, peak memory, cold/warm `load_study()`
latency, and on-disk bundle size:

* **`expression-only`** — `activity.rds` and `networks.txt` are skipped;
  `prepare_study_from_eset()` receives `activity_eset = NULL` and
  `networks_path = NULL`.
* **`full`** — historical default: expression + activity + networks.

Modes write into separate per-study subdirs
(`<cfg$output>/<studyID>/<mode>/`) so the bundles don't clobber. The
driver submits two 29-task arrays back-to-back, chained on
`ended(expression-only)` so the page cache and the parallel filesystem
aren't shared between modes (which would smear the warm-load timing),
then a small compare job that joins the two TSVs:

```sh
# Dry-run: print the bsub commands without submitting
./paper/portal/portal_studies_compare.sh --configs-dir paper/configs --dry-run

# Submit (29 expression-only -> 29 full -> 1 compare)
./paper/portal/portal_studies_compare.sh --configs-dir paper/configs \
    --mem 64000 --wall 8:00 --queue "standard priority"

# Re-run a subset (e.g. after fixing one study's input)
./paper/portal/portal_studies_compare.sh --configs-dir paper/configs \
    --only 2317,2327
```

Driver flags (same shape as `portal_studies_hpc.sh`):

| Flag | Default | Meaning |
| --- | --- | --- |
| `--configs-dir <dir>` | `paper/configs` | Folder of YAML configs. |
| `--queue "<q1 q2 ...>"` | `"standard priority"` | LSF queue list; tasks per mode are distributed across queues (one sub-array per queue, disjoint index ranges). |
| `--mem <MB>` | `64000` | Memory per task (higher default than the standard run because some studies need it for `full` mode). |
| `--cores <n>` | `4` | Cores per task. |
| `--wall <hh:mm>` | `8:00` | Wall-clock limit per task. |
| `--project <id>` | `scminer` | LSF charge code. |
| `--only <csv>` | (all) | Comma-separated study IDs to limit to. |
| `--dry-run` | — | Print the bsub commands without submitting. |

Per-mode TSVs land at:

```
paper/metrics/portal_studies_<studyID>_expression-only.tsv
paper/metrics/portal_studies_<studyID>_full.tsv
```

The dependent compare job produces the joined wide table:

```sh
paper/metrics/portal_studies_compare.tsv
```

One row per study, with the following metric triples
(`*_full`, `*_expr_only`, `delta_*` where the delta is `full − expr_only`,
so positive = TF/sig adds cost):

| Metric | Columns |
| --- | --- |
| Process wall time | `prepare_seconds_full`, `prepare_seconds_expr_only`, `delta_prepare_seconds` |
| Process peak memory (MB) | `prepare_peak_mb_full`, `prepare_peak_mb_expr_only`, `delta_prepare_peak_mb` |
| Cold load | `load_seconds_full`, `load_seconds_expr_only`, `delta_load_seconds` |
| Warm load | `load_seconds_warm_full`, `load_seconds_warm_expr_only`, `delta_load_seconds_warm` |
| Bundle size | `bundle_bytes_full`, `bundle_bytes_expr_only`, `delta_bundle_bytes` |
| Total output | `total_output_bytes_full`, `total_output_bytes_expr_only`, `delta_total_output_bytes` |

Identity / sanity columns: `studyID`, `n_cells_full`, `n_genes_full`,
`n_clusters_full`, `net_tf_edges_full`, `net_sig_edges_full`,
`status_full`, `status_expr_only`, `note_full`, `note_expr_only`.
If `net_tf_edges_full` / `net_sig_edges_full` are `0` for a study, the
`full` run silently dropped its network — the delta is meaningless for
that row and the YAML's `input.networks` is likely missing or
unreadable.

Re-run just the aggregator (e.g. after re-submitting one mode) without
touching LSF:

```sh
Rscript paper/portal/portal_compare.R \
    --expr-only-glob "paper/metrics/portal_studies_*_expression-only.tsv" \
    --full-glob      "paper/metrics/portal_studies_*_full.tsv" \
    --out            "paper/metrics/portal_studies_compare.tsv"
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
  failure in one mode does **not** block the other 28. Failed per-mode
  tasks leave no TSV; the compare aggregator handles missing rows with
  `NA` deltas.
* For 2317 / Covid650k, prime the dense `expression.rds` with
  [`sparseify_eset.sh`](portal/sparseify_eset.sh) and update its YAML
  before submitting the compare run — see [Priming large studies](#priming-large-studies-dense-backed-expressionsets)
  below.

### Tables

`paper/benchmarks/tables.R` regenerates two manuscript tables from
the same TSVs:

```sh
Rscript paper/benchmarks/tables.R
```

Outputs land under `paper/tables/`:

| Table | Source | Output |
| --- | --- | --- |
| Table 1 — per-study metrics for the 26+ real portal studies | `paper/metrics/portal_studies.tsv` (status == `ok`) | `paper/tables/portal_studies.{md,tsv}` |
| Table 2 — synthetic-sweep grid underlying figure 1 (28 configs once `figures.R` finishes) | `paper/metrics/bundle_scaling.tsv` | `paper/tables/figure1_scaling.{md,tsv}` |

The Markdown files ship a leading caption block and a GFM table that
pandoc converts to LaTeX (or DOCX) without further intervention.
Embed them in `scminer-viewer.md` with a pandoc include filter, or
paste inline. The TSV versions are clean for any downstream
post-processing.

### Local / sequential modes (development & reruns)

```sh
# Local interactive — one or many studies on the current machine:
./paper/portal/portal_studies.sh --configs-dir paper/configs
./paper/portal/portal_studies.sh --config       paper/configs/2327.yaml
./paper/portal/portal_studies.sh --configs-dir  paper/configs --only 2327,2326

# Single bsub (sequential walk, no array) -- wrapper with mem / wall / cores knobs:
./paper/portal/portal_studies_single.sh --configs-dir paper/configs
./paper/portal/portal_studies_single.sh --configs-dir paper/configs --mem 64000 --wall 12:00
./paper/portal/portal_studies_single.sh --configs-dir paper/configs --only 2327,2326

# Single bsub (minimal, env-only -- defaults to mem=32000):
CONFIGS_DIR=$(pwd)/paper/configs bsub < paper/portal/portal_studies.bsub
ONLY=2327,2326 CONFIGS_DIR=$(pwd)/paper/configs bsub < paper/portal/portal_studies.bsub

# Single bsub with mem override (bsub flag overrides #BSUB directive):
CONFIGS_DIR=$(pwd)/paper/configs bsub -R "rusage[mem=64000]" -M 64000 -W 12:00 \
    < paper/portal/portal_studies.bsub
```

`portal_studies.R` flags:

| Flag | Default | Meaning |
| --- | --- | --- |
| `--configs-dir <dir>` | — | Process every YAML in the directory. |
| `--config <yaml>`     | — | Process exactly one YAML (used by array tasks). |
| `--studies-root <dir>` | — | Walk per-study subfolders (legacy mode). |
| `--only a,b,c`        | (all) | CSV of study IDs to restrict to. |
| `--mode <name>`       | `full` | `full` reads expression + activity + networks (historical default). `expression-only` skips activity and networks; used by [`portal_studies_compare.sh`](portal/portal_studies_compare.sh) to isolate the cost of TF/sig data. Each mode writes into its own per-study subdir (`<out>/<studyID>/<mode>/`) so the bundles don't clobber. |
| `--out <tsv>`         | `paper/metrics/portal_studies.tsv` | Output TSV path. |
| `--output-root <dir>` | — | Force every study's bundle / shards / graph files into `<dir>/<studyID>/`, overriding the YAML's `output:` value. |
| `--scratch <dir>`     | tempfile | Fallback when neither `--output-root` nor the YAML's `output:` is set. |
| `--quiet`             | — | Suppress progress messages. |

### Metrics columns (`paper/metrics/portal_studies.tsv`)

| Group | Columns |
| --- | --- |
| Identification | `studyID`, `mode`, `n_cells`, `n_genes`, `n_clusters`, `out_dir` |
| **Inputs** | `expr_input_bytes`, `act_input_bytes`, `net_input_bytes`, `total_input_bytes` |
| **Outputs** | `bundle_bytes`, `shard_bytes`, `graph_bytes`, `total_output_bytes` |
| Wall time & memory | `prepare_seconds`, `prepare_peak_mb`, `load_seconds`, `load_seconds_warm` (re-call of `load_study()` in the same R session, isolates parse cost from cold-page-cache cost) |
| Gene fetch | `fetch_median`, `fetch_mean`, `fetch_max`, `n_fetched` |
| Networks | `net_tf_edges`, `net_sig_edges` |
| Status | `status` (`ok` / `skipped` / `error-too-large` / `error`), `note` |

Status values:

| `status` | When | TSV columns |
| --- | --- | --- |
| `ok` | Benchmark completed | All metric columns filled |
| `skipped` | YAML uses the legacy `input.genes` (raw matrix + `genes.csv`) layout that `prepare_study_from_eset()` can't consume | Metrics NA; `note` carries reason |
| `error-too-large` | Run-time failure attributable to R's `dgCMatrix` 2^31-1 nnz cap (typically inside `.reindex_rows()` for very large activity matrices) | Metrics NA; `note` carries the offending message |
| `error` | Any other failure | Metrics NA; `note` carries the error message |

The bsub task always exits 0 regardless of status (skip / cap / error are all "expected outcomes"), so the merge job's `ended(array)` dependency always fires.

### Resource sizing

Per-study runtime scales roughly linearly with cell count:

* < 10 k cells (Tex, Tregs, ground-truth set): ~ 2–20 min, ≤ 4 GB peak
* 10 k–100 k cells (PBMC, Hepatoblastoma, Glutamatergic): ~ 20–60 min,
  4–12 GB peak
* 100 k+ cells (Covid97k, GSE155446, HMC76k): ~ 1–3 h, 12–24 GB peak
* 650 k cells (Covid650k): bump to `--mem 64000 --wall 12:00`

A full 29-task array completes in roughly the wall time of the
largest single study (the array runs in parallel across hosts).

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
./paper/portal/sparseify_eset.sh \
    --in /research_jude/.../Studies/Covid650k/expression.rds \
    --verify-only --mem 400000 --wall 1:00
```

If the output reports `class: dgCMatrix`, the rds is already sparse —
skip to Step 3. If it reports `class: dgeMatrix` or `matrix`, proceed
to Step 2.

**Step 2 — submit the conversion job**:

```sh
./paper/portal/sparseify_eset.sh \
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
./paper/portal/sparseify_eset.sh \
    --in /research_jude/.../Studies/Covid650k/activity.rds \
    --mem 200000 --wall 4:00
```

**Step 3 — point the YAML at the converted files**:

```yaml
# paper/configs/2317.yaml
input:
    expression: /research_jude/.../Studies/Covid650k/expression.sparse.rds
    activity:   /research_jude/.../Studies/Covid650k/activity.sparse.rds
    networks:   /research_jude/.../Studies/Covid650k/networks.txt
```

**Step 4 — re-run the benchmark with ordinary memory** (the sparse
backing fits comfortably in 64 GB):

```sh
bash paper/portal/portal_studies_single.sh --configs-dir paper/configs \
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

* **Legacy raw-matrix layout (none currently in `configs/`).**
  YAMLs using the historical raw `.rds` matrix + separate `.genes.csv`
  layout (older mammary-gland exports, etc.) are not Biobase
  `ExpressionSet`s. `portal_studies.R` detects `input.genes` in the
  YAML and emits a `skipped` row (the bsub task still exits 0). The
  29 active configs no longer include any; convert to `ExpressionSet`
  on HPC before adding such YAMLs back.

* **R `dgCMatrix` 2^31 nnz cap.** `Matrix::sparseMatrix()` refuses
  more than 2^31-1 non-zero entries. `portal_studies.R` reports
  matrix nnz at load time but only warns; the run is allowed to
  proceed since scminerViewer's shard writer (`.write_graph_shards`)
  streams genes one at a time and never materializes a full sparse
  copy. The cap can still bite for *activity* matrices — they pass
  through `.reindex_rows()` which builds a master-shape sparse
  matrix — and such failures surface as `status = "error-too-large"`
  in the TSV. Expression matrices (e.g. ATRT at 18 k × 138 k,
  2.5 G nnz) bundle normally regardless of nnz.

* **Dense-backed source `ExpressionSet`s OOM during `readRDS()`.**
  Studies whose `exprs()` is stored as a `dgeMatrix` / base `matrix`
  rather than `dgCMatrix` can need 150–300 GB to deserialize, even
  when the on-disk rds is only a few GB. Symptom: `TERM_MEMLIMIT`,
  exit 143. Fix: see *Priming large studies (dense-backed
  ExpressionSets)* above — one-time `sparseify_eset.sh` job, then
  point the YAML at the new `.sparse.rds`.

## Figure (anatomy)

`figure1.pdf` is a 3 × 2 panel:

* **A** — Bundle size vs shard-tree size as a function of cell count,
  faceted by gene count.
* **B** — `load_study()` cold-start latency. Sub-second across the
  sweep because the bundle carries no matrix values.
* **C** — `gene_values()` first-fetch latency (median, error bar to
  max).
* **D** — `discover_studies()` scaling on a multi-study root.
* **E** — `prepare_study_data()` wall time vs cells.
* **F** — `prepare_study_data()` peak resident memory vs cells.

## Manuscript

Plain Markdown so it can be reviewed on GitHub directly or converted
to PDF / DOCX with pandoc:

```sh
# PDF -- xelatex (Unicode-safe). Page geometry / fonts come from the
# YAML frontmatter at the top of scminer-viewer.md, so no -V flags
# needed here. Edit those keys (paperwidth, left/right margins, etc.)
# to change the layout.
pandoc paper/scminer-viewer.md \
       --resource-path=paper \
       --pdf-engine=xelatex \
       -o paper/scminer-viewer.pdf

# DOCX (Unicode-safe; YAML geometry ignored, Word controls layout)
pandoc paper/scminer-viewer.md \
       --resource-path=paper \
       -o paper/scminer-viewer.docx
```

If `xelatex` is missing, install TinyTeX (no sudo, ~150 MB):

```sh
Rscript -e 'install.packages("tinytex"); tinytex::install_tinytex()'
which xelatex      # should resolve to ~/Library/TinyTeX/bin/.../xelatex on macOS
```

Bibliography is inline (no `.bib`) — short Applications-Note format
matches the SJARACNe template (Khatamian *et al.*, 2019).
