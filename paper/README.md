# scMINER Viewer — Applications Note

Short paper modelled after Khatamian *et al.* (*Bioinformatics*
35:2165, 2019). Source manuscript, reproducible benchmarks, and
LSF-driven HPC pipeline for the live scMINER Portal studies.

## Contents

| File | What it is |
| --- | --- |
| [`scminer-viewer.md`](scminer-viewer.md)         | Manuscript draft (~2 pages, Markdown). |
| [`methods.R`](methods.R)                         | Synthetic study generator + bench primitives. Sourced by `figures.R`. |
| [`figures.R`](figures.R)                         | Synthetic-sweep benchmark + real 2327 row; writes `figure1.{pdf,png}` and the scaling TSVs. |
| [`configs/`](configs/)                           | 34 YAML configs (one per study, named `<study.ID>.yaml`) pointing at HPC data. |
| [`portal_studies.R`](portal_studies.R)           | Loads each YAML, runs `prepare_study_from_eset → load_study → gene_values`, writes per-study metrics. |
| [`portal_studies.sh`](portal_studies.sh)         | Local / interactive wrapper around `portal_studies.R`. |
| [`portal_studies.bsub`](portal_studies.bsub)     | LSF script: single-job mode and job-array task body. |
| [`portal_studies_hpc.sh`](portal_studies_hpc.sh) | One-command HPC driver — submits the array (one task per study) + dependent merge job. |
| [`portal_studies_single.sh`](portal_studies_single.sh) | One-command HPC driver — submits a single sequential bsub with configurable mem/cores/wall. |
| [`portal_merge.R`](portal_merge.R)               | Concatenates per-study TSVs into `metrics/portal_studies.tsv`. |
| [`figures/`](figures/)                           | Generated `figure1.{pdf,png}`. |
| [`metrics/`](metrics/)                           | Generated TSVs: `bundle_scaling.tsv`, `discover_scaling.tsv`, `real_study.tsv`, `portal_studies.tsv`. |
| [`logs/`](logs/), [`hpc/`](hpc/)                 | Runtime artefacts (manifests, bsub stdout/stderr). |

## Synthetic-sweep benchmark (figure 1)

From the project root with `scminerViewer` installed:

```sh
Rscript paper/figures.R
```

Run time ≈ 1–2 min on a laptop; ≈ 5 s for the multi-study
discover sweep; < 1 s for the real 2327 row. Outputs land in
`paper/figures/` and `paper/metrics/`. Edit `SCALING_GRID` /
`DISCOVER_GRID` in `figures.R` to widen or shrink the sweep.

## Portal benchmark (34 real studies on HPC)

The `paper/configs/` folder holds **34 YAML configs**, one per study,
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
./paper/portal_studies_hpc.sh --configs-dir paper/configs
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

### Local / sequential modes (development & reruns)

```sh
# Local interactive — one or many studies on the current machine:
./paper/portal_studies.sh --configs-dir paper/configs
./paper/portal_studies.sh --config       paper/configs/2327.yaml
./paper/portal_studies.sh --configs-dir  paper/configs --only 2327,2326

# Single bsub (sequential walk, no array) -- wrapper with mem / wall / cores knobs:
./paper/portal_studies_single.sh --configs-dir paper/configs
./paper/portal_studies_single.sh --configs-dir paper/configs --mem 64000 --wall 12:00
./paper/portal_studies_single.sh --configs-dir paper/configs --only 2327,2326

# Single bsub (minimal, env-only -- defaults to mem=32000):
CONFIGS_DIR=$(pwd)/paper/configs bsub < paper/portal_studies.bsub
ONLY=2327,2326 CONFIGS_DIR=$(pwd)/paper/configs bsub < paper/portal_studies.bsub

# Single bsub with mem override (bsub flag overrides #BSUB directive):
CONFIGS_DIR=$(pwd)/paper/configs bsub -R "rusage[mem=64000]" -M 64000 -W 12:00 \
    < paper/portal_studies.bsub
```

`portal_studies.R` flags:

| Flag | Default | Meaning |
| --- | --- | --- |
| `--configs-dir <dir>` | — | Process every YAML in the directory. |
| `--config <yaml>`     | — | Process exactly one YAML (used by array tasks). |
| `--studies-root <dir>` | — | Walk per-study subfolders (legacy mode). |
| `--only a,b,c`        | (all) | CSV of study IDs to restrict to. |
| `--out <tsv>`         | `paper/metrics/portal_studies.tsv` | Output TSV path. |
| `--output-root <dir>` | — | Force every study's bundle / shards / graph files into `<dir>/<studyID>/`, overriding the YAML's `output:` value. |
| `--scratch <dir>`     | tempfile | Fallback when neither `--output-root` nor the YAML's `output:` is set. |
| `--quiet`             | — | Suppress progress messages. |

### Metrics columns (`paper/metrics/portal_studies.tsv`)

| Group | Columns |
| --- | --- |
| Identification | `studyID`, `n_cells`, `n_genes`, `n_clusters` |
| **Inputs** | `expr_input_bytes`, `act_input_bytes`, `net_input_bytes`, `total_input_bytes` |
| **Outputs** | `bundle_bytes`, `shard_bytes`, `graph_bytes`, `total_output_bytes` |
| Wall time & memory | `prepare_seconds`, `prepare_peak_mb`, `load_seconds` |
| Gene fetch | `fetch_median`, `fetch_mean`, `fetch_max`, `n_fetched` |
| Networks | `net_tf_edges`, `net_sig_edges` |

### Resource sizing

Per-study runtime scales roughly linearly with cell count:

* < 10 k cells (Tex, Tregs, ground-truth set): ~ 2–20 min, ≤ 4 GB peak
* 10 k–100 k cells (PBMC, Hepatoblastoma, Glutamatergic): ~ 20–60 min,
  4–12 GB peak
* 100 k+ cells (Covid97k, GSE155446, HMC76k): ~ 1–3 h, 12–24 GB peak
* 650 k cells (Covid650k): bump to `--mem 64000 --wall 12:00`

A full 34-task array completes in roughly the wall time of the
largest single study (the array runs in parallel across hosts).

### Known limitations

* `2343.yaml` – `2347.yaml` (KKYan mammary-gland nest substudies) use
  a legacy raw-matrix + `genes.csv` layout and are not `ExpressionSet`s.
  `portal_studies.R` stops with a clear error on those configs; the
  driver array task will fail fast. Convert to `ExpressionSet` on
  HPC before re-running.

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
pandoc paper/scminer-viewer.md \
       --resource-path=paper \
       -o paper/scminer-viewer.pdf
```

Bibliography is inline (no `.bib`) — short Applications-Note format
matches the SJARACNe template (Khatamian *et al.*, 2019).
