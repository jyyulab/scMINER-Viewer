# scMINER Viewer — Applications Note

Short paper modelled after Khatamian *et al.* (*Bioinformatics*
35:2165, 2019). Source manuscript and reproducible benchmarks.

## Contents

| File | What it is |
| --- | --- |
| [`scminer-viewer.md`](scminer-viewer.md)      | Manuscript draft (~2 pages, Markdown). |
| [`methods.R`](methods.R)                      | Benchmarking helpers — synthetic study generator, timing primitives, bundle / load / fetch / discover benchmarks. Sourced by `figures.R`. |
| [`figures.R`](figures.R)                      | Runs the full benchmark sweep + the real 2327 study, writes a multi-panel PDF/PNG figure and metrics TSVs. |
| [`figures/`](figures/)                        | Generated figure files: `figure1.pdf` + `figure1.png`. |
| [`metrics/`](metrics/)                        | Generated TSVs: `bundle_scaling.tsv`, `discover_scaling.tsv`, `real_study.tsv`. |

## Reproduce

From the project root, with `scminerViewer` installed (and the real
2327 study optionally present at `data/2327/2327.scminer.h5` for the
real-study row):

```sh
Rscript paper/figures.R
```

Run time on commodity laptop hardware: ≈ 1–2 min for the synthetic
sweep, ~ 5 s for the multi-study discover sweep, < 1 s for the real
2327 row. Outputs land in `paper/figures/` and `paper/metrics/`.

To customise the sweep, edit `SCALING_GRID` and `DISCOVER_GRID` in
`figures.R`.

## Figure (anatomy)

`figure1.pdf` is a 2 × 2 panel:

* **A** — Bundle size vs shard-tree size as a function of cell count,
  faceted by gene count. Lazy-mode bundles stay flat in MB while the
  underlying shard tree grows linearly with the matrix size.
* **B** — `load_study()` cold-start latency. Sub-second across the
  sweep because the bundle carries no matrix values.
* **C** — `gene_values()` first-fetch latency (median, error bar to
  max). One shard read per gene; tens of ms at typical sizes.
* **D** — `discover_studies()` scaling on a multi-study root. Linear
  in the number of studies; constant per-study cost ~ 50 ms.

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
