# scMINER Viewer

Standalone, offline viewer for [scMINER](https://github.com/jyyulab/scMINER)
single-cell studies. Two sibling packages share a small HDF5 bundle
format and the on-disk shard tree:

| Package | Lang. | What it does |
| --- | --- | --- |
| [`scminerViewer/`](scminerViewer/)   | R      | Prepares scMINER studies **from a Biobase `ExpressionSet`** (writes the graph-import layout + the `.scminer.h5` bundle) and serves a Shiny app + multi-study browser. |
| [`scminer_viewer/`](scminer_viewer/) | Python | Prepares scMINER studies **from an `AnnData` (.h5ad)** input *and* reads any `.scminer.h5` bundle (R- or Python-written); serves a Shiny-for-Python single-study viewer + multi-study browser. |

Twin writers, twin readers — bundles produced by either side are
byte-compatible and interchangeable. See
[`paper/figures/architecture.png`](paper/figures/architecture.png) for
the full data-flow diagram.

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

## 🚀 Try the demo

Three ways to see the viewer, in increasing order of "real":

| Path | Size | How |
| --- | --- | --- |
| **In-package demo** — 2327 (Tex) bundle + 200 curated gene shards (canonical T-cell / exhaustion markers). Already shipped inside the R package. | ~80 MB (installed with the package) | `scminerViewer::run_demo()` |
| **Pre-processed bundle** — full 2327 bundle + complete shard tree (every gene). Fastest way to browse all ~10k genes without running `prepare_study()`. | ~484 MB download | Download `2327-processed.tar.gz` from the [`demo-data-v1`](https://github.com/hzhou98/scMINER-Viewer/releases/tag/demo-data-v1) release, untar into `data/2327/`, then `run_app("data/2327/2327.scminer.h5")`. |
| **Run the full pipeline** — download the source `expression.rds` (22 MB), `activity.rds` (236 MB), `networks.txt` (57 MB), write a YAML, run `prepare_study()`. The closest thing to "I have my own scMINER data". | ~315 MB download + 1–3 min processing | See [Tutorial § A.1](book/03-tutorial.Rmd) (rendered: `docs/tutorial.html`). |

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
├── book/, docs/                       # bookdown sources + rendered site
├── data/                              # multi-study root (gitignored)
│   ├── input/2327/                    # raw RDS/TSV inputs (from demo-data-v1)
│   └── 2327/                          # per-study processed subfolder
│       ├── 2327.scminer.h5
│       ├── Cell/, Gene/, Network_*/
│       ├── study_meta/2327_study_meta.csv
│       ├── study_gene_{expression,tf,sig}/2327_*.csv
│       ├── expression_files/2327/{meta.csv, <letter>/<gene>.csv.gz}
│       └── activity_files/2327/{meta.csv, TF/<letter>/..., SIG/<letter>/...}
├── scminerViewer/                     # R package
│   ├── R/, tests/, inst/extdata/      # incl. example_config.yml + 200-gene mini demo
│   └── inst/scripts/{build_2327_bundle.R, build_demo_data.R}
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

**Get the 2327 demo inputs** — required for (B) below; skip if you'll
only use (Z), (A), or (D) against the in-package or already-built
bundle. The three source files are published as a GitHub Release
asset (315 MB total):

```sh
# Create the input dir under the (gitignored) data/ tree
mkdir -p data/input/2327

# Download the three source files into it
RELEASE_URL=https://sjcrh-my.sharepoint.com/:f:/r/personal/hzhou98_stjude_org/Documents/DevOps/release-assets?csf=1&web=1&e=uiRyM4
curl -L -o data/input/2327/expression.rds "$RELEASE_URL/expression.rds"   #  22 MB
curl -L -o data/input/2327/activity.rds   "$RELEASE_URL/activity.rds"     # 236 MB
curl -L -o data/input/2327/networks.txt   "$RELEASE_URL/networks.txt"     #  57 MB

# Optional but recommended — verify the 236 MB activity.rds especially
shasum -a 256 data/input/2327/{expression.rds,activity.rds,networks.txt}
# expression.rds  9c28047387a181263107cab3076d426aaa32fcc402e6a2af6b4fb1ec5b910960
# activity.rds    10b94c7826ce513c792ce0cdd57f6c1cad1adb3576f8f6f17779d35b6413c487
# networks.txt    9096a1c3604a93e19e2ac177728686ca0f45c6de48a9b27149fc6e6e9b8dc856
```

Prefer the pre-built bundle? Grab `2327-processed.tar.gz` (~484 MB)
from the same release, untar into `data/`, and skip straight to (D)
— see the table in [🚀 Try the demo](#-try-the-demo) above.

**Use it** — five things you can do (Z is the fastest):

```sh
# (Z) Zero-config demo — open the in-package 2327 mini-bundle (200 genes)
Rscript -e 'scminerViewer::run_demo(port = 8000)'

# (A) Build the 2327 bundle from an on-disk graph layout (no RDS needed)
Rscript scminerViewer/inst/scripts/build_2327_bundle.R
# → writes data/2327/2327.scminer.h5

# (B) Full prepare_study pipeline driven by a YAML config. config-2327.yml
#     ships at the repo root and points at the files you downloaded above
#     under data/input/2327/. Writes data/2327/2327.scminer.h5 + the full
#     shard tree.
Rscript -e 'scminerViewer::prepare_study("config-2327.yml")'

# (C) Multi-study browser — one URL per study, card-grid index of every
#     <sid>/<sid>.scminer.h5 under the root. Shards default to
#     dirname(bundle_path); pass shard_dir = ... only if you've split
#     the bundles and shard trees onto different parents.
Rscript -e 'scminerViewer::run_browser("data", port = 8000)'

# (D) Single-study app
Rscript -e 'scminerViewer::run_app("data/2327/2327.scminer.h5", port = 8000)'
```

Full multi-study walkthrough (download demo data → prep → browse → add
another) is in the bookdown [Tutorial chapter](book/03-tutorial.Rmd)
(rendered: `docs/tutorial.html`). Section A walks through downloading
the 2327 inputs from the [`demo-data-v1`](https://github.com/hzhou98/scMINER-Viewer/releases/tag/demo-data-v1)
release and running `prepare_study()` end-to-end.

The package also ships a fully-annotated YAML template at
`system.file("extdata", "example_config.yml", package = "scminerViewer")`
— copy it as the starting point for any new study.

See [`scminerViewer/README.md`](scminerViewer/README.md) for the YAML
config schema, the full exported API, and how to prepare a fresh study
from a Biobase `ExpressionSet`.

### Python side

```sh
python3 -m venv .venv && source .venv/bin/activate
pip install -e "scminer_viewer[dev]"          # viewer + test deps
# Add `[prepare]` if you also want to build bundles from AnnData:
pip install -e "scminer_viewer[prepare,dev]"

# Read a bundle (either R- or Python-written)
scminer-viewer info   data/2327/2327.scminer.h5
scminer-viewer run    data/2327/2327.scminer.h5 --port 8000
scminer-viewer browse data/                     --port 8000
scminer-viewer list   data/

# Or write one yourself from an AnnData + YAML config
scminer-viewer prepare config-2327.yml
```

See [`scminer_viewer/README.md`](scminer_viewer/README.md) for the
`Study` dataclass + programmatic API (`load_study`, `gene_values`,
`prepare_study`, `run_app`, `run_browser`).

## Exposing to a remote host (HPC / cluster nodes)

By default both `run` and `browse` bind to `127.0.0.1`, so only the
local machine can connect. When you're on an HPC compute node and want
to view from your workstation, you have two options.

**Option 1 — SSH tunnel (recommended).** Leave the viewer on
`127.0.0.1` and forward the port back to your workstation. Nothing on
the cluster's network is opened, and you piggy-back on SSH's auth:

```sh
# On the HPC node:
scminer-viewer run /path/to/2327.scminer.h5 --port 8000 --no-browser

# On your workstation (in another shell):
ssh -L 8000:localhost:8000 user@<hpc-node>
# Then open http://localhost:8000 in your local browser.
```

If your HPC requires logging into a head node first, chain a jump:
`ssh -L 8000:<compute-node>:8000 -J user@head user@<compute-node>`.

**Option 2 — bind to a public interface.** Use `--allow-remote`
(shortcut for `--host 0.0.0.0 --no-browser`) and connect directly:

```sh
# On the HPC node:
scminer-viewer run /path/to/2327.scminer.h5 --port 8000 --allow-remote

# On your workstation:
# open http://<hpc-node-hostname-or-ip>:8000
```

Two caveats: (1) the viewer has **no authentication** — anyone routable
to the node can read the study; (2) most clusters firewall arbitrary
ports inbound, so you may need to coordinate with admins or pick from
an allowed range. Prefer Option 1 unless you specifically need
unauthenticated access from a shared LAN.

`browse` takes the same flags:

```sh
scminer-viewer browse /path/to/data --port 8000 --allow-remote
```

## Bundle format (R ↔ Python contract)

Documented in full in [`IMPLEMENTATION.md`](IMPLEMENTATION.md). Highlights:

- Single HDF5 file, `bundleVersion = 1`.
- Strings are UTF-8.
- Three optional gene-index datasets under `/index/` enumerate which
  genes have shards in each matrix; the apps gate gene selection on
  these.
- An optional `/defaults/genes` dataset lists genes the app should
  auto-load on startup (mirrors the Vue portal's `preGenes`).
- The bundle and the shard tree are normally co-located:
  `load_study(bundle_path)` defaults `shard_dir = dirname(bundle_path)`.

## Status

- **R**: 39 `test_that()` blocks / 202 `expect_*()` assertions across 8
  test files (`scminerViewer/tests/testthat/`). Tests that need
  optional deps (`R.utils`, `ggsci`, etc.) skip gracefully.
- **Python**: 43 tests across 4 test files (`scminer_viewer/tests/`).
  Reader fixtures are built by shelling out to `Rscript`, so the
  Python reader is verified against bytes produced by the R writer
  (and vice versa — `test_prepare.py` round-trips a Python-written
  bundle through the same reader).
- **Bundle contract**: a `.scminer.h5` written by either side is read
  by both. `bundleVersion = 1`. See
  [`IMPLEMENTATION.md`](IMPLEMENTATION.md) for the field-by-field spec.
- **2327 bundle**: 77 MB at `data/2327/2327.scminer.h5`. 8464 cells ×
  9861 genes × 3 clusters; 9861 expression / 925 TF / 4708 SIG genes
  indexed; 743K total network edges; values fetched lazily from the
  shard tree.

## License

scMINER Viewer is released under the **Apache License, Version 2.0**.
See [`LICENSE`](LICENSE) for the full text and [`NOTICE`](NOTICE) for
attribution requirements. Both packages (`scminerViewer/`,
`scminer_viewer/`) carry symlinks back to the root `LICENSE` so each
distribution unit ships a copy.

    Copyright 2026 St. Jude Children's Research Hospital

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

