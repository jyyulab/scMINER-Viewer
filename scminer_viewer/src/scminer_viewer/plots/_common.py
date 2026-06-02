"""Helpers shared across plot modules."""

from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import pandas as pd

from ..data import Study

DEFAULT_COLORS = [
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
]


def cluster_color_map(study: Study) -> dict[str, str]:
    colors = list(study.clusters["color"])
    if any(c is None or (isinstance(c, str) and not c) for c in colors):
        colors = [DEFAULT_COLORS[i % len(DEFAULT_COLORS)]
                  for i in range(len(colors))]
    return dict(zip(study.clusters.index, colors))


def cells_mask(
    study: Study,
    active_clusters: Optional[Iterable[str]],
    cell_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Intersect the cluster filter with an optional per-cell mask.

    `cell_mask` (boolean ndarray of length n_cells) is typically supplied
    by the downsampling reactive.
    """
    if active_clusters is None:
        base = np.ones(study.n_cells, dtype=bool)
    else:
        active = set(active_clusters)
        base = study.cells["cellType"].isin(active).to_numpy()
    if cell_mask is not None:
        base = base & np.asarray(cell_mask, dtype=bool)
    return base


def aggregate_by_cluster(
    study: Study,
    genes: list[str],
    relationship: str,
    active_clusters: Optional[list[str]] = None,
    cell_mask: Optional[np.ndarray] = None,
) -> Optional[tuple[pd.DataFrame, pd.DataFrame]]:
    """Return (mean, pct_expressing) per gene per cluster (lazy reads).

    Each gene's row is loaded from disk via `study.gene_values(...)`.
    Optional `cell_mask` restricts the per-cluster aggregation to the
    sampled subset of cells (used by the Sampling % control).
    """
    if not genes:
        return None
    present_genes: list[str] = []
    rows: list[np.ndarray] = []
    for g in genes:
        vals = study.gene_values(g, relationship)
        if vals is None:
            continue
        present_genes.append(g)
        rows.append(vals)
    if not present_genes:
        return None

    clusters = active_clusters or list(study.clusters.index)
    cell_types = study.cells["cellType"].to_numpy()
    sample = (np.asarray(cell_mask, dtype=bool)
              if cell_mask is not None
              else np.ones(len(cell_types), dtype=bool))
    means = np.zeros((len(present_genes), len(clusters)), dtype=np.float64)
    pcts = np.zeros_like(means)
    for i, vals in enumerate(rows):
        for j, cluster in enumerate(clusters):
            mask = (cell_types == cluster) & sample
            sub = vals[mask]
            sub = sub[~np.isnan(sub)]
            n = sub.size
            if n == 0:
                continue
            means[i, j] = float(sub.mean())
            pcts[i, j] = float((sub != 0).sum()) / n
    return (
        pd.DataFrame(means, index=present_genes, columns=clusters),
        pd.DataFrame(pcts, index=present_genes, columns=clusters),
    )


# Injected into the page once to enforce equal-aspect ratio on every plotly
# scatter plot after box zoom or initial render. Uses MutationObserver to
# attach to each .js-plotly-plot element as plotly creates it. Mirrors the
# htmlwidgets::onRender approach used in the R package.
ASPECT_CORRECT_JS = r"""
(function () {
  var attached = new WeakSet();

  // Only correct scatter/scattergl plots whose both axes are continuous.
  // Heatmap, violin, and bubble plots are excluded:
  //   - heatmap/violin: non-scatter trace types -> allScatter false
  //   - bubble: scatter traces but x-axis is categorical (cluster names)
  function isEmbeddingPlot(el) {
    if (!el._fullData || !el._fullData.length) return false;
    // Opt-out: plots that set layout.meta.scvFill stretch to fill their
    // panel (the cluster overview) instead of locking to equal aspect.
    var meta = el._fullLayout && el._fullLayout.meta;
    if (meta && meta.scvFill) return false;
    var allScatter = el._fullData.every(function (t) {
      return t.type === 'scattergl' || t.type === 'scatter';
    });
    if (!allScatter) return false;
    var xa = el._fullLayout && el._fullLayout.xaxis;
    var ya = el._fullLayout && el._fullLayout.yaxis;
    if (!xa || !ya) return false;
    return xa.type !== 'category' && ya.type !== 'category';
  }

  function correctAspect(el) {
    if (!isEmbeddingPlot(el)) return null;
    var xa = el._fullLayout.xaxis;
    var ya = el._fullLayout.yaxis;
    if (!xa._length || !ya._length) return null;
    var xUpp = Math.abs(xa.range[1] - xa.range[0]) / xa._length;
    var yUpp = Math.abs(ya.range[1] - ya.range[0]) / ya._length;
    var tol  = 1e-6 * Math.max(xUpp, yUpp);
    if (Math.abs(xUpp - yUpp) <= tol) return null;
    var upd = {};
    if (xUpp > yUpp) {
      var halfY = xUpp * ya._length / 2;
      var midY  = (ya.range[0] + ya.range[1]) / 2;
      upd['yaxis.range[0]'] = midY - halfY;
      upd['yaxis.range[1]'] = midY + halfY;
    } else {
      var halfX = yUpp * xa._length / 2;
      var midX  = (xa.range[0] + xa.range[1]) / 2;
      upd['xaxis.range[0]'] = midX - halfX;
      upd['xaxis.range[1]'] = midX + halfX;
    }
    return upd;
  }

  function attachTo(el) {
    if (attached.has(el)) return;
    attached.add(el);
    var busy = false;
    function run() {
      if (busy) return;
      var upd = correctAspect(el);
      if (!upd) return;
      busy = true;
      Plotly.relayout(el, upd).then(function () { busy = false; });
    }
    el.on('plotly_afterplot', run);
    el.on('plotly_relayout', function (ed) {
      if (!Object.keys(ed).some(function (k) {
        return /\.range\[/.test(k);
      })) return;
      run();
    });
  }

  function scan() {
    document.querySelectorAll('.js-plotly-plot').forEach(attachTo);
  }

  new MutationObserver(scan).observe(document.documentElement, {
    childList: true, subtree: true,
    attributes: true, attributeFilter: ['class'],
  });
  scan();
})();
"""


def empty_figure(title: str = ""):
    import plotly.graph_objects as go
    fig = go.Figure()
    if title:
        fig.update_layout(title=title)
    fig.update_layout(margin=dict(l=20, r=20, t=40, b=20))
    return fig
