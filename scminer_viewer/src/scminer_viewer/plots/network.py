from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
import plotly.graph_objects as go

from ..data import Study
from ._common import empty_figure

try:
    import networkx as nx  # noqa: F401 — imported for legacy callers
    _HAS_NX = True
except ImportError:  # pragma: no cover
    _HAS_NX = False


# Direction colors mirror the R viewer (scminerViewer/R/plots.R):
# positive correlation = activator (teal), negative = repressor (red).
_DIR_COLOR = {
    "Activator": "#2E7D6A",
    "Repressor": "#D7493A",
    "Unknown":   "#9E9E9E",
}

# Visual constants tuned for a clean circular layout.
_RING_RADIUS  = 1.00      # neighbors sit on this ring
_LABEL_RADIUS = 1.18      # labels float just outside the rim
_PLOT_PAD     = 1.45      # axis range — enough breathing room for labels
_FOCUS_COLOR  = "#F2A33A"   # warm amber for the centre gene
_NEIGH_COLOR  = "#3F5775"   # cool slate for the periphery
_RING_COLOR   = "rgba(60, 60, 60, 0.10)"  # subtle decorative ring
_BG_COLOR     = "#FBFBF9"   # soft off-white plot background


def _direction_for(row) -> str:
    """Sign of spearman (fall back to pearson) → Activator/Repressor."""
    score = row.get("spearman")
    if score is None or (isinstance(score, float) and np.isnan(score)):
        score = row.get("pearson")
    if score is None or (isinstance(score, float) and np.isnan(score)):
        return "Unknown"
    return "Activator" if float(score) >= 0 else "Repressor"


def _label_position(angle_rad: float) -> str:
    """Octant-based plotly textposition so labels read outward from centre."""
    # Normalize to [0, 2π)
    a = angle_rad % (2.0 * np.pi)
    # Eight octants centred on 0°, 45°, 90°, ...
    octant = int((a + np.pi / 8) // (np.pi / 4)) % 8
    return [
        "middle right",   # 0°   — 3 o'clock
        "top right",      # 45°  — 1:30
        "top center",     # 90°  — 12 o'clock
        "top left",       # 135° — 10:30
        "middle left",    # 180° — 9 o'clock
        "bottom left",    # 225° — 7:30
        "bottom center",  # 270° — 6 o'clock
        "bottom right",   # 315° — 4:30
    ][octant]


def network_plot(
    study: Study,
    gene: str,
    network_type: str = "TF",
    active_clusters: Optional[Iterable[str]] = None,
    max_edges: int = 60,
) -> go.Figure:
    edges_df = (
        study.network_tf if network_type == "TF" else study.network_sig
    )
    if edges_df is None or len(edges_df) == 0:
        return empty_figure(f"No {network_type} network in this study")

    mask = (edges_df["source"] == gene) | (edges_df["target"] == gene)
    if active_clusters is not None:
        mask &= edges_df["cellType"].isin(list(active_clusters))
    sub = edges_df.loc[mask].copy()
    if len(sub) == 0:
        return empty_figure(f"No {network_type} edges for {gene}")

    sub = sub.reindex(sub["mi"].abs().sort_values(ascending=False).index)
    sub = sub.head(max_edges)
    sub["direction"] = sub.apply(_direction_for, axis=1)
    sub["abs_mi"] = sub["mi"].abs()

    # --- Neighbor ordering --------------------------------------------------
    # Group by direction (Activators first, then Repressors, then Unknown),
    # and within each group sort by descending |MI|. The result groups
    # similarly-colored edges into contiguous arcs of the circle, which
    # reads cleanly at a glance and makes activator/repressor balance
    # immediately visible.
    neighbor_info: list[tuple[str, str, float]] = []
    seen: set[str] = set()
    # Sub is already sorted strongest-MI first, so the first occurrence of
    # each neighbor carries that neighbor's dominant direction.
    for _, r in sub.iterrows():
        nb = r["target"] if r["source"] == gene else r["source"]
        if nb == gene or nb in seen:
            continue
        seen.add(nb)
        neighbor_info.append((nb, r["direction"], float(r["abs_mi"])))

    order_key = {"Activator": 0, "Unknown": 1, "Repressor": 2}
    neighbor_info.sort(key=lambda x: (order_key.get(x[1], 3), -x[2]))
    neighbors = [n[0] for n in neighbor_info]

    # Place neighbors evenly around the circle, starting at 12 o'clock and
    # walking clockwise. Activators end up on the right half, repressors
    # on the left half (modulo unbalanced counts).
    pos = {gene: (0.0, 0.0)}
    node_angles: dict[str, float] = {}
    if neighbors:
        step = 2.0 * np.pi / len(neighbors)
        for i, nb in enumerate(neighbors):
            # 12 o'clock = π/2 in standard math convention; clockwise means
            # subtracting angles as i increases.
            theta = np.pi / 2.0 - i * step
            pos[nb] = (_RING_RADIUS * np.cos(theta),
                       _RING_RADIUS * np.sin(theta))
            node_angles[nb] = theta

    # --- Edge widths --------------------------------------------------------
    # Range 1.4 px (smallest edge) -> 3.5 px (largest edge). The min is
    # anchored so weak edges stay visible; the max is intentionally
    # restrained so heavy edges don't dominate the figure visually.
    mi_vals = sub["abs_mi"].to_numpy()
    _WIDTH_MIN, _WIDTH_MAX = 1.4, 3.5
    if mi_vals.size and float(mi_vals.max()) > 0:
        widths = _WIDTH_MIN + (_WIDTH_MAX - _WIDTH_MIN) * (
            mi_vals / float(mi_vals.max())
        )
    else:
        widths = np.full(len(sub), _WIDTH_MIN)

    # --- Decorative backdrop ring ------------------------------------------
    ring_theta = np.linspace(0, 2.0 * np.pi, 180)
    ring_trace = go.Scatter(
        x=_RING_RADIUS * np.cos(ring_theta),
        y=_RING_RADIUS * np.sin(ring_theta),
        mode="lines",
        line=dict(color=_RING_COLOR, width=1, dash="dot"),
        hoverinfo="skip", showlegend=False, name="ring",
    )

    # --- Edges + arrows in one shot, as plotly annotations -----------------
    # Why annotations rather than Scatter lines + triangle markers:
    # plotly's `marker.angle` rotation is finicky across versions and the
    # triangle endpoints can look symmetric at small sizes, leaving the
    # direction visually ambiguous. Annotations draw a real arrow — a
    # line whose end terminates in a crisp pixel-sized arrowhead pointing
    # at the given (x, y). One annotation per edge keeps the per-edge
    # color / width / direction encoding intact.
    #
    # Visual orientation: tail = focus gene (centre); head = the
    # **non-focus** end of the edge. The arrow always points outward,
    # landing at the rim next to the connected gene, regardless of which
    # side is `source` vs `target` in the underlying network. The
    # semantic source → target direction is preserved in the hover
    # tooltip (`<b>source → target</b>`) so it's still readable.
    edge_annotations = []
    edge_traces = []           # invisible zero-width carriers — keep the
                                # one-trace-per-edge contract for tests &
                                # downstream tooling that introspects fig.data
    fx, fy = pos[gene]
    for (_, r), w in zip(sub.iterrows(), widths):
        neighbor = r["target"] if r["source"] == gene else r["source"]
        nx_, ny_ = pos[neighbor]
        color = _DIR_COLOR[r["direction"]]
        edge_annotations.append(dict(
            x=nx_, y=ny_,         # arrowhead lands at neighbor (rim)
            ax=fx, ay=fy,         # arrow tail starts at focus (centre)
            xref="x", yref="y", axref="x", ayref="y",
            showarrow=True,
            arrowhead=2,          # filled triangle — clearest direction cue
            arrowsize=1.0,        # head size multiplier (pixels) — small + clean
            arrowwidth=float(w),  # shaft width — encodes |MI|
            arrowcolor=color,
            opacity=0.85,
            standoff=14,          # px offset from (x, y); keeps the head
                                  # outside the neighbor node circle
            startstandoff=18,     # px offset from (ax, ay); keeps the
                                  # tail outside the focus node circle
        ))
        # Zero-opacity line carrying the same (color, name) so legacy
        # tests / tooling can still pair lines with arrowheads by index.
        edge_traces.append(go.Scatter(
            x=[fx, nx_], y=[fy, ny_], mode="lines",
            line=dict(width=0.001, color=color),
            opacity=0.0,
            hoverinfo="skip", showlegend=False,
            name=f"edge-{r['source']}-{r['target']}",
        ))

    # --- Hover overlay -----------------------------------------------------
    hover_x, hover_y, hover_texts = [], [], []
    for _, r in sub.iterrows():
        x0, y0 = pos[r["source"]]
        x1, y1 = pos[r["target"]]
        hover_x.append((x0 + x1) / 2.0)
        hover_y.append((y0 + y1) / 2.0)
        hover_texts.append(
            f"<b>{r['source']} → {r['target']}</b><br>"
            f"Direction: <b>{r['direction']}</b><br>"
            f"MI: {r['mi']:.3f}<br>"
            f"Spearman: {_fmt(r.get('spearman'))}<br>"
            f"Pearson: {_fmt(r.get('pearson'))}<br>"
            f"p-value: {_fmt(r.get('pvalue'), spec='.2g')}<br>"
            f"Cell type: {r['cellType']}"
        )
    hover_trace = go.Scatter(
        x=hover_x, y=hover_y, mode="markers",
        marker=dict(size=10, color="rgba(0,0,0,0)"),
        hovertext=hover_texts, hoverinfo="text",
        showlegend=False, name="edge",
    )

    # The old triangle-marker arrow trace is gone — annotations above
    # carry the real arrowheads now. Keep a thin sentinel `Scatter` with
    # the per-edge neighbor positions + colors so downstream code (tests,
    # legacy callers) that introspect a trace named "arrow" still finds
    # the metadata it expects. The markers are fully transparent so they
    # don't paint anything on the figure.
    arrow_x, arrow_y, arrow_colors = [], [], []
    for _, r in sub.iterrows():
        # Arrowhead is drawn at the **non-focus end** of the edge (rim).
        neighbor = r["target"] if r["source"] == gene else r["source"]
        nx_, ny_ = pos[neighbor]
        arrow_x.append(nx_)
        arrow_y.append(ny_)
        arrow_colors.append(_DIR_COLOR[r["direction"]])
    arrow_trace = go.Scatter(
        x=arrow_x, y=arrow_y, mode="markers",
        marker=dict(size=1, color=arrow_colors, opacity=0.0),
        hoverinfo="skip", showlegend=False, name="arrow",
    )

    # --- Nodes -------------------------------------------------------------
    # Two separate traces — focus and neighbors — so the focus has its
    # own marker styling (size + accent ring) without disturbing the
    # neighbor labels' radial positioning.
    focus_x, focus_y = pos[gene]
    focus_trace = go.Scatter(
        x=[focus_x], y=[focus_y], mode="markers+text",
        text=[f"<b>{gene}</b>"],
        textposition="middle center",
        textfont=dict(size=13, color="#1f1f1f"),
        marker=dict(size=38, color=_FOCUS_COLOR,
                    line=dict(width=2, color="#ffffff"),
                    opacity=0.95),
        hoverinfo="text", hovertext=[f"<b>{gene}</b> (focus)"],
        showlegend=False, name="focus-node",
    )

    if neighbors:
        n_xs   = [pos[n_][0] for n_ in neighbors]
        n_ys   = [pos[n_][1] for n_ in neighbors]
        # Labels actually rendered just outside the ring via separate
        # text-only trace so the marker isn't pushed off-circle.
        n_textpos = [_label_position(node_angles[n_]) for n_ in neighbors]
        neighbor_trace = go.Scatter(
            x=n_xs, y=n_ys, mode="markers",
            marker=dict(size=15, color=_NEIGH_COLOR,
                        line=dict(width=1.5, color="#ffffff"),
                        opacity=0.92),
            hoverinfo="text",
            hovertext=[f"<b>{n_}</b>" for n_ in neighbors],
            showlegend=False, name="neighbor-nodes",
        )
        label_xs = [_LABEL_RADIUS * np.cos(node_angles[n_]) for n_ in neighbors]
        label_ys = [_LABEL_RADIUS * np.sin(node_angles[n_]) for n_ in neighbors]
        label_trace = go.Scatter(
            x=label_xs, y=label_ys, mode="text",
            text=neighbors,
            textposition=n_textpos,
            textfont=dict(size=11, color="#1f1f1f"),
            hoverinfo="skip", showlegend=False, name="neighbor-labels",
        )
    else:
        neighbor_trace = label_trace = go.Scatter(x=[], y=[], showlegend=False)

    # --- Legend traces (zero-length scatters that show up in the legend) ---
    legend_labels = {
        "Activator": "Activator (+)",
        "Repressor": "Repressor (−)",
        "Unknown":   "Unknown",
    }
    legend_traces = [
        go.Scatter(
            x=[None], y=[None], mode="lines",
            line=dict(color=_DIR_COLOR[d], width=4),
            name=legend_labels[d],
            showlegend=True, hoverinfo="skip",
        )
        for d in ("Activator", "Repressor", "Unknown")
        if d in sub["direction"].unique()
    ]

    fig = go.Figure(data=[
        ring_trace,
        *edge_traces,
        arrow_trace,
        hover_trace,
        neighbor_trace,
        label_trace,
        focus_trace,
        *legend_traces,
    ])

    title = (
        f"<b>{network_type} network for {gene}</b>"
        f"<br><span style='font-size:11px;color:#666'>"
        f"{len(sub)} edges · width ∝ |MI| · color = direction"
        f"</span>"
    )
    fig.update_layout(
        title=dict(text=title, x=0.5, xanchor="center", y=0.97),
        xaxis=dict(
            showticklabels=False, showgrid=False, zeroline=False,
            range=[-_PLOT_PAD, _PLOT_PAD],
            scaleanchor="y", scaleratio=1, constrain="domain",
            fixedrange=True,
        ),
        yaxis=dict(
            showticklabels=False, showgrid=False, zeroline=False,
            range=[-_PLOT_PAD, _PLOT_PAD],
            fixedrange=True,
        ),
        margin=dict(l=20, r=20, t=70, b=20),
        plot_bgcolor=_BG_COLOR,
        paper_bgcolor="#ffffff",
        legend=dict(
            orientation="h", yanchor="bottom", y=-0.02,
            xanchor="center", x=0.5,
            bgcolor="rgba(255,255,255,0)",
        ),
        # Square canvas so the circle stays a circle at any sizing.
        width=720, height=720,
        autosize=False,
        # Real arrow annotations — pixel-sized arrowheads that
        # render crisply regardless of zoom.
        annotations=edge_annotations,
    )
    return fig


def _fmt(value, spec: str = ".3f") -> str:
    """Format a scalar that may be NaN / None as a fixed string."""
    if value is None:
        return "—"
    try:
        f = float(value)
    except (TypeError, ValueError):
        return str(value)
    if np.isnan(f):
        return "—"
    return format(f, spec)
