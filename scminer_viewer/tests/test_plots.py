"""Smoke-test the plot helpers."""

from __future__ import annotations

import plotly.graph_objects as go

from scminer_viewer import load_study
from scminer_viewer.plots import (
    bubble_plot,
    cluster_plot,
    feature_plot,
    heatmap_plot,
    network_plot,
    violin_plot,
)


def test_cluster_plot(fixture_bundle):
    s = load_study(fixture_bundle)
    fig = cluster_plot(s, active_clusters=list(s.clusters.index),
                      dot_size=4, show_labels=True)
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == len(s.clusters)


def test_feature_plot(fixture_bundle):
    s = load_study(fixture_bundle)
    fig = feature_plot(s, gene=s.genes[0], dot_size=4)
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 1


def test_violin_plot(fixture_bundle):
    s = load_study(fixture_bundle)
    fig = violin_plot(s, gene=s.genes[0])
    assert isinstance(fig, go.Figure)
    assert len(fig.data) >= 1


def test_heatmap_plot(fixture_bundle):
    s = load_study(fixture_bundle)
    fig = heatmap_plot(s, genes=list(s.genes[:5]))
    assert isinstance(fig, go.Figure)
    heatmap_traces = [t for t in fig.data if t.type == "heatmap"]
    # One annotation row (cluster colors) + one main heatmap (genes x cells).
    assert len(heatmap_traces) == 2
    # Legend traces — one invisible scatter marker per cluster.
    legend_traces = [t for t in fig.data
                     if t.type == "scatter" and getattr(t, "showlegend", False)]
    assert len(legend_traces) == len(s.clusters)
    # Main heatmap has one row per gene, one column per cell (not per cluster).
    main = next(t for t in heatmap_traces if len(t.y) > 1)
    assert len(main.y) == 5
    assert len(main.x) == s.n_cells


def test_bubble_plot(fixture_bundle):
    s = load_study(fixture_bundle)
    fig = bubble_plot(s, genes=list(s.genes[:5]))
    assert isinstance(fig, go.Figure)


def test_network_plot(fixture_bundle):
    s = load_study(fixture_bundle)
    # Pick a gene that is referenced in the network (deterministic via seed)
    src = s.network_tf["source"].iloc[0]
    fig = network_plot(s, gene=src, network_type="TF")
    assert isinstance(fig, go.Figure)


def test_network_plot_centers_focus_gene(fixture_bundle):
    """Focus gene must sit at the origin; neighbor markers on the ring."""
    import math
    s = load_study(fixture_bundle)
    src = s.network_tf["source"].iloc[0]
    fig = network_plot(s, gene=src, network_type="TF")

    focus_trace = next(t for t in fig.data
                       if getattr(t, "name", "") == "focus-node")
    assert list(focus_trace.x) == [0.0]
    assert list(focus_trace.y) == [0.0]

    neighbor_trace = next(t for t in fig.data
                          if getattr(t, "name", "") == "neighbor-nodes")
    for x, y in zip(neighbor_trace.x, neighbor_trace.y):
        assert math.isclose(math.hypot(x, y), 1.0, abs_tol=1e-9), \
            f"neighbor off the unit circle (r={math.hypot(x, y)})"

    # Figure is square / equal-aspect so the circle stays a circle.
    assert fig.layout.xaxis.scaleanchor == "y"
    assert fig.layout.width == fig.layout.height


def test_network_plot_line_color_matches_arrow_color(fixture_bundle):
    """For every edge: line trace's color must equal the arrow marker's color.

    Each edge gets its own go.Scatter line trace; the arrow trace carries
    per-marker colors. We pair them by position.
    """
    s = load_study(fixture_bundle)
    src = s.network_tf["source"].iloc[0]
    fig = network_plot(s, gene=src, network_type="TF")

    line_traces = [
        t for t in fig.data
        if getattr(t, "mode", None) == "lines"
        and (t.name or "").startswith("edge-")
    ]
    arrow_trace = next(t for t in fig.data if getattr(t, "name", "") == "arrow")
    arrow_colors = list(arrow_trace.marker.color)
    assert len(line_traces) == len(arrow_colors), (
        f"line trace count {len(line_traces)} != arrow count {len(arrow_colors)}"
    )
    for i, (lt, ac) in enumerate(zip(line_traces, arrow_colors)):
        assert lt.line.color.lower() == ac.lower(), (
            f"edge {i}: line color {lt.line.color} != arrow color {ac}"
        )


def test_network_plot_uses_annotation_arrows(fixture_bundle):
    """Each edge must have a real plotly arrow annotation (not just a marker).

    Locks the design choice: arrows are drawn as `layout.annotations`
    with `showarrow=True` and `arrowhead >= 1`, so the head is a crisp
    pixel-sized arrowhead instead of a possibly-symmetric triangle
    marker. Per-annotation color must match the matching edge's color.
    """
    s = load_study(fixture_bundle)
    src = s.network_tf["source"].iloc[0]
    fig = network_plot(s, gene=src, network_type="TF")

    line_traces = [t for t in fig.data
                   if (t.name or "").startswith("edge-")]
    annotations = list(fig.layout.annotations or [])
    arrow_annots = [a for a in annotations if a.showarrow]

    # One arrow annotation per edge.
    assert len(arrow_annots) == len(line_traces), (
        f"arrow annotations {len(arrow_annots)} != line traces "
        f"{len(line_traces)}"
    )
    # Every annotation uses a real arrowhead.
    for a in arrow_annots:
        assert a.arrowhead is not None and a.arrowhead >= 1
        assert a.arrowwidth and a.arrowwidth > 0
        assert a.arrowcolor

    # Color per annotation matches the paired line trace.
    for lt, a in zip(line_traces, arrow_annots):
        assert lt.line.color.lower() == a.arrowcolor.lower(), (
            f"line color {lt.line.color} != arrow color {a.arrowcolor}"
        )


def test_network_plot_arrows_at_neighbor_end(fixture_bundle):
    """Each arrow's tail is at the focus (centre), head at the non-focus end.

    Locks the visual orientation: arrow always lands at the rim next to
    the connected gene, never at the centre. Tail must be at (0, 0).
    """
    import math
    s = load_study(fixture_bundle)
    src = s.network_tf["source"].iloc[0]
    fig = network_plot(s, gene=src, network_type="TF")
    annotations = [a for a in (fig.layout.annotations or []) if a.showarrow]

    for a in annotations:
        # Tail at the focus gene (origin).
        assert math.isclose(a.ax, 0.0, abs_tol=1e-12)
        assert math.isclose(a.ay, 0.0, abs_tol=1e-12)
        # Head out at the rim (on the unit circle, modulo standoff which
        # is in pixels not data units — so the data coordinate stays at r=1).
        assert math.isclose(math.hypot(a.x, a.y), 1.0, abs_tol=1e-9), (
            f"arrow head not at rim (r={math.hypot(a.x, a.y):.4f})"
        )


def test_network_plot_shows_mi_and_direction(fixture_bundle):
    """Hover overlay must surface MI + Direction + correlations."""
    s = load_study(fixture_bundle)
    src = s.network_tf["source"].iloc[0]
    fig = network_plot(s, gene=src, network_type="TF")

    # Find the edge-midpoint hover trace.
    hover_traces = [
        t for t in fig.data
        if getattr(t, "hovertext", None) and "Direction" in (t.hovertext[0] or "")
    ]
    assert hover_traces, "no edge hover trace carrying Direction"
    sample = hover_traces[0].hovertext[0]
    assert "MI:" in sample
    assert "Direction:" in sample
    assert "Spearman:" in sample or "Pearson:" in sample
    assert ("Activator" in sample) or ("Repressor" in sample) or ("Unknown" in sample)

    # Title surfaces |MI| width encoding.
    assert "width" in (fig.layout.title.text or "").lower()
    assert "MI" in (fig.layout.title.text or "")

    # Legend has at least one direction entry visible.
    legend_names = [t.name for t in fig.data
                    if getattr(t, "showlegend", False)]
    assert any(("Activator" in (n or "")) or ("Repressor" in (n or ""))
               for n in legend_names)
