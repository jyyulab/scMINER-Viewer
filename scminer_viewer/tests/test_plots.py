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
    assert len(fig.data) == 1


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
