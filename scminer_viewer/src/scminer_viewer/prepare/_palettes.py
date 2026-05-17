"""Ported ggsci color palettes for cluster auto-fill.

R's `ggsci` package ships ~16 journal-themed palettes; the R port of
`prepare_study()` lets users pick one by name (default `npg`). The
Python prep doesn't pull ggsci in, so these hex tables are
hand-copied from the ggsci 4.2.0 sources to match the bytes the R
writer produces for the same `cluster_palette` choice.

Each palette is the **alpha-stripped** form (#RRGGBB only) — the R
writer also strips the alpha byte that ggsci's `pal_*()` adds.
"""

from __future__ import annotations

# Color table per ggsci palette. Lengths match what the corresponding
# ggsci::pal_*() returns when called with no `n` argument.
PALETTES: dict[str, list[str]] = {
    "npg": [
        "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
        "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
    ],
    "aaas": [
        "#3B4992", "#EE0000", "#008B45", "#631879", "#008280",
        "#BB0021", "#5F559B", "#A20056", "#808180", "#1B1919",
    ],
    "lancet": [
        "#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F",
        "#FDAF91", "#AD002A", "#ADB6B6", "#1B1919",
    ],
    "nejm": [
        "#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1",
        "#6F99AD", "#FFDC91", "#EE4C97",
    ],
    "jama": [
        "#374E55", "#DF8F44", "#00A1D5", "#B24745", "#79AF97",
        "#6A6599", "#80796B",
    ],
    "jco": [
        "#0073C2", "#EFC000", "#868686", "#CD534C", "#7AA6DC",
        "#003C67", "#8F7700", "#3B3B3B", "#A73030", "#4A6990",
    ],
    "ucscgb": [
        "#FF0000", "#FF9900", "#FFCC00", "#00FF00", "#6699FF",
        "#CC33FF", "#99991E", "#999999", "#FF00CC", "#CC0000",
        "#FFCCCC", "#FFFF00", "#CCFF00", "#358000", "#0000CC",
        "#99CCFF", "#00FFFF", "#CCFFFF", "#9900CC", "#CC99FF",
        "#996600", "#666600", "#666666", "#CCCCCC", "#79CC3D",
        "#CCCC99",
    ],
    "d3": [
        "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
        "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
    ],
    "locuszoom": [
        "#D43F3A", "#EEA236", "#5CB85C", "#46B8DA", "#357EBD",
        "#9632B8", "#B8B8B8",
    ],
    "igv": [
        "#5050FF", "#CE3D32", "#749B58", "#F0E685", "#466983",
        "#BA6338", "#5DB1DD", "#802268", "#6BD76B", "#D595A7",
        "#924822", "#837B8D", "#C75127", "#D58F5C", "#7A65A5",
        "#E4AF69", "#3B1B53", "#CDDEB7", "#612A79", "#AE1F63",
        "#E7C76F", "#5A655E", "#CC9900", "#99CC00", "#A9A9A9",
        "#CC9900",
    ],
    "uchicago": [
        "#800000", "#767676", "#FFA319", "#8A9045", "#155F83",
        "#C16622", "#8F3931", "#58593F", "#350E20",
    ],
    "startrek": [
        "#CC0C00", "#5C88DA", "#84BD00", "#FFCD00", "#7C878E",
        "#00B5E2", "#00AF66",
    ],
    "tron": [
        "#FF410D", "#6EE2FF", "#F7C530", "#95CC5E", "#D0DFE6",
        "#F79D1E", "#748AA6",
    ],
    "futurama": [
        "#FF6F00", "#C71000", "#008EA0", "#8A4198", "#5A9599",
        "#FF6348", "#84D7E1", "#FF95A8", "#3D3B25", "#ADE2D0",
        "#1A5354", "#3F4041",
    ],
    "rickandmorty": [
        "#FAFD7C", "#82491E", "#24325F", "#B7E4F9", "#FB6467",
        "#526E2D", "#E762D7", "#E89242", "#FAE48B", "#A6EEE6",
        "#917C5D", "#69C8EC",
    ],
    "simpsons": [
        "#FED439", "#709AE1", "#8A9197", "#D2AF81", "#FD7446",
        "#D5E4A2", "#197EC0", "#F05C3B", "#46732E", "#71D0F5",
        "#370335", "#075149", "#C80813", "#91331F", "#1A9993",
        "#FD8CC1",
    ],
}

DEFAULT_PALETTE = "npg"


def cluster_colors(n: int, palette: str = DEFAULT_PALETTE) -> list[str]:
    """Return `n` hex colors from the named palette.

    Unknown palette names fall back to `npg` with a stderr warning
    (matching the R behavior).
    """
    if n <= 0:
        return []
    name = palette.lower() if palette else DEFAULT_PALETTE
    colors = PALETTES.get(name)
    if colors is None:
        import warnings
        warnings.warn(
            f"Unknown cluster palette '{palette}'; falling back to npg. "
            f"Valid names: {', '.join(sorted(PALETTES))}",
            UserWarning,
            stacklevel=2,
        )
        colors = PALETTES[DEFAULT_PALETTE]
    # Recycle if the requested count exceeds the palette length, just
    # like R's rep_len().
    return [colors[i % len(colors)] for i in range(n)]
