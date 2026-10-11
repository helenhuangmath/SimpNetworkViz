"""Shared color logic so static plots and interactive HTML look the same."""
from __future__ import annotations

from typing import Dict, Tuple

import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import to_hex, to_rgb

# Nature Publishing Group-like palette followed by NEJM/Lancet-like colors.
PALETTE = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
    "#7E6148", "#B09C85", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#EE4C97", "#0072B5",
    "#BC3C29", "#925E9F", "#FDAF91", "#AD002A",
]
# Used for pathway hubs in gene-pathway plots when genes are colored red/blue.
HUB_PALETTE = ["#00A087", "#7876B1", "#E18727", "#4DBBD5", "#925E9F", "#20854E", "#B09C85", "#6F99AD",
               "#EE4C97", "#7E6148", "#91D1C2", "#8491B4"]
NEUTRAL = "#BDBDBD"
CONTEXT = "#ECECEC"  # pathways in a library network that are not enriched
# Qualitative colors for compared conditions (Dark2 + Set1 extras); distinct from the cluster palette.
CONDITION_PALETTE = ["#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#377EB8",
                     "#E41A1C", "#984EA3", "#4DAF4A", "#FF7F00", "#A65628", "#F781BF", "#999999", "#17BECF"]
UP_CMAP, DOWN_CMAP, DIVERGING_CMAP, SEQ_CMAP = "Reds", "Blues", "RdBu_r", "YlOrRd"
COLOR_MODES = ("auto", "cluster", "direction", "significance", "nes", "logfc", "condition")


def lighten(color: str, amount: float) -> str:
    r, g, b = to_rgb(color)
    return to_hex((r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount))


def darken(color: str, amount: float) -> str:
    r, g, b = to_rgb(color)
    return to_hex((r * (1 - amount), g * (1 - amount), b * (1 - amount)))


def group_color(gid: int) -> str:
    return NEUTRAL if gid <= 0 else PALETTE[(gid - 1) % len(PALETTE)]


def _vals(net, attr):
    return {n: d.get(attr) for n, d in net.graph.nodes(data=True) if d.get(attr) is not None and not pd.isna(d.get(attr))}


def condition_color(i: int) -> str:
    return CONDITION_PALETTE[i % len(CONDITION_PALETTE)]


def resolve_color_by(net, color_by: str = "auto") -> str:
    if color_by not in COLOR_MODES and not _vals(net, color_by):
        raise ValueError(f"color_by must be one of {COLOR_MODES} or a numeric node attribute")
    if color_by == "condition" and not net.conditions:
        raise ValueError("color_by='condition' needs a comparison network (compare_conditions)")
    if color_by != "auto":
        return color_by
    if net.kind in ("gene", "gene-pathway"):
        if _vals(net, "logfc"):
            return "logfc"
        return "condition" if net.conditions else "cluster"
    if net.conditions:
        return "condition"
    dirs = {d.get("direction") for _, d in net.graph.nodes(data=True)}
    if dirs & {"up", "down"}:
        return "direction"
    if _vals(net, "nes"):
        return "nes"
    return "cluster"


def _wedge_colors(net, mode: str, vmax: float):
    """Per-condition slice colors for comparison networks: ``{node: [hex | None]}``.

    ``mode="condition"``: slice = condition color, darker with significance.
    ``mode="direction"``: slice = red (up) / blue (down) shaded by -log10 FDR.
    Slices of conditions in which the pathway is not enriched stay ``None`` (white).
    """
    conds = net.conditions
    out = {}
    for n, d in net.graph.nodes(data=True):
        if d.get("highlight") is False:
            continue
        slices = []
        if d.get("type") == "pathway":
            for i, c in enumerate(conds):
                f = d.get(f"fdr__{c}")
                active = d.get("conditions") is None or c in str(d.get("conditions", "")).split(";")
                if f is None or pd.isna(f) or not active:
                    slices.append(None)
                    continue
                t = min(_nlog(f) / vmax, 1.0) if vmax else 1.0
                direction = d.get(f"direction__{c}", "")
                if mode == "direction" and direction in ("up", "down"):
                    cmap = UP_CMAP if direction == "up" else DOWN_CMAP
                    slices.append(to_hex(colormaps[cmap](0.35 + 0.6 * t)))
                elif mode == "direction":
                    slices.append(NEUTRAL)
                else:
                    slices.append(lighten(condition_color(i), 0.55 * (1 - t)))
        else:  # gene: slice = condition in which the gene is a hit
            hit = set(str(d.get("conditions", "")).split(";"))
            slices = [condition_color(i) if c in hit else None for i, c in enumerate(conds)]
        if any(slices):
            out[n] = slices
    return out


def _nlog(x) -> float:
    return float(-np.log10(max(float(x), 1e-300)))


def node_colors(net, color_by: str = "auto") -> Tuple[Dict[str, str], dict]:
    """Return ``{node: hex color}`` and a legend description."""
    mode = resolve_color_by(net, color_by)
    G = net.graph
    colors: Dict[str, str] = {}
    legend: dict = {"mode": mode}
    cluster_col = {n: group_color(net.groups.get(n, 0)) for n in G}

    if mode == "condition":
        sig = {n: _nlog(v) for n, v in _vals(net, "fdr").items()}
        vmax = max(max(sig.values()) if sig else 1.0, 1.0)
        legend.update(conditions=list(net.conditions), vmin=0.0, vmax=vmax, label="-log10 FDR")
        wedges = _wedge_colors(net, "condition", vmax)
        legend["wedges"] = wedges
        for n, d in G.nodes(data=True):
            w = wedges.get(n)
            if w and any(w):
                # Base color = the most significant condition (used where pies cannot be drawn).
                best = min((i for i, c in enumerate(w) if c), key=lambda i: d.get(f"fdr__{net.conditions[i]}") or 0.0)
                colors[n] = w[best]
            elif d.get("type") == "pathway" and net.kind == "gene-pathway":
                colors[n] = cluster_col[n]
            else:
                colors[n] = CONTEXT if d.get("highlight") is False else NEUTRAL
        return colors, legend
    if mode == "cluster":
        colors = {n: (CONTEXT if G.nodes[n].get("highlight") is False else c) for n, c in cluster_col.items()}
    elif mode in ("direction", "significance"):
        sig = {n: _nlog(v) for n, v in _vals(net, "fdr").items()}
        vmax = max(sig.values()) if sig else 1.0
        vmax = max(vmax, 1.0)
        legend.update(vmin=0.0, vmax=vmax, label="-log10 FDR")
        if mode == "direction" and net.conditions:
            legend["conditions"] = list(net.conditions)
            legend["wedges"] = _wedge_colors(net, "direction", vmax)
        for n, d in G.nodes(data=True):
            t = sig.get(n, vmax * 0.5) / vmax
            if d.get("highlight") is False:
                colors[n] = CONTEXT
            elif mode == "significance":
                colors[n] = to_hex(colormaps[SEQ_CMAP](0.25 + 0.7 * t)) if n in sig else cluster_col[n]
            else:
                direction = d.get("direction")
                if direction == "up":
                    colors[n] = to_hex(colormaps[UP_CMAP](0.35 + 0.6 * t))
                elif direction == "down":
                    colors[n] = to_hex(colormaps[DOWN_CMAP](0.35 + 0.6 * t))
                elif direction == "mixed":
                    colors[n] = "#B48EAD"  # purple: up in some conditions, down in others
                else:
                    colors[n] = NEUTRAL if d.get("type") == "pathway" else cluster_col[n]
    else:
        attr = mode
        vals = _vals(net, attr)
        arr = np.abs(np.array(list(vals.values()), dtype=float)) if vals else np.array([1.0])
        vmax = float(np.percentile(arr, 95)) or 1.0
        legend.update(vmin=-vmax, vmax=vmax, label={"logfc": "log2 fold change", "nes": "NES"}.get(attr, attr))
        cmap = colormaps[DIVERGING_CMAP]
        for n, d in G.nodes(data=True):
            if d.get("highlight") is False:
                colors[n] = CONTEXT
            elif n in vals:
                t = (np.clip(vals[n], -vmax, vmax) + vmax) / (2 * vmax)
                colors[n] = to_hex(cmap(0.08 + 0.84 * t))
            elif d.get("type") == "pathway" and net.kind == "gene-pathway":
                gid = net.groups.get(n, 0)
                colors[n] = "#8C8C8C" if gid <= 0 else HUB_PALETTE[(gid - 1) % len(HUB_PALETTE)]
            else:
                colors[n] = "#E0E0E0"
    return colors, legend


def hulls_colored(net, mode: str, style: str = "auto") -> bool:
    """Whether cluster outlines use cluster colors (else neutral grey)."""
    # Cluster colors by default, so pathway and gene figures of the same
    # clusters match; "neutral" gives the classic grey AutoAnnotate look.
    return style != "neutral"


def hull_colors(net, mode: str, style: str = "auto") -> Dict[int, Tuple[str, str]]:
    """``{group: (fill, edge)}`` for cluster outlines.

    ``style``: "auto"/"cluster" (cluster colors) or "neutral" (grey).
    """
    out = {}
    colored = hulls_colored(net, mode, style)
    # In a library network, clusters without any enriched pathway stay grey.
    has_context = any(d.get("highlight") is False for _, d in net.graph.nodes(data=True))
    lit = {g for n, g in net.groups.items() if net.graph.nodes[n].get("highlight") is not False} if has_context else None
    for gid in net.group_labels:
        if colored and (lit is None or gid in lit):
            base = group_color(gid)
            out[gid] = (lighten(base, 0.86), lighten(base, 0.25))
        else:
            out[gid] = ("#F3F1EC", "#B7B2A6")
    return out
