"""Publication-style static network plots (PNG/PDF/SVG) with matplotlib."""
from __future__ import annotations

from collections import defaultdict
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colormaps, patheffects
from matplotlib.colors import to_rgba
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Circle, Polygon, Wedge
from scipy.spatial import ConvexHull

from .layout import HULL_PAD, LABEL_FONTSIZE, NODE_LABEL_SIZE, NODE_LABEL_WRAP, label_width_chars, wrap_label
from .style import (DIVERGING_CMAP, DOWN_CMAP, SEQ_CMAP, UP_CMAP, condition_color, darken, group_color, hull_colors,
                    hulls_colored, lighten, node_colors)

plt.rcParams["pdf.fonttype"] = 42  # editable text in Illustrator
plt.rcParams["svg.fonttype"] = "none"

#: All figure text uses Arial, regular weight. Liberation Sans (metrically
#: identical to Arial) is the fallback on systems without Arial.
FONT_FAMILY = ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = FONT_FAMILY + [f for f in plt.rcParams["font.sans-serif"] if f not in FONT_FAMILY]
plt.rcParams["font.weight"] = "normal"
plt.rcParams["axes.titleweight"] = "normal"
plt.rcParams["axes.labelweight"] = "normal"

# Font sizes (pt) for figure text other than cluster/node labels (see layout.py).
TITLE_SIZE = 18.0
LEGEND_TITLE_SIZE = 13.0
LEGEND_TEXT_SIZE = 11.0
LEGEND_BAND = 1.35  # inches reserved for the legend row


def hull_polygon(points: np.ndarray, radii: np.ndarray, pad: float = HULL_PAD, n: int = 24, extra=None) -> np.ndarray:
    """Smooth outline around circles: convex hull of each circle grown by ``pad``.

    ``extra``: additional points (e.g. node-label corners) the outline must enclose.
    """
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ring = np.stack([np.cos(ang), np.sin(ang)], axis=1)
    pts = np.concatenate([p + ring * (r + pad) for p, r in zip(points, radii)])
    if extra is not None and len(extra):
        # Round the label corners too so the outline stays smooth.
        small = ring * pad * 0.6
        pts = np.concatenate([pts] + [np.asarray(e) + small for e in extra])
    try:
        hull = ConvexHull(pts)
        return pts[hull.vertices]
    except Exception:  # degenerate (collinear) input
        return pts


def resolve_node_labels(net, node_labels="auto") -> set:
    """Which nodes get a text label (``auto``: small networks, genes, cnet)."""
    G = net.graph
    if node_labels in (None, False, "none"):
        return set()
    if node_labels in (True, "all"):
        return set(G.nodes)
    if node_labels == "top":
        best = {}
        for n, d in G.nodes(data=True):
            gid = net.groups.get(n, 0)
            f = d.get("fdr")
            s = -np.log10(max(f, 1e-300)) if f is not None and f == f else d.get("size_value", 0)
            if gid not in best or s > best[gid][0]:
                best[gid] = (s, n)
        return {n for _, n in best.values()}
    if isinstance(node_labels, (list, set, tuple)):
        return set(node_labels) & set(G.nodes)
    # auto
    n = G.number_of_nodes()
    if net.kind == "pathway":
        return set(G.nodes) if n <= 30 else set()
    if net.kind == "gene-pathway":
        return set(G.nodes)
    return set(G.nodes) if n <= 200 else set()


def _label_box(net, n, size):
    """Corners (inches) of the text box drawn under node ``n``."""
    from .layout import node_label_lines

    x, y = net.pos[n]
    r = net.radius[n]
    lines = node_label_lines(net, n)
    half_w = max(len(l) for l in lines) * 0.56 * size / 72 / 2
    bottom = y - r - 0.02 - len(lines) * 1.1 * size / 72
    return [(x - half_w, y - r), (x + half_w, y - r), (x - half_w, bottom), (x + half_w, bottom)]


def _geometry(net, hulls_on: bool, show_cluster_labels: bool, cluster_label_size: float, margin: float = 0.25,
              labelled=(), node_label_size: float = NODE_LABEL_SIZE):
    """Cluster outlines, label anchors and data extents (inches) of a network."""
    G, pos, rad = net.graph, net.pos, net.radius
    boxes = {n: _label_box(net, n, node_label_size) for n in labelled}
    members = defaultdict(list)
    for n in G.nodes:
        members[net.groups.get(n, 0)].append(n)
    hull_polys, label_specs = {}, []
    xs, ys = [], []
    for n in G.nodes:
        x, y = pos[n]
        r = rad[n]
        xs += [x - r, x + r]
        ys += [y - r, y + r]
    for b in boxes.values():
        xs += [p[0] for p in b]
        ys += [p[1] for p in b]
    for gid, nodes in members.items():
        if gid == 0:
            continue
        P = np.array([pos[n] for n in nodes])
        R = np.array([rad[n] for n in nodes])
        if hulls_on:
            poly = hull_polygon(P, R, extra=[c for n in nodes if n in boxes for c in boxes[n]])
            hull_polys[gid] = poly
            xs += list(poly[:, 0])
            ys += list(poly[:, 1])
            top = poly[:, 1].max()
            cx = (poly[:, 0].max() + poly[:, 0].min()) / 2
        else:
            top = (P[:, 1] + R).max() + 0.05
            cx = P[:, 0].mean()
        text = net.group_labels.get(gid, "")
        if show_cluster_labels and text:
            text = wrap_label(text, label_width_chars(len(nodes)))
            label_specs.append((gid, cx, top + 0.03, text))
            ys.append(top + 0.03 + len(text.split("\n")) * 1.25 * cluster_label_size / 72)
    extent = (min(xs) - margin, max(xs) + margin, min(ys) - margin, max(ys) + margin)
    return hull_polys, label_specs, extent


class _Panel:
    """Everything needed to draw one network into one axes."""

    def __init__(self, net, color_by="auto", node_labels="auto", show_hulls="auto", show_cluster_labels="auto",
                 hull_style="auto", edge_color="#8F8F8F", edge_alpha=0.45, cluster_label_size=LABEL_FONTSIZE,
                 node_label_size=None):
        self.net = net
        self.colors, self.leg = node_colors(net, color_by)
        self.mode = self.leg["mode"]
        hulls_on = (net.layout_style == "clustered") if show_hulls == "auto" else bool(show_hulls)
        if show_cluster_labels == "auto":
            show_cluster_labels = net.layout_style == "clustered"
        self.hull_style = hull_style
        self.hcolors = hull_colors(net, self.mode, hull_style)
        self.to_label = net.label_nodes if node_labels == "auto" else resolve_node_labels(net, node_labels)
        self.hull_polys, self.label_specs, self.extent = _geometry(
            net, hulls_on, show_cluster_labels, cluster_label_size, labelled=self.to_label,
            node_label_size=node_label_size or NODE_LABEL_SIZE)
        x0, x1, y0, y1 = self.extent
        self.W, self.H = x1 - x0, y1 - y0
        self.node_labels = node_labels
        self.edge_color, self.edge_alpha = edge_color, edge_alpha
        self.cluster_label_size, self.node_label_size = cluster_label_size, node_label_size

    def draw(self, ax, f: float):
        net, G, pos, rad, colors = self.net, self.net.graph, self.net.pos, self.net.radius, self.colors
        x0, x1, y0, y1 = self.extent
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ax.axis("off")
        fs = lambda pt: max(pt * f, 3.0)
        lw = lambda w: w * f

        for gid, poly in self.hull_polys.items():
            fill, edge = self.hcolors.get(gid, ("#F2F2F2", "#BBBBBB"))
            ax.add_patch(Polygon(poly, closed=True, facecolor=fill, edgecolor=edge, linewidth=lw(1.3), zorder=0,
                                 joinstyle="round"))

        if G.number_of_edges():
            segs, widths, ecols, alphas = [], [], [], []
            ws = np.array([d.get("weight", 1.0) for _, _, d in G.edges(data=True)], dtype=float)
            wmin, wmax = ws.min(), ws.max()
            for (u, v, d), w in zip(G.edges(data=True), ws):
                segs.append([pos[u], pos[v]])
                t = (w - wmin) / (wmax - wmin) if wmax > wmin else 0.5
                widths.append(lw(0.4 + 2.2 * t) if net.kind == "pathway" else lw(0.3 + 1.2 * t))
                gu, gv = net.groups.get(u, 0), net.groups.get(v, 0)
                # Edges between clusters are drawn fainter so modules stay readable.
                alphas.append(self.edge_alpha if (gu == gv and gu) or net.kind == "gene-pathway" else self.edge_alpha * 0.45)
                if net.kind == "gene-pathway":
                    p = u if G.nodes[u].get("type") == "pathway" else v
                    ecols.append(lighten(colors[p], 0.3))
                else:
                    ecols.append(self.edge_color)
            rgba = [to_rgba(c, a) for c, a in zip(ecols, alphas)]
            ax.add_collection(LineCollection(segs, linewidths=widths, colors=rgba, zorder=1, capstyle="round"))

        order = sorted(G.nodes, key=lambda n: -rad[n])
        wedges = self.leg.get("wedges") or {}
        plain = [n for n in order if n not in wedges]
        circles = [Circle(pos[n], rad[n]) for n in plain]
        fc = [colors[n] for n in plain]
        ec = [darken(colors[n], 0.35) for n in plain]
        ax.add_collection(PatchCollection(circles, facecolors=fc, edgecolors=ec, linewidths=lw(0.7), zorder=2))
        if wedges:
            # Pie nodes: one slice per condition, clockwise from 12 o'clock; empty slices white.
            patches, wfc = [], []
            for n in order:
                if n not in wedges:
                    continue
                k = len(wedges[n])
                for i, c in enumerate(wedges[n]):
                    patches.append(Wedge(pos[n], rad[n], 90 - 360 * (i + 1) / k, 90 - 360 * i / k))
                    wfc.append(c or "#FFFFFF")
            ax.add_collection(PatchCollection(patches, facecolors=wfc, edgecolors="#5A5A5A", linewidths=lw(0.4), zorder=2))
            ax.add_collection(PatchCollection([Circle(pos[n], rad[n]) for n in order if n in wedges], facecolors="none",
                                              edgecolors="#444444", linewidths=lw(0.8), zorder=3))

        stroke = [patheffects.withStroke(linewidth=lw(2.2), foreground="white")]
        colored_hulls = hulls_colored(net, self.mode, self.hull_style)
        for gid, cx, y, text in self.label_specs:
            base = self.hcolors.get(gid, ("", "#555555"))[1]
            col = darken(base, 0.45) if colored_hulls else "#3A3A3A"
            ax.text(cx, y, text, ha="center", va="bottom", fontsize=fs(self.cluster_label_size),
                    color=col, zorder=5, linespacing=1.05, path_effects=stroke)
        size = self.node_label_size or NODE_LABEL_SIZE
        for n in self.to_label:
            d = G.nodes[n]
            x, y = pos[n]
            if d.get("type") == "gene":
                ax.text(x, y - rad[n] - 0.015, d.get("label", n), ha="center", va="top", fontsize=fs(size),
                        color="#222222", zorder=4, path_effects=stroke)
            else:
                ax.text(x, y - rad[n] - 0.02, wrap_label(d.get("label", n), NODE_LABEL_WRAP), ha="center", va="top",
                        fontsize=fs(size), color="#222222", zorder=4, linespacing=1.0, path_effects=stroke)


def plot_network(
    net,
    path: Optional[str] = None,
    color_by: str = "auto",
    node_labels="auto",
    show_hulls="auto",
    show_cluster_labels="auto",
    hull_style: str = "auto",
    title: Optional[str] = None,
    legend: bool = True,
    max_width: float = 16.0,
    max_height: float = 20.0,
    dpi: int = 300,
    edge_color: str = "#8F8F8F",
    edge_alpha: float = 0.45,
    cluster_label_size: float = LABEL_FONTSIZE,
    node_label_size: Optional[float] = None,
    font_family: Optional[str] = None,
):
    """Draw a :class:`~simpnetworkviz.network.Network` and optionally save it.

    Parameters
    ----------
    color_by : "auto", "cluster", "direction" (red up / blue down, shade = FDR),
        "significance" (-log10 FDR), "nes", "logfc" or any numeric node attribute
    node_labels : "auto", "all", "none", "top" (best node per cluster) or a list of node ids
    show_hulls : draw cluster outlines ("auto": on for clustered layouts)
    hull_style : "auto", "cluster" (outline in cluster color) or "neutral" (grey)
    path : output file; the extension decides the format (png, pdf, svg, ...)
    """
    if font_family:
        plt.rcParams["font.family"] = font_family
    panel = _Panel(net, color_by, node_labels, show_hulls, show_cluster_labels, hull_style, edge_color, edge_alpha,
                   cluster_label_size, node_label_size)
    W, H = panel.W, panel.H
    title_band = 0.6 if title else 0.0
    legend_band = LEGEND_BAND if legend else 0.0
    total_h = H + title_band + legend_band
    W_fig = W
    for _ in range(2):  # second pass widens the figure when the legend row needs more room
        f = min(1.0, max_width / W_fig, max_height / total_h)
        fig = plt.figure(figsize=(W_fig * f, total_h * f), facecolor="white")
        ax = fig.add_axes([(W_fig - W) / 2 / W_fig, legend_band / total_h, W / W_fig, H / total_h])
        panel.draw(ax, f)
        if title:
            fig.text(0.5, 1 - (title_band * 0.5) / total_h, title, ha="center", va="center", fontsize=max(TITLE_SIZE * f, 3.0))
        if not legend:
            break
        used = _draw_legend(fig, net, panel.leg, W_fig, legend_band, [0, 0, 1, legend_band / total_h], f)
        if used <= W_fig + 1e-6:
            break
        plt.close(fig)
        W_fig = used
    if path:
        fig.savefig(path, dpi=dpi, facecolor="white")
    return fig


def plot_pathway_gene_network(
    pathway_net,
    gene_net,
    path: Optional[str] = None,
    color_by: str = "auto",
    gene_color_by: str = "auto",
    titles=("Pathway network", "Gene network"),
    panel_labels=("a", "b"),
    arrangement: str = "auto",
    legend: bool = True,
    max_width: float = 18.0,
    max_height: float = 22.0,
    dpi: int = 300,
    **panel_kwargs,
):
    """Two-panel figure: pathway network (a) and its gene network (b).

    Both panels share cluster ids, so each cluster keeps the same outline color
    and label in both panels, making it easy to see which genes drive which
    pathway theme. ``arrangement`` is "horizontal", "vertical" or "auto"
    (whichever gives the more compact figure).
    """
    pa = _Panel(pathway_net, color_by, hull_style="cluster", **panel_kwargs)
    pb = _Panel(gene_net, gene_color_by, hull_style="cluster", **panel_kwargs)
    title_band = 0.6 if any(titles) or any(panel_labels) else 0.0
    legend_band = LEGEND_BAND if legend else 0.0
    gap = 0.5
    hA, hB = pa.H + title_band + legend_band, pb.H + title_band + legend_band
    horiz = (pa.W + gap + pb.W, max(hA, hB))
    vert = (max(pa.W, pb.W), hA + gap + hB)
    if arrangement == "auto":
        # Prefer the arrangement whose shape fits the page box best.
        score = lambda wh: min(max_width / wh[0], max_height / wh[1])
        arrangement = "horizontal" if score(horiz) >= score(vert) else "vertical"
    W, Htot = horiz if arrangement == "horizontal" else vert
    f = min(1.0, max_width / W, max_height / Htot)
    fig = plt.figure(figsize=(W * f, Htot * f), facecolor="white")

    def place(panel, x_left, y_top, slot_w, title, letter, legend_y=None):
        # All positions in inches from the bottom-left; converted to figure fractions.
        h = panel.H + title_band + legend_band
        x_off = x_left + (slot_w - panel.W) / 2
        frac = lambda x, y, w, hh: [x / W, y / Htot, w / W, hh / Htot]
        y_bottom = y_top - h
        legend_y = y_bottom if legend_y is None else legend_y
        ax = fig.add_axes(frac(x_off, y_bottom + legend_band, panel.W, panel.H))
        panel.draw(ax, f)
        if title_band:
            ty = (y_top - title_band / 2) / Htot
            fig.text(x_left / W + 0.06 / W, ty, letter or "", ha="left", va="center", fontsize=max((TITLE_SIZE + 2) * f, 3.0))
            fig.text((x_left + slot_w / 2) / W, ty, title or "", ha="center", va="center", fontsize=max(TITLE_SIZE * f, 3.0))
        if legend:
            _draw_legend(fig, panel.net, panel.leg, panel.W, legend_band, frac(x_off, legend_y, panel.W, legend_band), f)

    if arrangement == "horizontal":
        # Networks top-aligned, legends share the bottom row.
        place(pa, 0.0, Htot, pa.W, titles[0], panel_labels[0], legend_y=0.0)
        place(pb, pa.W + gap, Htot, pb.W, titles[1], panel_labels[1], legend_y=0.0)
    else:
        place(pa, 0.0, Htot, W, titles[0], panel_labels[0])
        place(pb, 0.0, Htot - hA - gap, W, titles[1], panel_labels[1])
    if path:
        fig.savefig(path, dpi=dpi, facecolor="white")
    return fig


def _draw_legend(fig, net, leg, W, band, rect, f):
    """Legend row: node-size key, color scale and (pathway networks) edge key."""
    lax = fig.add_axes(rect)
    lax.set_xlim(0, W)
    lax.set_ylim(0, band)
    lax.set_aspect("equal")
    lax.axis("off")
    fs = lambda pt: max(pt * f, 3.0)
    pt = lambda size: size / 72  # points -> inches
    G = net.graph
    x = 0.3
    y_title = band - 0.08
    y_row = y_title - pt(LEGEND_TITLE_SIZE) * 1.4 - 0.3  # centre line of the keys

    def heading(x, text):
        lax.text(x, y_title, text, fontsize=fs(LEGEND_TITLE_SIZE), va="top", ha="left")

    # size legend -------------------------------------------------------------
    size_nodes = [n for n, d in G.nodes(data=True) if d.get("type") == ("gene" if net.kind == "gene" else "pathway")]
    if size_nodes:
        vals = sorted(size_nodes, key=lambda n: G.nodes[n].get("size_value", 0))
        seen, uniq = set(), []
        for n in [vals[0], vals[len(vals) // 2], vals[-1]]:
            v = int(round(G.nodes[n].get("size_value", 0)))
            if v not in seen:
                seen.add(v)
                uniq.append((n, v))
        heading(x, {"pathway": "Gene-set size", "gene": "# pathways", "gene-pathway": "# genes"}[net.kind])
        x0 = x
        for n, v in uniq:
            r = net.radius[n]
            x += max(r, 0.12)
            lax.add_patch(Circle((x, y_row), r, facecolor="#E6E6E6", edgecolor="#7A7A7A", linewidth=0.8 * f))
            lax.text(x, y_row - 0.26, str(v), fontsize=fs(LEGEND_TEXT_SIZE), ha="center", va="top")
            x += max(r, 0.12) + 0.28
        x = max(x, x0 + len("Gene-set size") * 0.56 * pt(LEGEND_TITLE_SIZE)) + 0.45

    # color legend ------------------------------------------------------------
    mode = leg["mode"]
    grad = np.linspace(0, 1, 128)[None, :]
    bar_w, bar_h = 1.8, 0.17

    def bar(x, cmap, lo, hi, label, ticks, y):
        img = colormaps[cmap](lo + (hi - lo) * grad)
        lax.imshow(img, extent=[x, x + bar_w, y - bar_h / 2, y + bar_h / 2], aspect="auto", zorder=1)
        if label:
            lax.text(x - 0.08, y, label, fontsize=fs(LEGEND_TEXT_SIZE), ha="right", va="center")
        for tx, tl in ticks:
            lax.text(x + tx * bar_w, y - bar_h / 2 - 0.05, tl, fontsize=fs(LEGEND_TEXT_SIZE), ha="center", va="top")

    conds = leg.get("conditions")
    if conds and mode in ("condition", "direction"):
        # Pie key: slice order (clockwise from 12 o'clock) with the condition names.
        heading(x, "Conditions (pie slices)")
        k = len(conds)
        r = 0.2
        cx, cy = x + r + 0.05, y_row
        for i, c in enumerate(conds):
            col = condition_color(i) if mode == "condition" else "#D0D0D0"
            lax.add_patch(Wedge((cx, cy), r, 90 - 360 * (i + 1) / k, 90 - 360 * i / k, facecolor=col,
                                edgecolor="#444444", linewidth=0.6 * f))
            ang = np.deg2rad(90 - 360 * (i + 0.5) / k)
            lax.text(cx + 0.62 * r * np.cos(ang), cy + 0.62 * r * np.sin(ang), str(i + 1), fontsize=fs(LEGEND_TEXT_SIZE - 3),
                     ha="center", va="center", color="#111111")
        x = cx + r + 0.15
        per_col = 3
        col_w = max(len(f"{i + 1} {c}") for i, c in enumerate(conds)) * 0.56 * pt(LEGEND_TEXT_SIZE) + 0.25
        for i, c in enumerate(conds):
            colx = x + (i // per_col) * col_w
            rowy = y_row + (1 - (i % per_col)) * pt(LEGEND_TEXT_SIZE) * 1.5
            if mode == "condition":
                lax.add_patch(Circle((colx + 0.06, rowy), 0.05, facecolor=condition_color(i), edgecolor="#444444", linewidth=0.5 * f))
                lax.text(colx + 0.16, rowy, f"{i + 1}  {c}", fontsize=fs(LEGEND_TEXT_SIZE), va="center", ha="left")
            else:
                lax.text(colx, rowy, f"{i + 1}  {c}", fontsize=fs(LEGEND_TEXT_SIZE), va="center", ha="left")
        x += ((len(conds) - 1) // per_col + 1) * col_w + 0.3
        if mode == "condition":
            heading(x, "Shade: -log10 FDR")
            lax.text(x, y_row, f"light = 0, full = {leg['vmax']:.1f}\nwhite slice = not enriched", fontsize=fs(LEGEND_TEXT_SIZE),
                     va="center", ha="left", color="#555555")
            x += 2.2
    if mode == "direction":
        heading(x, "-log10 FDR")
        x += 0.6
        bar(x, UP_CMAP, 0.35, 0.95, "Up", [], y=y_row + 0.14)
        bar(x, DOWN_CMAP, 0.35, 0.95, "Down", [(0, "0"), (1, f"{leg['vmax']:.1f}")], y=y_row - 0.1)
        x += bar_w + 0.6
    elif mode == "significance":
        heading(x, "-log10 FDR")
        bar(x + 0.1, SEQ_CMAP, 0.25, 0.95, "", [(0, "0"), (1, f"{leg['vmax']:.1f}")], y=y_row)
        x += bar_w + 0.6
    elif mode not in ("cluster", "condition"):
        heading(x, leg.get("label", mode))
        v = leg["vmax"]
        bar(x + 0.1, DIVERGING_CMAP, 0.08, 0.92, "", [(0, f"{-v:.1f}"), (0.5, "0"), (1, f"{v:.1f}")], y=y_row)
        x += bar_w + 0.6
    if any(d.get("highlight") is False for _, d in G.nodes(data=True)):
        heading(x, "Grey nodes")
        lax.text(x, y_row, "library pathways\nnot enriched", fontsize=fs(LEGEND_TEXT_SIZE), va="center", ha="left", color="#555555")
        x += 1.6
    if net.kind == "pathway" and G.number_of_edges():
        heading(x, "Edge width")
        lax.text(x, y_row, "gene-set overlap\n(similarity)", fontsize=fs(LEGEND_TEXT_SIZE), va="center", ha="left",
                 color="#555555")
        x += 1.6
    return x + 0.2  # width (inches) the legend row needs


def plot_cluster_barplot(net, path: Optional[str] = None, max_clusters: int = 20, width: float = 7.0, dpi: int = 300):
    """Bar chart of cluster annotations ranked by best -log10 FDR (companion panel)."""
    ct = net.cluster_table()
    ct = ct[ct["best_fdr"].notna()].head(max_clusters)
    if ct.empty:
        return None
    vals = -np.log10(ct["best_fdr"].clip(lower=1e-300)).values
    fig, ax = plt.subplots(figsize=(width, 0.42 * len(ct) + 1.0))
    ax.barh(range(len(ct)), vals, color=[group_color(g) for g in ct["cluster"]], height=0.7)
    ax.set_yticks(range(len(ct)))
    ax.set_yticklabels([f"{l} ({n})" for l, n in zip(ct["label"], ct["n_nodes"])], fontsize=LEGEND_TEXT_SIZE + 1)
    ax.invert_yaxis()
    ax.set_xlabel("-log10 FDR (best pathway in cluster)", fontsize=LEGEND_TITLE_SIZE)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=LEGEND_TEXT_SIZE)
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=dpi, facecolor="white")
    return fig
