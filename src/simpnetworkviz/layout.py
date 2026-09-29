"""Cluster-aware layouts in physical units (1 unit ≈ 1 inch at scale 1)."""
from __future__ import annotations

import textwrap
from collections import defaultdict
from typing import Dict, List, Tuple

import networkx as nx
import numpy as np

# Geometry constants shared with plotting (inches / points).
HULL_PAD = 0.12          # gap between node edge and cluster outline
BLOCK_GAP = 0.35         # gap between packed clusters
LABEL_FONTSIZE = 14.0    # cluster label font size (pt)
LABEL_LINE = 1.25 * LABEL_FONTSIZE / 72
LABEL_CHAR = 0.56 * LABEL_FONTSIZE / 72
NODE_LABEL_SIZE = 9.0    # node label font size (pt)
NODE_LABEL_WRAP = 16     # characters per line for pathway node labels


def wrap_label(text: str, width: int = 24) -> str:
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False)) or text


def label_width_chars(n_nodes: int) -> int:
    return int(min(34, max(16, 12 + 3 * np.sqrt(n_nodes))))


def node_radii(net) -> Dict[str, float]:
    """Node radius (inches) from ``size_value`` using a sqrt scale."""
    G = net.graph
    out = {}
    for kind in ("pathway", "gene"):
        nodes = [n for n, d in G.nodes(data=True) if d.get("type", "pathway") == kind]
        if not nodes:
            continue
        rmin, rmax = (0.07, 0.22) if kind == "pathway" else (0.05, 0.13)
        if net.kind == "gene-pathway" and kind == "pathway":
            rmin, rmax = 0.10, 0.24
        v = np.sqrt(np.array([max(float(G.nodes[n].get("size_value", 1) or 1), 1) for n in nodes]))
        t = (v - v.min()) / (v.max() - v.min()) if v.max() > v.min() else np.full(len(v), 0.5)
        for n, ti in zip(nodes, t):
            out[n] = float(rmin + (rmax - rmin) * ti)
    return out


def node_label_lines(net, n) -> List[str]:
    d = net.graph.nodes[n]
    text = str(d.get("label", n))
    return [text] if d.get("type") == "gene" else wrap_label(text, NODE_LABEL_WRAP).split("\n")


def collision_radii(net, labelled) -> Dict[str, float]:
    """Radius used to keep nodes apart; grows to fit node labels when shown."""
    out = {}
    char_w = 0.56 * NODE_LABEL_SIZE / 72
    for n, r in net.radius.items():
        if n in labelled:
            lines = node_label_lines(net, n)
            half_w = max(len(l) for l in lines) * char_w / 2 + 0.03
            out[n] = max(r + 0.02, half_w, r + len(lines) * 1.1 * NODE_LABEL_SIZE / 72 / 2 + 0.03)
        else:
            out[n] = r
    return out


def _remove_overlaps(P: np.ndarray, r: np.ndarray, gap: float, iters: int = 200) -> np.ndarray:
    P = P.copy()
    n = len(P)
    if n < 2:
        return P
    rng = np.random.default_rng(0)
    for _ in range(iters):
        moved = False
        d = P[:, None, :] - P[None, :, :]
        dist = np.sqrt((d ** 2).sum(-1)) + np.eye(n)
        need = r[:, None] + r[None, :] + gap
        over = (dist < need) & ~np.eye(n, dtype=bool)
        if not over.any():
            break
        for i, j in zip(*np.nonzero(np.triu(over))):
            vec = P[i] - P[j]
            dd = np.linalg.norm(vec)
            if dd < 1e-9:
                vec = rng.normal(size=2)
                dd = np.linalg.norm(vec)
            push = (need[i, j] - dd) / 2 + 1e-4
            P[i] += vec / dd * push
            P[j] -= vec / dd * push
            moved = True
        if not moved:
            break
    return P


def _local_layout(G: nx.Graph, nodes: List[str], radius: Dict[str, float], seed: int) -> np.ndarray:
    """``radius`` is the collision radius (node + label room)."""
    """Lay out one cluster, centered at the origin, with no overlapping nodes."""
    n = len(nodes)
    r = np.array([radius[x] for x in nodes])
    if n == 1:
        return np.zeros((1, 2))
    sub = G.subgraph(nodes)
    spacing = 2 * r.mean() + 0.16
    if sub.number_of_edges() == 0:
        # Unconnected nodes (e.g. singletons): tidy grid, most significant first.
        cols = int(np.ceil(np.sqrt(n * 1.6)))
        step = 2 * r.max() + 0.14
        P = np.array([[(i % cols) * step, -(i // cols) * step] for i in range(n)], dtype=float)
        return P - P.mean(axis=0)
    if n == 2:
        P = np.array([[-1.0, 0.0], [1.0, 0.0]])
    elif nx.is_connected(sub) and n <= 80:
        P = nx.kamada_kawai_layout(sub, weight=None)
        P = np.array([P[x] for x in nodes])
    else:
        pos = nx.spring_layout(sub, weight="weight", seed=seed, iterations=300, k=1.2 / np.sqrt(n))
        P = np.array([pos[x] for x in nodes])
    P = P - P.mean(axis=0)
    # Rotate so the long axis is horizontal (wide clusters pack better with labels on top).
    if n >= 3:
        _, vecs = np.linalg.eigh(np.cov(P.T))
        P = P @ vecs[:, ::-1]
    # Scale so the median nearest-neighbour distance equals the target spacing.
    d = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1)) + np.eye(n) * 1e9
    nn = np.median(d.min(axis=1))
    P = P * (spacing / nn if nn > 1e-9 else 1.0)
    P = _remove_overlaps(P, r, gap=0.06)
    return P - (P.max(axis=0) + P.min(axis=0)) / 2


def clustered_layout(net, seed: int = 0, aspect: float = 1.5) -> Dict[str, Tuple[float, float]]:
    """Lay out each cluster separately and pack clusters in rows.

    Mirrors Cytoscape AutoAnnotate's "layout network to prevent cluster
    overlap": every cluster gets its own compact region with room for its label.
    Unclustered nodes are placed in a grid after the clusters.
    """
    G = net.graph
    coll = collision_radii(net, getattr(net, "label_nodes", set()))
    members = defaultdict(list)
    for n in G.nodes:
        members[net.groups.get(n, 0)].append(n)

    def sig(n):
        f = G.nodes[n].get("fdr")
        return -np.log10(max(f, 1e-300)) if f is not None and f == f else 0.0

    blocks = []
    for gid in sorted(members, key=lambda g: (g == 0, g)):
        nodes = sorted(members[gid], key=lambda n: (-sig(n), -G.nodes[n].get("size_value", 0), str(n)))
        P = _local_layout(G, nodes, coll, seed)
        r = np.array([coll[x] for x in nodes])
        pad = HULL_PAD if gid else 0.02
        xmin, xmax = (P[:, 0] - r).min() - pad, (P[:, 0] + r).max() + pad
        ymin, ymax = (P[:, 1] - r).min() - pad, (P[:, 1] + r).max() + pad
        w, h = xmax - xmin, ymax - ymin
        label_h = 0.0
        if gid and net.group_labels.get(gid):
            text = wrap_label(net.group_labels[gid], label_width_chars(len(nodes)))
            lines = text.split("\n")
            label_h = len(lines) * LABEL_LINE + 0.08
            w = max(w, max(len(l) for l in lines) * LABEL_CHAR + 0.1)
        blocks.append({"gid": gid, "nodes": nodes, "P": P, "cx": (xmin + xmax) / 2, "cy": (ymin + ymax) / 2,
                       "w": w, "h": h + label_h, "body_h": h})

    total_area = sum((b["w"] + BLOCK_GAP) * (b["h"] + BLOCK_GAP) for b in blocks)
    row_width = max(max(b["w"] for b in blocks), np.sqrt(total_area * aspect))

    rows, cur, cur_w = [], [], 0.0
    for b in blocks:
        if cur and cur_w + b["w"] > row_width:
            rows.append(cur)
            cur, cur_w = [], 0.0
        cur.append(b)
        cur_w += b["w"] + BLOCK_GAP
    if cur:
        rows.append(cur)

    pos = {}
    y_top = 0.0
    for row in rows:
        row_h = max(b["h"] for b in row)
        row_w = sum(b["w"] for b in row) + BLOCK_GAP * (len(row) - 1)
        x = -row_w / 2
        for b in row:
            label_h = b["h"] - b["body_h"]
            # center body horizontally in its slot, top-align below the label
            bx = x + b["w"] / 2
            by = y_top - label_h - b["body_h"] / 2
            for n, p in zip(b["nodes"], b["P"]):
                pos[n] = (float(p[0] - b["cx"] + bx), float(p[1] - b["cy"] + by))
            x += b["w"] + BLOCK_GAP
        y_top -= row_h + BLOCK_GAP
    return pos


def spring_layout(net, seed: int = 0) -> Dict[str, Tuple[float, float]]:
    """Force-directed layout of the whole graph (used for gene-pathway networks)."""
    G = net.graph
    nodes = list(G.nodes)
    n = len(nodes)
    if n == 1:
        return {nodes[0]: (0.0, 0.0)}
    pos = nx.spring_layout(G, weight="weight", seed=seed, iterations=500, k=1.5 / np.sqrt(n))
    P = np.array([pos[x] for x in nodes])
    coll = collision_radii(net, getattr(net, "label_nodes", set()))
    r = np.array([coll[x] for x in nodes])
    P -= P.mean(axis=0)
    span = np.sqrt(n) * (2 * r.mean() + 0.35)
    ext = np.abs(P).max() or 1.0
    P = P / ext * span / 2
    P = _remove_overlaps(P, r, gap=0.08)
    return {x: (float(p[0]), float(p[1])) for x, p in zip(nodes, P)}
