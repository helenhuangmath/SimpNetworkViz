"""A fixed 2-D map of a whole gene-set library, with enrichment results overlaid.

The map is computed once per library from gene-set similarity (UMAP when the
``umap-learn`` package is installed, otherwise t-SNE from scikit-learn),
clustered into regions with automatic labels, and saved as a CSV. Every later
analysis is drawn on the same coordinates, so results from different
experiments, conditions or cell clusters are directly comparable: the question
"where in the pathway landscape do my hits fall?" has the same answer space
every time.
"""
from __future__ import annotations

import warnings
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import sparse

from .genesets import GeneSets, normalize_gene
from .network import clean_pathway_name

MAP_COLUMNS = ["pathway", "x", "y", "set_size", "cluster", "cluster_label"]


# ---------------------------------------------------------------------------
# Building the map
# ---------------------------------------------------------------------------


def _tfidf_matrix(gene_sets: GeneSets):
    names = list(gene_sets)
    genes = sorted({g for v in gene_sets.values() for g in v})
    gi = {g: i for i, g in enumerate(genes)}
    rows, cols = [], []
    for r, name in enumerate(names):
        for g in set(gene_sets[name]):
            rows.append(r)
            cols.append(gi[g])
    M = sparse.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=(len(names), len(genes)))
    df = np.asarray(M.sum(axis=0)).ravel()
    idf = np.log((1 + len(names)) / (1 + df)) + 1.0  # genes in few sets weigh more
    M = M.multiply(idf.astype(np.float32)).tocsr()
    norm = np.sqrt(np.asarray(M.multiply(M).sum(axis=1)).ravel())
    M = sparse.diags(1.0 / np.maximum(norm, 1e-12)) @ M
    return names, M.tocsr()


def _svd(M, k: int, seed: int) -> np.ndarray:
    from scipy.sparse.linalg import svds

    k = int(min(k, min(M.shape) - 1))
    rng = np.random.default_rng(seed)
    U, S, _ = svds(M, k=k, v0=rng.standard_normal(min(M.shape)))
    order = np.argsort(-S)
    X = (U[:, order] * S[order]).astype(np.float32)
    return X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)


def _knn(X: np.ndarray, k: int, chunk: int = 1024):
    """k nearest neighbours by cosine similarity (rows of ``X`` are unit length)."""
    n = len(X)
    idx = np.empty((n, k), dtype=np.int64)
    sim = np.empty((n, k), dtype=np.float32)
    for s in range(0, n, chunk):
        S = X[s:s + chunk] @ X.T
        S[np.arange(S.shape[0]), np.arange(s, s + S.shape[0])] = -np.inf
        part = np.argpartition(-S, k, axis=1)[:, :k]
        vals = np.take_along_axis(S, part, axis=1)
        o = np.argsort(-vals, axis=1)
        idx[s:s + chunk] = np.take_along_axis(part, o, axis=1)
        sim[s:s + chunk] = np.take_along_axis(vals, o, axis=1)
    return idx, np.clip(sim, 0, 1)


def _knn_graph(idx, sim):
    n, k = idx.shape
    rows = np.repeat(np.arange(n), k)
    W = sparse.csr_matrix((sim.ravel(), (rows, idx.ravel())), shape=(n, n))
    return W.maximum(W.T)


def _import_umap():
    """``import umap`` that survives a broken TensorFlow install.

    ``umap`` imports TensorFlow at import time for its optional parametric
    variant; on some machines TensorFlow is installed but fails to load (for
    example an x86 build on Apple silicon), which would otherwise make plain
    UMAP unavailable as well. Blocking the TensorFlow import for that one
    import is harmless: parametric UMAP is not used here.
    """
    import sys

    try:
        import umap
        return umap
    except Exception:
        if "tensorflow" in sys.modules or "umap" in sys.modules:
            raise
        sys.modules["tensorflow"] = None  # makes ``import tensorflow`` raise ImportError
        try:
            import umap
            return umap
        finally:
            sys.modules.pop("tensorflow", None)


def _embed(X: np.ndarray, method: str, n_neighbors: int, min_dist: float, seed: int):
    if method in ("auto", "umap"):
        try:
            umap = _import_umap()

            emb = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric="cosine", random_state=seed).fit_transform(X)
            return np.asarray(emb, dtype=float), "umap"
        except Exception as exc:  # not installed, or a broken native build
            if method == "umap":
                raise ImportError(f"UMAP is not available ({exc}). pip install umap-learn, or use method='tsne'/'spectral'.") from exc
    if method in ("auto", "tsne"):
        try:
            from sklearn.manifold import TSNE

            perp = int(min(30, max(5, len(X) // 4)))
            emb = TSNE(n_components=2, metric="cosine", init="pca", perplexity=perp, random_state=seed).fit_transform(X)
            return np.asarray(emb, dtype=float), "tsne"
        except Exception as exc:
            if method == "tsne":
                raise ImportError(f"t-SNE needs scikit-learn ({exc}).") from exc
    raise ImportError("Building a pathway map needs umap-learn (recommended) or scikit-learn (t-SNE): "
                      "pip install umap-learn   or   pip install 'simpnetworkviz[map]'")


def _map_clusters(X: np.ndarray, n_neighbors: int, resolution: float, seed: int) -> np.ndarray:
    import networkx as nx

    idx, sim = _knn(X, n_neighbors)
    G = nx.Graph()
    G.add_nodes_from(range(len(X)))
    for i in range(len(X)):
        for j, s in zip(idx[i], sim[i]):
            if s > 0:
                G.add_edge(int(i), int(j), weight=float(s))
    comms = nx.community.louvain_communities(G, weight="weight", resolution=resolution, seed=seed)
    comms = sorted(comms, key=len, reverse=True)
    lab = np.zeros(len(X), dtype=int)
    for c, members in enumerate(comms, start=1):
        for m in members:
            lab[m] = c
    return lab


def build_pathway_map(
    gene_sets: GeneSets,
    method: str = "auto",
    min_size: int = 10,
    max_size: int = 500,
    n_svd: int = 50,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    cluster_resolution: float = 1.0,
    min_cluster_label: int = 5,
    label_words: int = 3,
    seed: int = 0,
) -> pd.DataFrame:
    """Embed every gene set of a library in 2-D and group the map into labelled regions.

    Gene sets become TF-IDF-weighted gene vectors, reduced to ``n_svd`` dimensions
    by truncated SVD, then embedded with ``method``: ``"umap"`` (needs
    ``umap-learn``), ``"tsne"`` (needs scikit-learn), ``"spectral"`` (no extra
    dependency) or ``"auto"`` (the first available in that order). Regions are
    Louvain communities of the k-nearest-neighbour graph, labelled like network
    clusters. Returns a table with :data:`MAP_COLUMNS`; ``attrs["method"]`` records
    the embedding used. Build it once per library and save it with :func:`save_map`.
    """
    gs = {k: v for k, v in gene_sets.items() if min_size <= len(set(v)) <= max_size}
    if len(gs) < 5:
        raise ValueError("Need at least 5 gene sets within the size limits to build a map.")
    names, M = _tfidf_matrix(gs)
    X = _svd(M, n_svd, seed)
    k = int(min(n_neighbors, len(names) - 1))
    emb, used = _embed(X, method, k, min_dist, seed)
    clusters = _map_clusters(X, k, cluster_resolution, seed)
    df = pd.DataFrame({"pathway": names, "x": emb[:, 0], "y": emb[:, 1],
                       "set_size": [len(set(gs[n])) for n in names], "cluster": clusters})
    labels = region_labels(df, min_cluster_label, label_words)
    df["cluster_label"] = df["cluster"].map(labels).fillna("")
    df.loc[df["cluster_label"] == "", "cluster"] = 0
    df.attrs["method"] = used
    return df[MAP_COLUMNS]


# ---------------------------------------------------------------------------
# One map per category (MSigDB collections, several GMT files, custom groups)
# ---------------------------------------------------------------------------

MAP_METHODS = ("tsne", "umap")


def split_by_category(gene_sets: GeneSets, categories: Optional[Dict[str, str]] = None,
                      min_sets: int = 1) -> Dict[str, GeneSets]:
    """Split a library into ``{category: gene_sets}``.

    ``categories`` (``{pathway: category}``) takes precedence; otherwise the
    category is read from the pathway name with
    :func:`~simpnetworkviz.network.pathway_category`, which knows the MSigDB
    collection prefixes (HALLMARK, GOBP/GOCC/GOMF, KEGG, REACTOME, WP, PID,
    BIOCARTA, HP, MIR ...) and the ``"Library: pathway"`` names produced by
    :func:`~simpnetworkviz.genesets.load_gene_sets` for several sources.
    Categories with fewer than ``min_sets`` sets are dropped. Categories are
    returned largest first.
    """
    from .network import pathway_category

    out: Dict[str, GeneSets] = defaultdict(dict)
    for name, genes in gene_sets.items():
        out[pathway_category(name, categories)][name] = genes
    kept = {c: gs for c, gs in out.items() if len(gs) >= min_sets}
    return dict(sorted(kept.items(), key=lambda kv: -len(kv[1])))


def build_category_maps(
    gene_sets: GeneSets,
    methods: Sequence[str] = MAP_METHODS,
    categories: Optional[Dict[str, str]] = None,
    min_sets: int = 20,
    **map_kwargs,
) -> Dict[str, Dict[str, pd.DataFrame]]:
    """One pathway map per category and embedding method.

    The library is split with :func:`split_by_category` (MSigDB collections by
    default, e.g. a full ``msigdb.v2025.1.Hs.symbols.gmt`` gives Hallmark, GO BP,
    GO CC, GO MF, KEGG, Reactome, WikiPathways ... as separate maps) and
    :func:`build_pathway_map` is run on each category with every method in
    ``methods`` (``"tsne"`` and ``"umap"`` by default). Categories with fewer
    than ``min_sets`` sets after the size filter are skipped with a warning,
    as is a method whose package is not installed. ``map_kwargs`` are passed
    to :func:`build_pathway_map` (``min_size``, ``max_size``, ``n_neighbors``,
    ``cluster_resolution``, ``seed`` ...).

    Returns ``{category: {method: map}}``. Within a category the t-SNE and UMAP
    maps share their regions and labels (regions are clustered in similarity
    space), so only the coordinates differ.
    """
    min_size = map_kwargs.get("min_size", 10)
    max_size = map_kwargs.get("max_size", 500)
    min_sets = max(int(min_sets), 5)  # build_pathway_map needs at least 5 sets
    groups = split_by_category(gene_sets, categories, min_sets=1)
    maps: Dict[str, Dict[str, pd.DataFrame]] = {}
    unavailable = set()
    for cat, gs in groups.items():
        n_ok = sum(1 for v in gs.values() if min_size <= len(set(v)) <= max_size)
        if n_ok < min_sets:
            warnings.warn(f"Category {cat!r}: {n_ok} gene sets within {min_size}-{max_size} genes (< {min_sets}), skipped.",
                          stacklevel=2)
            continue
        maps[cat] = {}
        for method in methods:
            if method in unavailable:
                continue
            try:
                pmap = build_pathway_map(gs, method=method, **map_kwargs)
            except ImportError as exc:
                warnings.warn(f"{method}: {exc}", stacklevel=2)
                unavailable.add(method)
                continue
            pmap.attrs["category"] = cat
            maps[cat][method] = pmap
        if not maps[cat]:
            del maps[cat]
    if not maps:
        raise ValueError("No category map could be built: check the gene-set names, the size limits and that "
                         "umap-learn or scikit-learn is installed.")
    return maps


def category_slug(category: str) -> str:
    """File-name-safe form of a category name (``"GO biological process"`` -> ``"go_biological_process"``)."""
    import re

    slug = re.sub(r"[^A-Za-z0-9]+", "_", str(category)).strip("_").lower()
    return slug or "category"


def save_category_maps(maps: Dict[str, Dict[str, pd.DataFrame]], outdir, plot: bool = True,
                       dpi: int = 300, formats: Sequence[str] = ("pdf",)) -> pd.DataFrame:
    """Write every map of :func:`build_category_maps` to ``outdir``.

    Files are ``<category>_map_<method>.csv`` (reuse with ``--map``) and, with
    ``plot=True``, ``<category>_map_<method>.<format>``. An index table
    ``category_maps.csv`` (category, method, n_sets, n_regions, csv, figure) is
    written too and returned.
    """
    import os

    os.makedirs(outdir, exist_ok=True)
    rows = []
    for cat, per_method in maps.items():
        for method, pmap in per_method.items():
            stem = os.path.join(outdir, f"{category_slug(cat)}_map_{method}")
            save_map(pmap, stem + ".csv")
            from .maphtml import write_map_html

            write_map_html(pmap, stem + ".html", title=f"{cat} ({method.upper()})")
            fig_path = ""
            if plot:
                import matplotlib.pyplot as plt

                fig = plot_pathway_map(pmap, None, None, title=f"{cat} ({method.upper()})", dpi=dpi)
                for f in formats:
                    fig.savefig(f"{stem}.{f}", dpi=dpi, facecolor="white", bbox_inches="tight")
                    fig_path = fig_path or f"{stem}.{f}"
                plt.close(fig)
            rows.append({"category": cat, "method": method, "n_sets": len(pmap),
                         "n_regions": int(pmap.loc[pmap["cluster"] > 0, "cluster"].nunique()),
                         "csv": stem + ".csv", "figure": fig_path})
    index = pd.DataFrame(rows, columns=["category", "method", "n_sets", "n_regions", "csv", "figure"])
    index.to_csv(os.path.join(outdir, "category_maps.csv"), index=False)
    return index


def plot_category_maps(maps: Dict[str, Dict[str, pd.DataFrame]], path=None, table: Optional[pd.DataFrame] = None,
                       fdr_cutoff: Optional[float] = 0.05, panel: float = 4.5, dpi: int = 300, label_top: int = 0):
    """One figure with a panel per category (rows) and embedding method (columns).

    ``table`` (an enrichment or comparison table) highlights the enriched
    pathways on each panel, as :func:`plot_pathway_map` does on a single map;
    ``label_top`` names the most significant hits per panel.
    """
    import matplotlib.pyplot as plt

    cats = list(maps)
    methods = sorted({m for per in maps.values() for m in per}, key=lambda m: (m != "tsne", m))
    nrow, ncol = len(cats), len(methods)
    fig, axes = plt.subplots(nrow, ncol, figsize=(panel * ncol, panel * nrow * 0.95), squeeze=False, facecolor="white")
    for i, cat in enumerate(cats):
        for j, method in enumerate(methods):
            ax = axes[i][j]
            pmap = maps[cat].get(method)
            if pmap is None:
                ax.axis("off")
                continue
            plot_pathway_map(pmap, table, None, fdr_cutoff=fdr_cutoff, label_top=label_top, title=cat,
                             width=panel, dpi=dpi, ax=ax)
            if table is not None and len(ax.get_legend_handles_labels()[0]) and (i, j) != (0, ncol - 1):
                ax.get_legend().remove()  # one legend per figure is enough
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=dpi, facecolor="white", bbox_inches="tight")
    return fig


def _repel_labels(ax, texts, iterations: int = 300, pad_px: float = 3.0, damping: float = 0.7) -> None:
    """Nudge overlapping text boxes apart (display space), keeping them inside the axes.

    A small, dependency-free version of what ``adjustText`` does: on every pass
    each pair of overlapping boxes is pushed apart along the axis with the
    smaller overlap (vertical moves preferred, so labels stay near their
    region). Annotations keep their arrow anchored to the original point; only
    the text box is measured, not the arrow.
    """
    from matplotlib.text import Text

    texts = [t for t in texts if t.get_visible()]
    if len(texts) < 2:
        return
    fig = ax.figure
    fig.canvas.draw()  # finalize autoscaling and aspect so display coordinates are the real ones
    renderer = fig.canvas.get_renderer()
    ax_box = ax.get_window_extent(renderer)
    inv = ax.transData.inverted()
    pt = 72.0 / fig.dpi  # pixels -> points, for annotation offsets

    def bbox(t):
        return Text.get_window_extent(t, renderer).padded(pad_px)

    for _ in range(iterations):
        boxes = [bbox(t) for t in texts]
        shifts = np.zeros((len(texts), 2))
        moved = False
        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                a, b = boxes[i], boxes[j]
                if not a.overlaps(b):
                    continue
                dx = min(a.x1 - b.x0, b.x1 - a.x0)
                dy = min(a.y1 - b.y0, b.y1 - a.y0)
                if dy * 1.5 <= dx:
                    sign = 1 if (a.y0 + a.y1) >= (b.y0 + b.y1) else -1
                    shifts[i, 1] += sign * dy / 2
                    shifts[j, 1] -= sign * dy / 2
                else:
                    sign = 1 if (a.x0 + a.x1) >= (b.x0 + b.x1) else -1
                    shifts[i, 0] += sign * dx / 2
                    shifts[j, 0] -= sign * dx / 2
                moved = True
        if not moved:
            break
        for t, (sx, sy), b in zip(texts, shifts * damping, boxes):
            if sx == 0 and sy == 0:
                continue
            sx += min(0.0, ax_box.x1 - (b.x1 + sx)) + max(0.0, ax_box.x0 - (b.x0 + sx))
            sy += min(0.0, ax_box.y1 - (b.y1 + sy)) + max(0.0, ax_box.y0 - (b.y0 + sy))
            if hasattr(t, "xyann"):
                t.xyann = (t.xyann[0] + sx * pt, t.xyann[1] + sy * pt)
            else:
                X, Y = ax.transData.transform(t.get_position())
                t.set_position(inv.transform((X + sx, Y + sy)))


def _words(name: str) -> Dict[str, float]:
    """``{word: relative position in the name}`` for the informative words of a pathway name."""
    import re

    from .network import STOPWORDS, _KEEP_UPPER

    toks = [t for t in re.split(r"[\s,/()\[\]:;]+", clean_pathway_name(name)) if t]
    out: Dict[str, float] = {}
    for i, t in enumerate(toks):
        key = t.lower().strip("-")
        if key and key not in STOPWORDS and len(key) > 2 and not key.isdigit():
            out.setdefault(t if t.upper() in _KEEP_UPPER else key, i / max(len(toks) - 1, 1))
    return out


def region_labels(pmap: pd.DataFrame, min_size: int = 5, max_words: int = 3) -> Dict[int, str]:
    """Label map regions with words that are frequent inside the region and rare outside.

    Plain word frequency (as used for small network clusters) picks generic GO
    words such as "response" or "signaling" for regions of hundreds of terms;
    weighting each word by its inverse frequency across the whole library gives
    the words that set the region apart.
    """
    words = {p: _words(p) for p in pmap["pathway"]}
    df_all: Counter = Counter()
    for ws in words.values():
        df_all.update({w.lower() for w in ws})
    N = len(words)
    labels = {}
    for c, grp in pmap[pmap["cluster"] > 0].groupby("cluster"):
        if len(grp) < min_size:
            continue
        inside: Counter = Counter()
        display, pos = {}, defaultdict(list)
        for p in grp["pathway"]:
            for w, rel in words[p].items():
                inside[w.lower()] += 1
                display.setdefault(w.lower(), w)
                pos[w.lower()].append(rel)
        need = max(2, 0.08 * len(grp))
        score = {w: (n / len(grp)) * np.log(N / df_all[w]) for w, n in inside.items() if n >= need}
        if not score:  # nothing shared: name the region after its largest gene set
            biggest = grp.sort_values("set_size", ascending=False)["pathway"].iloc[0]
            labels[int(c)] = clean_pathway_name(biggest)
            continue
        top = sorted(score, key=lambda w: -score[w])[:max_words]
        top.sort(key=lambda w: np.mean(pos[w]))  # keep the order words have in the names
        label = " ".join(display[w] for w in top)
        labels[int(c)] = label[:1].upper() + label[1:]
    return labels


def save_map(pmap: pd.DataFrame, path) -> str:
    out = pmap.copy()
    out.insert(0, "method", pmap.attrs.get("method", ""))
    out.to_csv(path, index=False)
    return str(path)


def load_map(path) -> pd.DataFrame:
    df = pd.read_csv(path)
    method = str(df["method"].iloc[0]) if "method" in df.columns and len(df) else ""
    df = df[[c for c in MAP_COLUMNS if c in df.columns]].copy()
    df["cluster_label"] = df["cluster_label"].fillna("")
    df.attrs["method"] = method
    return df


# ---------------------------------------------------------------------------
# Overlaying results
# ---------------------------------------------------------------------------


def map_hits(pmap: pd.DataFrame, table: pd.DataFrame, fdr_cutoff: Optional[float] = 0.05) -> pd.DataFrame:
    """Enriched pathways of an enrichment/comparison table joined to their map position.

    Returns one row per (pathway, condition) with ``x, y, cluster, cluster_label,
    fdr, direction, condition``. Pathways absent from the map (outside its size
    limits or not in the library) are dropped.
    """
    from .compare import condition_names

    df = table.copy()
    if "condition" not in df.columns:
        df["condition"] = ""
    if fdr_cutoff is not None:
        df = df[df["fdr"].notna() & (df["fdr"] <= fdr_cutoff)]
    df = df.sort_values("fdr").drop_duplicates(["pathway", "condition"])
    hits = df.merge(pmap, on="pathway", how="inner")
    cols = ["pathway", "condition", "fdr", "direction", "x", "y", "cluster", "cluster_label"]
    cols += [c for c in ("nes", "overlap", "genes") if c in hits.columns]
    hits = hits[cols].reset_index(drop=True)
    hits.attrs["conditions"] = condition_names(table) or [""]
    return hits


def region_summary(pmap: pd.DataFrame, hits: pd.DataFrame) -> pd.DataFrame:
    """How many enriched pathways fall into each map region, per condition.

    Columns: ``cluster, cluster_label, n_pathways_in_region`` and one count column
    per condition, plus ``fraction_<condition>`` (hits / region size).
    """
    conds = hits.attrs.get("conditions") or sorted(hits["condition"].unique())
    size = pmap[pmap["cluster"] > 0].groupby("cluster").size()
    labels = pmap.groupby("cluster")["cluster_label"].first()
    rows = []
    for c in sorted(size.index):
        row = {"cluster": c, "cluster_label": labels[c], "n_pathways_in_region": int(size[c])}
        sub = hits[hits["cluster"] == c]
        for cond in conds:
            n = int((sub["condition"].astype(str) == str(cond)).sum())
            key = cond or "hits"
            row[f"n_{key}"] = n
            row[f"fraction_{key}"] = round(n / size[c], 3)
        rows.append(row)
    out = pd.DataFrame(rows)
    count_cols = [c for c in out.columns if c.startswith("n_") and c != "n_pathways_in_region"]
    return out.assign(_t=out[count_cols].sum(axis=1)).sort_values("_t", ascending=False).drop(columns="_t").reset_index(drop=True)


def plot_pathway_map(
    pmap: pd.DataFrame,
    table: Optional[pd.DataFrame] = None,
    path=None,
    fdr_cutoff: Optional[float] = 0.05,
    label_top: int = 12,
    show_region_labels: bool = True,
    title: Optional[str] = None,
    width: float = 11.0,
    dpi: int = 300,
    ax=None,
    color_by_category: bool = False,
    categories: Optional[Dict[str, str]] = None,
):
    """Draw the library map, with the enriched pathways of ``table`` highlighted.

    Every pathway is a light grey point; region labels sit at region centres.
    ``color_by_category=True`` colors the points by source collection instead
    (:func:`~simpnetworkviz.network.pathway_category`, or ``categories`` =
    ``{pathway: category}``), for a map built from several databases.
    Enriched pathways are drawn on top: one condition = filled circles sized by
    -log10 FDR; several conditions = small pie charts, one slice per condition,
    as in the comparison networks. The ``label_top`` most significant hits are
    named.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Wedge

    from .plotting import LEGEND_TEXT_SIZE, LEGEND_TITLE_SIZE
    from .style import SEQ_CMAP, condition_color, darken

    from matplotlib import colormaps

    hits = map_hits(pmap, table, fdr_cutoff) if table is not None else pmap.iloc[0:0]
    conds = list(hits.attrs.get("conditions", [""])) if table is not None else []
    multi = len(conds) > 1

    span = max(pmap["x"].max() - pmap["x"].min(), pmap["y"].max() - pmap["y"].min())
    if ax is None:
        fig, ax = plt.subplots(figsize=(width, width * 0.85), facecolor="white")
    else:
        fig = ax.figure
    small = width < 7  # panel of a grid: smaller markers and text
    cat_handles = []
    if color_by_category:
        from .network import pathway_category
        from .style import lighten

        cats = pmap["pathway"].map(lambda n: pathway_category(n, categories))
        order = list(cats.value_counts().index)
        for i, c in enumerate(order):
            m = (cats == c).values
            col = lighten(condition_color(i), 0.35)
            ax.scatter(pmap["x"][m], pmap["y"][m], s=5 if small else 10, c=col, edgecolors="none", zorder=1, alpha=0.9)
            cat_handles.append(plt.Line2D([], [], marker="o", linestyle="", markersize=6, markerfacecolor=col,
                                          markeredgecolor="none", label=f"{c} ({int(m.sum())})"))
    else:
        ax.scatter(pmap["x"], pmap["y"], s=4 if small else 9, c="#D9D9D9", edgecolors="none", zorder=1)
    region_texts = []
    if show_region_labels:
        import textwrap

        for c, grp in pmap[pmap["cluster"] > 0].groupby("cluster"):
            lab = grp["cluster_label"].iloc[0]
            if lab:
                t = ax.text(grp["x"].median(), grp["y"].median(), textwrap.fill(lab, 16 if small else 20),
                            ha="center", va="center", fontsize=6.5 if small else 9, color="#3A3A3A", zorder=5, linespacing=1.05,
                            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.72))
                region_texts.append(t)
    if len(hits):
        sig = -np.log10(hits["fdr"].clip(lower=1e-300))
        vmax = max(float(sig.max()), 1.0)
        r0 = span * 0.006
        by_pw = defaultdict(dict)
        for _, h in hits.iterrows():
            by_pw[h["pathway"]][str(h["condition"])] = (float(h["x"]), float(h["y"]), float(-np.log10(max(h["fdr"], 1e-300))))
        for pw, per in by_pw.items():
            x, y = next(iter(per.values()))[:2]
            best = max(v[2] for v in per.values())
            r = r0 * (1.0 + 1.6 * best / vmax)
            if multi:
                k = len(conds)
                for i, c in enumerate(conds):
                    col = condition_color(i) if c in per else "#FFFFFF"
                    ax.add_patch(Wedge((x, y), r, 90 - 360 * (i + 1) / k, 90 - 360 * i / k, facecolor=col,
                                       edgecolor="#333333", linewidth=0.4, zorder=4))
            else:
                col = colormaps[SEQ_CMAP](0.3 + 0.65 * best / vmax)
                ax.add_patch(Circle((x, y), r, facecolor=col, edgecolor=darken("#C94A2A", 0.3), linewidth=0.5, zorder=4))
        top = hits.sort_values("fdr").drop_duplicates("pathway").head(label_top)
        hit_texts = []
        import textwrap

        offsets = [(10, 10), (-10, 10), (10, -10), (-10, -10)]  # spread the starting positions around the points
        for i, (_, h) in enumerate(top.iterrows()):
            dx, dy = offsets[i % 4]
            a = ax.annotate(textwrap.fill(clean_pathway_name(h["pathway"]), 28), (h["x"], h["y"]), xytext=(dx, dy),
                            textcoords="offset points", ha="left" if dx > 0 else "right", va="bottom" if dy > 0 else "top",
                            fontsize=6.5 if small else 8, color="#1A1A1A", zorder=7, linespacing=1.05,
                            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#BBBBBB", lw=0.4, alpha=0.9),
                            arrowprops=dict(arrowstyle="-", color="#777777", lw=0.5, shrinkA=0, shrinkB=2))
            hit_texts.append(a)
    else:
        hit_texts = []
    ax.set_aspect("equal")
    ax.axis("off")
    _repel_labels(ax, region_texts + hit_texts)
    method = pmap.attrs.get("method", "")
    sub = f"{len(pmap)} gene sets" + (f", {method.upper()} on gene-set similarity" if method else "")
    if title:
        ax.set_title(title, fontsize=9 if small else 14, pad=6 if small else 12)
    ax.text(0.0, -0.02, sub, transform=ax.transAxes, fontsize=LEGEND_TEXT_SIZE - (3 if small else 1), color="#666666", va="top")
    if cat_handles:
        # A figure legend anchored to the axes: it survives the condition legend below and is
        # included in bbox_inches="tight" (an axes legend re-added with add_artist is not).
        fig.legend(handles=cat_handles, title="Collection", loc="upper left", bbox_transform=ax.transAxes,
                   bbox_to_anchor=(1.0, 0.55 if (len(hits) and multi) else 1.0),
                   frameon=False, fontsize=LEGEND_TEXT_SIZE - 1, title_fontsize=LEGEND_TITLE_SIZE - 2)
    if len(hits):
        if multi:
            handles = [plt.Line2D([], [], marker="o", linestyle="", markersize=7, markerfacecolor=condition_color(i),
                                  markeredgecolor="#333333", label=f"{i + 1}  {c}") for i, c in enumerate(conds)]
            ax.legend(handles=handles, title="Enriched in", loc="upper left",
                      bbox_to_anchor=(1.0, 1.0), frameon=False, fontsize=LEGEND_TEXT_SIZE - 1, title_fontsize=LEGEND_TITLE_SIZE - 2)
        else:
            ax.text(1.0, -0.02, f"colored = enriched (FDR ≤ {fdr_cutoff}); size and shade = -log10 FDR (max {vmax:.1f})",
                    transform=ax.transAxes, fontsize=LEGEND_TEXT_SIZE - 1, color="#666666", va="top", ha="right")
    if path:
        fig.savefig(path, dpi=dpi, facecolor="white", bbox_inches="tight")
    return fig
