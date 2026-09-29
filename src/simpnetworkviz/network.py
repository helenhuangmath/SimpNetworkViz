"""Build pathway (EnrichmentMap-style) and gene networks."""
from __future__ import annotations

import re
import warnings
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import networkx as nx
import numpy as np
import pandas as pd

from .enrichment import read_enrichment, split_genes
from .genesets import GeneSets, normalize_gene

# ---------------------------------------------------------------------------
# Network container
# ---------------------------------------------------------------------------


@dataclass
class Network:
    """A laid-out network ready for plotting/export.

    Attributes
    ----------
    graph : networkx graph; node attributes include ``label``, ``type``
        (``pathway``/``gene``), ``size_value``, ``cluster`` and, where
        available, ``fdr``, ``direction``, ``nes``, ``logfc``, ``genes``.
    kind : ``"pathway"``, ``"gene"`` or ``"gene-pathway"``
    groups : node -> cluster id (0 = not in a cluster)
    group_labels : cluster id -> annotation text
    pos : node -> (x, y) in layout units (≈ inches at scale 1)
    """

    graph: nx.Graph
    kind: str
    groups: Dict[str, int]
    group_labels: Dict[int, str]
    pos: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    radius: Dict[str, float] = field(default_factory=dict)
    layout_style: str = "clustered"
    label_nodes: set = field(default_factory=set)

    # -- tables -------------------------------------------------------------
    def node_table(self) -> pd.DataFrame:
        rows = []
        for n, d in self.graph.nodes(data=True):
            row = {"id": n}
            row.update({k: v for k, v in d.items() if not isinstance(v, (list, set, tuple))})
            row["cluster"] = self.groups.get(n, 0)
            row["cluster_label"] = self.group_labels.get(row["cluster"], "")
            if n in self.pos:
                row["x"], row["y"] = self.pos[n]
            rows.append(row)
        df = pd.DataFrame(rows)
        order = df["cluster"].replace(0, 10**9)
        return df.assign(_o=order).sort_values(["_o", "id"]).drop(columns="_o").reset_index(drop=True)

    def edge_table(self) -> pd.DataFrame:
        rows = [{"source": u, "target": v, **d} for u, v, d in self.graph.edges(data=True)]
        return pd.DataFrame(rows, columns=["source", "target", "weight"] if not rows else None)

    def cluster_table(self) -> pd.DataFrame:
        members = defaultdict(list)
        for n, g in self.groups.items():
            if g:
                members[g].append(n)
        rows = []
        for g in sorted(members):
            nodes = members[g]
            fdrs = [self.graph.nodes[n].get("fdr") for n in nodes]
            fdrs = [f for f in fdrs if f is not None and not pd.isna(f)]
            rows.append({
                "cluster": g,
                "label": self.group_labels.get(g, ""),
                "n_nodes": len(nodes),
                "best_fdr": min(fdrs) if fdrs else np.nan,
                "members": "; ".join(self.graph.nodes[n].get("name", n) for n in nodes),
            })
        return pd.DataFrame(rows)

    # -- convenience wrappers --------------------------------------------------
    def plot(self, path=None, **kwargs):
        from .plotting import plot_network
        return plot_network(self, path, **kwargs)

    def to_html(self, path, **kwargs):
        from .export import write_html
        return write_html(self, path, **kwargs)

    def to_cytoscape(self, path, **kwargs):
        from .export import write_cyjs
        return write_cyjs(self, path, **kwargs)

    def to_graphml(self, path):
        from .export import write_graphml
        return write_graphml(self, path)

    def save(self, outdir, prefix: str = "network", formats: Sequence[str] = ("png", "pdf", "html", "cyjs", "graphml", "csv"), **plot_kwargs):
        from .export import save_all
        return save_all(self, outdir, prefix, formats, **plot_kwargs)


# ---------------------------------------------------------------------------
# Pathway names and cluster annotation
# ---------------------------------------------------------------------------

_PREFIX = re.compile(
    r"^(HALLMARK|KEGG|REACTOME|GOBP|GOCC|GOMF|GO|WP|BIOCARTA|PID|NABA|WIKIPATHWAYS|BIOPLANET|KEGG_MEDICUS(_REFERENCE)?)[_:\s]+",
    re.I,
)
_SUFFIXES = [
    r"\s*\((GO|R-HSA|R-MMU|WP|hsa|mmu)[:\-_]?[\w\-.]*\)\s*$",
    r"\s+R-(HSA|MMU)-\d+\s*$",
    r"\s+WP\d+\s*$",
    r"\s+hsa\d{5}\s*$",
    r"\s*-\s*(Homo sapiens|Mus musculus)\s*\(\w+\)\s*$",
    r"\s+GO:\d{7}\s*$",
]
_KEEP_UPPER = {
    "DNA", "RNA", "MRNA", "TRNA", "RRNA", "SNRNA", "ATP", "GTP", "NADH", "TNF", "TNFA", "NFKB", "NF-KB", "JAK", "STAT",
    "MTOR", "MTORC1", "MTORC2", "PI3K", "AKT", "MAPK", "ERK", "JNK", "TGF", "TGFB", "TGF-BETA", "WNT", "ECM", "UV", "UPR",
    "KRAS", "MYC", "HIF", "HIF1", "NK", "MHC", "HLA", "TCR", "BCR", "ER", "TCA", "RHO", "RAS", "EGFR", "VEGF", "IL", "P53",
    "G1", "G2", "G2M", "E2F", "II", "III", "I", "NOTCH", "HEDGEHOG", "PD-1", "PD-L1", "CD4", "CD8", "SUMO", "HDR", "NHEJ",
}
STOPWORDS = {
    "of", "the", "and", "in", "to", "by", "via", "for", "from", "on", "with", "a", "an", "as", "at", "or", "into",
    "involved", "process", "processes", "pathway", "pathways", "positive", "negative", "regulation", "regulated",
    "cellular", "hallmark", "homo", "sapiens", "human", "mus", "musculus", "mouse", "biological", "system", "activity",
    "type", "other", "mediated", "dependent", "containing", "complex", "through", "during", "cell", "cells",
}


def clean_pathway_name(name: str) -> str:
    """Make database pathway names human readable.

    ``HALLMARK_INTERFERON_GAMMA_RESPONSE`` -> ``Interferon gamma response``;
    ``cell cycle (GO:0007049)`` -> ``Cell cycle``.
    """
    s = str(name).strip()
    for suf in _SUFFIXES:
        s = re.sub(suf, "", s, flags=re.I)
    s = _PREFIX.sub("", s)
    s = s.replace("_", " ").strip()
    s = re.sub(r"\s+", " ", s)
    letters = [c for c in s if c.isalpha()]
    if letters and sum(c.isupper() for c in letters) / len(letters) > 0.8:
        words = []
        for w in s.split(" "):
            core = w.strip("(),")
            if core.upper() in _KEEP_UPPER or any(ch.isdigit() for ch in core):
                words.append(w)
            else:
                words.append(w.lower())
        s = " ".join(words)
    return s[:1].upper() + s[1:] if s else s


def _tokens(name: str) -> List[str]:
    return [t for t in re.split(r"[\s,/()\[\]:;]+", name) if t]


def annotate_cluster(names: Sequence[str], weights: Optional[Sequence[float]] = None, max_words: int = 3,
                     method: str = "auto") -> str:
    """Summarize a cluster of pathway names into a short label.

    ``method="words"`` picks the most frequent informative words (similar to
    Cytoscape AutoAnnotate's WordCloud), ordered as they typically appear in the
    names. ``"representative"`` returns the most significant pathway name.
    ``"auto"`` uses words for clusters of 3+ pathways when they are shared by at
    least two pathways, otherwise the representative name.
    """
    clean = [clean_pathway_name(n) for n in names]
    if weights is None:
        weights = [1.0] * len(clean)
    rep = clean[int(np.argmax(weights))] if clean else ""
    if method == "representative" or not clean:
        return rep
    if method == "auto" and len(clean) < 3:
        return rep
    doc_freq: Counter = Counter()
    score: Counter = Counter()
    positions = defaultdict(list)
    display = {}
    for name, w in zip(clean, weights):
        toks = _tokens(name)
        seen = set()
        for i, t in enumerate(toks):
            key = t.lower().strip("-")
            if key in STOPWORDS or len(key) < 2 or key.isdigit():
                continue
            if key not in seen:
                seen.add(key)
                doc_freq[key] += 1
                score[key] += 1 + 0.1 * w
                positions[key].append(i / max(len(toks) - 1, 1))
                display.setdefault(key, t if t.upper() in _KEEP_UPPER or any(c.isdigit() for c in t) else key)
    if not score:
        return rep
    ranked = sorted(score, key=lambda k: (-doc_freq[k], -score[k]))
    chosen = [k for k in ranked[:max_words] if doc_freq[k] >= 2]
    if not chosen:
        return rep
    chosen.sort(key=lambda k: np.mean(positions[k]))
    label = " ".join(display[k] for k in chosen)
    return label[:1].upper() + label[1:]


# ---------------------------------------------------------------------------
# Similarity + clustering
# ---------------------------------------------------------------------------


def similarity(a: set, b: set, metric: str = "combined", combined_k: float = 0.5) -> float:
    if not a or not b:
        return 0.0
    inter = len(a & b)
    if inter == 0:
        return 0.0
    jac = inter / len(a | b)
    ovl = inter / min(len(a), len(b))
    if metric == "jaccard":
        return jac
    if metric == "overlap":
        return ovl
    if metric == "combined":
        return combined_k * jac + (1 - combined_k) * ovl
    raise ValueError(f"Unknown similarity metric {metric!r} (use jaccard, overlap or combined)")


def mcl_clusters(G: nx.Graph, inflation: float = 2.0, weight: str = "weight", max_iter: int = 100) -> List[set]:
    """Markov clustering (MCL), as used by Cytoscape AutoAnnotate/clusterMaker."""
    nodes = list(G.nodes)
    if not nodes:
        return []
    A = nx.to_numpy_array(G, nodelist=nodes, weight=weight)
    # Self-loops weighted like the strongest edge (standard MCL practice) so
    # that loops do not overwhelm weak similarity edges.
    loops = A.max(axis=1)
    A = A + np.diag(np.where(loops > 0, loops, 1.0))
    M = A / A.sum(axis=0, keepdims=True)
    for _ in range(max_iter):
        prev = M
        M = M @ M
        M = np.power(M, inflation)
        M[M < 1e-8] = 0
        M = M / M.sum(axis=0, keepdims=True)
        if np.allclose(M, prev, atol=1e-8):
            break
    clusters, assigned = [], set()
    for i in range(len(nodes)):
        members = {nodes[j] for j in np.nonzero(M[i] > 1e-6)[0]} - assigned
        if M[i].sum() > 1e-6 and members:
            clusters.append(members)
            assigned |= members
    for n in nodes:
        if n not in assigned:
            clusters.append({n})
    return clusters


def cluster_graph(G: nx.Graph, method: str = "mcl", inflation: float = 2.0, resolution: float = 1.0, seed: int = 0) -> List[set]:
    if method == "mcl":
        return mcl_clusters(G, inflation=inflation)
    if method == "louvain":
        return [set(c) for c in nx.community.louvain_communities(G, weight="weight", resolution=resolution, seed=seed)]
    if method in ("components", "connected"):
        return [set(c) for c in nx.connected_components(G)]
    if method in ("none", None):
        return [{n} for n in G.nodes]
    raise ValueError(f"Unknown cluster method {method!r} (use mcl, louvain, components or none)")


def _rank_groups(clusters: List[set], score, min_cluster_size: int) -> Tuple[Dict[str, int], List[List[str]]]:
    """Order clusters by size then significance; clusters < min size get id 0."""
    clusters = sorted(clusters, key=lambda c: (-len(c), -max(score(n) for n in c), sorted(c)[0]))
    groups, ordered = {}, []
    gid = 0
    for c in clusters:
        members = sorted(c, key=lambda n: -score(n))
        if len(c) >= min_cluster_size:
            gid += 1
            ordered.append(members)
            for n in members:
                groups[n] = gid
        else:
            for n in members:
                groups[n] = 0
    return groups, ordered


def _neg_log10(x) -> float:
    if x is None or pd.isna(x):
        return 0.0
    return float(-np.log10(max(float(x), 1e-300)))


# ---------------------------------------------------------------------------
# Pathway network
# ---------------------------------------------------------------------------


def build_pathway_network(
    table,
    gene_sets: Optional[GeneSets] = None,
    fdr_cutoff: Optional[float] = 0.05,
    max_pathways: Optional[int] = 100,
    similarity_metric: str = "combined",
    similarity_cutoff: float = 0.375,
    use_full_gene_sets: bool = True,
    cluster_method: str = "mcl",
    mcl_inflation: float = 2.0,
    min_cluster_size: int = 2,
    label_method: str = "auto",
    label_words: int = 3,
    cluster_labels: Optional[Dict[int, str]] = None,
    layout: str = "clustered",
    node_labels="auto",
    seed: int = 0,
) -> Network:
    """Build an EnrichmentMap-style pathway network.

    Nodes are pathways (size = gene-set size, attributes carry FDR/direction),
    edges connect pathways whose gene sets overlap (similarity >= cutoff), and
    pathways are clustered (MCL by default) and annotated with a short label.

    ``table`` is an enrichment table (DataFrame or path; see
    :func:`read_enrichment`). ``gene_sets`` supplies the full pathway gene
    lists used for similarity (as EnrichmentMap does); without it, the genes in
    the table are used.
    """
    df = table if isinstance(table, pd.DataFrame) and "genes" in table.columns and "fdr" in table.columns else read_enrichment(table)
    df = df.copy()
    if fdr_cutoff is not None and df["fdr"].notna().any():
        df = df[df["fdr"].isna() | (df["fdr"] <= fdr_cutoff)]
    df = df.sort_values("fdr", na_position="last")
    if max_pathways:
        df = df.head(max_pathways)
    if df.empty:
        raise ValueError("No pathways left to plot. Relax the FDR cutoff (e.g. --fdr 0.25) or check your input.")

    G = nx.Graph()
    sets: Dict[str, set] = {}
    dup = df["pathway"].duplicated(keep=False)
    for (_, row), is_dup in zip(df.iterrows(), dup):
        direction = str(row.get("direction", "") or "")
        node = f"{row['pathway']} [{direction}]" if is_dup and direction else str(row["pathway"])
        if node in G:
            continue
        hits = split_genes(row.get("genes", ""))
        full = gene_sets.get(row["pathway"]) if gene_sets else None
        members = set(full) if (use_full_gene_sets and full) else set(hits)
        if not members and full:
            members = set(full)
        sets[node] = members
        set_size = row.get("set_size")
        if set_size is None or pd.isna(set_size):
            set_size = len(full) if full else len(members)
        G.add_node(
            node,
            name=str(row["pathway"]),
            label=clean_pathway_name(row["pathway"]),
            type="pathway",
            fdr=None if pd.isna(row.get("fdr", np.nan)) else float(row["fdr"]),
            pvalue=None if pd.isna(row.get("pvalue", np.nan)) else float(row["pvalue"]),
            nes=None if pd.isna(row.get("nes", np.nan)) else float(row["nes"]),
            direction=direction,
            set_size=int(set_size),
            overlap=len(hits),
            genes=";".join(hits),
            size_value=float(set_size),
        )
    nodes = list(G.nodes)
    for i, a in enumerate(nodes):
        for b in nodes[i + 1:]:
            s = similarity(sets[a], sets[b], similarity_metric)
            if s >= similarity_cutoff and s > 0:
                G.add_edge(a, b, weight=round(s, 4), shared=len(sets[a] & sets[b]))

    score = lambda n: _neg_log10(G.nodes[n].get("fdr")) + 1e-6 * G.nodes[n]["size_value"]
    clusters = cluster_graph(G, cluster_method, inflation=mcl_inflation, seed=seed)
    groups, ordered = _rank_groups(clusters, score, min_cluster_size)
    labels = {}
    for gid, members in enumerate(ordered, start=1):
        labels[gid] = annotate_cluster(
            [G.nodes[n]["name"] for n in members], [score(n) for n in members], label_words, label_method
        )
    if cluster_labels:
        labels.update({int(k): v for k, v in cluster_labels.items()})
    for n in G.nodes:
        G.nodes[n]["cluster"] = groups[n]
        G.nodes[n]["cluster_label"] = labels.get(groups[n], "")

    net = Network(G, "pathway", groups, labels)
    _layout(net, layout, seed, node_labels)
    return net


# ---------------------------------------------------------------------------
# Gene networks
# ---------------------------------------------------------------------------


def _pathway_genes(pnet: Network, gene_sets: Optional[GeneSets], nodes: Iterable[str]) -> Dict[str, List[str]]:
    out = {}
    for n in nodes:
        d = pnet.graph.nodes[n]
        genes = [g for g in d.get("genes", "").split(";") if g]
        if not genes and gene_sets:
            genes = list(gene_sets.get(d["name"], []))
        out[n] = genes
    return out


def _select_pathways(pnet: Network, clusters: Optional[Iterable[int]], pathways: Optional[Iterable[str]], max_pathways: Optional[int]) -> List[str]:
    G = pnet.graph
    if pathways:
        wanted = set(pathways)
        nodes = [n for n in G if n in wanted or G.nodes[n]["name"] in wanted]
    else:
        nodes = list(G)
        if clusters:
            keep = {int(c) for c in clusters}
            nodes = [n for n in nodes if pnet.groups.get(n, 0) in keep]
    nodes.sort(key=lambda n: -_neg_log10(G.nodes[n].get("fdr")))
    if max_pathways and len(nodes) > max_pathways:
        # Round-robin over clusters so every theme stays represented.
        by_group = defaultdict(list)
        for n in nodes:
            by_group[pnet.groups.get(n, 0)].append(n)
        order = sorted(by_group, key=lambda g: (g == 0, g))
        picked = []
        while len(picked) < max_pathways and any(by_group.values()):
            for g in order:
                if by_group[g] and len(picked) < max_pathways:
                    picked.append(by_group[g].pop(0))
        nodes = picked
    return nodes


def fetch_string_edges(genes: Sequence[str], species: int = 9606, required_score: int = 400, timeout: float = 60) -> List[Tuple[str, str, float]]:
    """Protein-protein interactions among ``genes`` from the STRING API."""
    import requests

    resp = requests.post(
        "https://string-db.org/api/tsv/network",
        data={"identifiers": "\r".join(genes), "species": species, "required_score": required_score,
              "caller_identity": "simpnetworkviz"},
        timeout=timeout,
    )
    resp.raise_for_status()
    lines = resp.text.strip().splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    ia, ib, isc = header.index("preferredName_A"), header.index("preferredName_B"), header.index("score")
    wanted = {normalize_gene(g) for g in genes}
    edges = {}
    for line in lines[1:]:
        f = line.split("\t")
        a, b = normalize_gene(f[ia]), normalize_gene(f[ib])
        if a in wanted and b in wanted and a != b:
            key = tuple(sorted((a, b)))
            edges[key] = max(edges.get(key, 0), float(f[isc]))
    return [(a, b, s) for (a, b), s in edges.items()]


def build_gene_network(
    pathway_network: Network,
    gene_sets: Optional[GeneSets] = None,
    mode: str = "shared",
    clusters: Optional[Iterable[int]] = None,
    pathways: Optional[Iterable[str]] = None,
    gene_values: Optional[Dict[str, Optional[float]]] = None,
    max_pathways: Optional[int] = None,
    max_genes: int = 150,
    min_shared: int = 1,
    string_species: int = 9606,
    string_score: int = 400,
    layout: Optional[str] = None,
    node_labels="auto",
    seed: int = 0,
) -> Network:
    """Build a gene-level network from (a subset of) a pathway network.

    mode
        ``"shared"``  genes are linked when they belong to the same pathway(s);
                      edge weight = number of shared pathways (offline).
        ``"string"``  protein-protein interactions from STRING (needs internet).
        ``"cnet"``    bipartite gene-pathway "concept network" (like
                      clusterProfiler's cnetplot).

    Genes are grouped by the pathway cluster they mostly belong to, so the gene
    modules carry the same annotation as the pathway network. Genes come from
    the pathway "genes" column (your query genes); for pathway lists without
    query genes the full gene sets are used. ``gene_values`` (e.g. log2FC)
    colors the gene nodes.
    """
    if mode not in {"shared", "string", "cnet"}:
        raise ValueError("mode must be 'shared', 'string' or 'cnet'")
    if max_pathways is None:
        max_pathways = 15 if mode == "cnet" else None
    pnodes = _select_pathways(pathway_network, clusters, pathways, max_pathways)
    if not pnodes:
        raise ValueError("No pathways selected for the gene network.")
    pgenes = _pathway_genes(pathway_network, gene_sets, pnodes)
    gene_values = {normalize_gene(k): v for k, v in (gene_values or {}).items()}

    membership: Dict[str, List[str]] = defaultdict(list)
    for p in pnodes:
        for g in pgenes[p]:
            membership[g].append(p)
    if not membership:
        raise ValueError("Selected pathways have no genes. Provide gene sets (--gmt/--library) or a gene column.")

    def gene_rank(g):
        v = gene_values.get(g)
        return (len(membership[g]), abs(v) if v is not None and not pd.isna(v) else 0.0)

    genes = sorted(membership, key=gene_rank, reverse=True)
    if max_genes and len(genes) > max_genes:
        genes = genes[:max_genes]
    gene_group = {}
    for g in genes:
        counts = Counter(pathway_network.groups.get(p, 0) for p in membership[g])
        best = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0] == 0, kv[0]))[0][0]
        gene_group[g] = best

    G = nx.Graph()
    for g in genes:
        v = gene_values.get(g)
        G.add_node(
            g, name=g, label=g, type="gene",
            logfc=None if v is None or pd.isna(v) else float(v),
            n_pathways=len(membership[g]),
            pathways="; ".join(pathway_network.graph.nodes[p]["label"] for p in membership[g]),
            size_value=float(len(membership[g])),
        )
    groups = dict(gene_group)
    labels = dict(pathway_network.group_labels)

    if mode == "cnet":
        for p in pnodes:
            d = pathway_network.graph.nodes[p]
            G.add_node(p, **{k: v for k, v in d.items()})
            G.nodes[p]["size_value"] = float(max(len([g for g in pgenes[p] if g in G]), 1))
            groups[p] = pathway_network.groups.get(p, 0)
            for g in pgenes[p]:
                if g in gene_group:
                    G.add_edge(p, g, weight=1.0)
        kind = "gene-pathway"
        default_layout = "spring"
    elif mode == "shared":
        gl = list(genes)
        sets = {g: set(membership[g]) for g in gl}
        for i, a in enumerate(gl):
            for b in gl[i + 1:]:
                k = len(sets[a] & sets[b])
                if k >= min_shared:
                    G.add_edge(a, b, weight=float(k), shared_pathways=k)
        kind = "gene"
        default_layout = "clustered"
    else:
        try:
            edges = fetch_string_edges(genes, species=string_species, required_score=string_score)
        except Exception as exc:  # network errors, API changes
            raise RuntimeError(f"Could not query STRING ({exc}). Use mode='shared' for an offline gene network.") from exc
        for a, b, s in edges:
            G.add_edge(a, b, weight=s, string_score=s)
        kind = "gene"
        default_layout = "clustered"

    for n in G.nodes:
        G.nodes[n]["cluster"] = groups.get(n, 0)
        G.nodes[n]["cluster_label"] = labels.get(groups.get(n, 0), "")
    used = {groups[n] for n in G.nodes}
    labels = {k: v for k, v in labels.items() if k in used}
    net = Network(G, kind, groups, labels)
    _layout(net, layout or default_layout, seed, node_labels)
    return net


def _layout(net: Network, style: str, seed: int, node_labels="auto") -> None:
    from .layout import clustered_layout, node_radii, spring_layout
    from .plotting import resolve_node_labels

    net.radius = node_radii(net)
    net.label_nodes = resolve_node_labels(net, node_labels)
    if style == "clustered":
        net.pos = clustered_layout(net, seed=seed)
    elif style == "spring":
        net.pos = spring_layout(net, seed=seed)
    else:
        raise ValueError("layout must be 'clustered' or 'spring'")
    net.layout_style = style
