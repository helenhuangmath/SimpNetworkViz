"""Over-representation analysis and reading enrichment tables from other tools."""
from __future__ import annotations

import re
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from .genesets import GeneSets, normalize_gene

#: Standard columns of an enrichment table used throughout the package.
STANDARD_COLUMNS = [
    "pathway", "overlap", "set_size", "query_size", "background_size",
    "fold_enrichment", "pvalue", "fdr", "direction", "genes",
]


def bh_fdr(pvalues) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def enrich(
    genes: Iterable[str],
    gene_sets: GeneSets,
    background: Optional[Iterable[str]] = None,
    min_size: int = 5,
    max_size: int = 1000,
    min_overlap: int = 2,
    direction: Optional[str] = None,
) -> pd.DataFrame:
    """Hypergeometric over-representation analysis with BH correction.

    Parameters
    ----------
    genes : query genes (e.g. differentially expressed genes)
    gene_sets : ``{pathway: [genes]}``, see :func:`load_gene_sets`
    background : universe of genes; defaults to all genes in ``gene_sets``.
        Using the genes actually measured in your experiment is recommended.
    min_size, max_size : gene-set size limits (after restricting to background)
    min_overlap : pathways with fewer query genes are not reported
    direction : optional label stored in the ``direction`` column ("up"/"down")
    """
    if background is None:
        universe = set().union(*[set(v) for v in gene_sets.values()]) if gene_sets else set()
    else:
        universe = {normalize_gene(g) for g in background}
    query = {normalize_gene(g) for g in genes} & universe
    N, n = len(universe), len(query)
    rows = []
    for name, members in gene_sets.items():
        members_u = set(members) & universe
        K = len(members_u)
        if K < min_size or K > max_size:
            continue
        hits = sorted(query & members_u)
        k = len(hits)
        if k < max(min_overlap, 1):
            continue
        p = float(hypergeom.sf(k - 1, N, K, n))
        rows.append({
            "pathway": name,
            "overlap": k,
            "set_size": K,
            "query_size": n,
            "background_size": N,
            "fold_enrichment": (k / n) / (K / N) if n and K else np.nan,
            "pvalue": p,
            "direction": direction or "",
            "genes": ";".join(hits),
        })
    df = pd.DataFrame(rows, columns=[c for c in STANDARD_COLUMNS if c != "fdr"])
    # Correct for every gene set tested, not only those that passed min_overlap.
    n_tested = sum(1 for v in gene_sets.values() if min_size <= len(set(v) & universe) <= max_size)
    if len(df):
        padded = np.concatenate([df["pvalue"].values, np.ones(max(n_tested - len(df), 0))])
        df["fdr"] = bh_fdr(padded)[: len(df)]
    else:
        df["fdr"] = []
    return df[STANDARD_COLUMNS].sort_values(["fdr", "pvalue"]).reset_index(drop=True)


def enrich_updown(up: Iterable[str], down: Iterable[str], gene_sets: GeneSets, **kwargs) -> pd.DataFrame:
    """Run :func:`enrich` separately for up- and down-regulated genes and stack the results."""
    res = [enrich(up, gene_sets, direction="up", **kwargs), enrich(down, gene_sets, direction="down", **kwargs)]
    res = [r for r in res if len(r)]
    if not res:
        return enrich([], gene_sets, **kwargs)
    return pd.concat(res).sort_values(["fdr", "pvalue"]).reset_index(drop=True)


def split_by_sign(gene_values: Dict[str, Optional[float]]):
    """Split ``{gene: log2FC}`` into up/down lists. Returns ``None`` if no values."""
    vals = {g: v for g, v in gene_values.items() if v is not None and not np.isnan(v)}
    if not vals:
        return None
    up = [g for g, v in vals.items() if v > 0]
    down = [g for g, v in vals.items() if v < 0]
    return up, down


# ---------------------------------------------------------------------------
# Reading enrichment tables produced by other tools
# ---------------------------------------------------------------------------

_ALIASES = {
    "pathway": ["pathway", "term", "description", "term_name", "name", "pathway_name",
                "gene_set", "geneset", "go_term", "native", "id", "term_id"],
    "fdr": ["fdr", "adjusted p-value", "adjusted_pvalue", "p.adjust", "padj", "qvalue", "q-value",
            "fdr q-val", "fdr_q_val", "benjamini", "p_adj", "adj_p", "adj.p.val", "q.value", "fdr_bh"],
    "pvalue": ["pvalue", "p-value", "p.value", "p_value", "pval", "nom p-val", "p", "raw p-value"],
    "genes": ["genes", "geneid", "gene_id", "intersections", "overlapping genes", "core_enrichment",
              "leading_edge", "leadingedge", "gene_list", "lead_genes", "hits", "members", "genelist"],
    "nes": ["nes", "normalized enrichment score", "es", "enrichmentscore", "zscore", "z-score"],
    "direction": ["direction", "regulation", "sign", "updown", "up_down", "group", "cluster"],
    "set_size": ["set_size", "setsize", "size", "term_size", "gs_size", "list_size", "bgratio"],
    "overlap": ["overlap", "count", "intersection_size", "gene_count", "n_genes", "generatio"],
}


def _find(cols: Dict[str, str], key: str, used: set) -> Optional[str]:
    for alias in _ALIASES[key]:
        if alias in cols and cols[alias] not in used:
            return cols[alias]
    return None


def _ratio_or_int(x):
    """'12/300' -> 12 (numerator); '12' -> 12."""
    if pd.isna(x):
        return np.nan
    s = str(x)
    m = re.match(r"\s*(\d+)\s*/\s*(\d+)", s)
    if m:
        return int(m.group(1))
    try:
        return float(s)
    except ValueError:
        return np.nan


def _ratio_denominator(x):
    if pd.isna(x):
        return np.nan
    m = re.match(r"\s*(\d+)\s*/\s*(\d+)", str(x))
    return int(m.group(2)) if m else np.nan


def split_genes(value) -> List[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, (list, tuple, set)):
        return [normalize_gene(str(g)) for g in value if str(g).strip()]
    s = str(value).strip().strip("[]")
    toks = re.split(r"[;/,|\s]+", s)
    return [normalize_gene(t.strip("'\"")) for t in toks if t.strip("'\"")]


def read_enrichment(source, sep: Optional[str] = None) -> pd.DataFrame:
    """Read an enrichment table and map it to the standard columns.

    Works with tables from this package, Enrichr, clusterProfiler (``enrichGO``,
    ``enrichKEGG``, ``GSEA``), g:Profiler, GSEA (``NES``/``FDR q-val``), Metascape-like
    exports and simple two-column ``pathway, fdr`` tables. ``source`` may be a
    path or a :class:`pandas.DataFrame`.
    """
    if isinstance(source, pd.DataFrame):
        raw = source.copy()
    else:
        path = str(source)
        if sep is None:
            sep = "," if path.lower().endswith(".csv") else None
        if path.lower().endswith((".xlsx", ".xls")):
            raw = pd.read_excel(path)
        else:
            raw = pd.read_csv(path, sep=sep, engine="python")
    cols = {str(c).strip().lower(): c for c in raw.columns}
    used: set = set()
    mapping = {}
    # Order matters: resolve most specific keys first.
    for key in ["genes", "nes", "direction", "set_size", "overlap", "pvalue", "fdr", "pathway"]:
        col = _find(cols, key, used)
        if col is not None:
            mapping[key] = col
            used.add(col)
    # clusterProfiler has both ID and Description: prefer Description as name.
    if "description" in cols and mapping.get("pathway") != cols["description"]:
        mapping["pathway"] = cols["description"]
    if "pathway" not in mapping:
        mapping["pathway"] = raw.columns[0]
    if "fdr" not in mapping and "pvalue" in mapping:
        mapping["fdr"] = mapping["pvalue"]

    out = pd.DataFrame({"pathway": raw[mapping["pathway"]].astype(str).str.strip()})
    out["fdr"] = pd.to_numeric(raw[mapping["fdr"]], errors="coerce") if "fdr" in mapping else np.nan
    out["pvalue"] = pd.to_numeric(raw[mapping["pvalue"]], errors="coerce") if "pvalue" in mapping else out["fdr"]
    out["genes"] = raw[mapping["genes"]].apply(lambda v: ";".join(split_genes(v))) if "genes" in mapping else ""
    if "overlap" in mapping:
        out["overlap"] = raw[mapping["overlap"]].apply(_ratio_or_int)
    else:
        out["overlap"] = out["genes"].apply(lambda s: len(s.split(";")) if s else np.nan)
    if "set_size" in mapping:
        col = raw[mapping["set_size"]]
        out["set_size"] = col.apply(_ratio_or_int)
    elif "overlap" in mapping:
        # Enrichr "Overlap" column is "k/K"
        out["set_size"] = raw[mapping["overlap"]].apply(_ratio_denominator)
    else:
        out["set_size"] = np.nan
    out["nes"] = pd.to_numeric(raw[mapping["nes"]], errors="coerce") if "nes" in mapping else np.nan
    if "direction" in mapping:
        out["direction"] = raw[mapping["direction"]].astype(str).str.lower().map(_norm_direction)
    elif out["nes"].notna().any():
        out["direction"] = np.where(out["nes"] > 0, "up", np.where(out["nes"] < 0, "down", ""))
    else:
        out["direction"] = ""
    out = out[out["pathway"].str.len() > 0]
    return out.sort_values("fdr", na_position="last").reset_index(drop=True)


def _norm_direction(s: str) -> str:
    s = s.strip().lower()
    if s in {"up", "upregulated", "up-regulated", "+", "1", "pos", "positive", "increased", "activated"}:
        return "up"
    if s in {"down", "downregulated", "down-regulated", "-", "-1", "neg", "negative", "decreased", "suppressed"}:
        return "down"
    return s if s not in {"nan", "none"} else ""


def read_pathway_list(path) -> List[str]:
    """Read pathway names, one per line."""
    with open(path, encoding="utf-8") as fh:
        return [l.strip() for l in fh if l.strip() and not l.startswith("#")]


def table_from_pathways(pathways: Iterable[str], gene_sets: GeneSets, query: Optional[Iterable[str]] = None) -> pd.DataFrame:
    """Build a minimal enrichment-like table for a user supplied list of pathways.

    Pathway names are matched exactly first, then case-insensitively ignoring
    underscores/punctuation, then by their cleaned name (``"G2M checkpoint"``
    matches ``HALLMARK_G2M_CHECKPOINT``). ``query`` genes (optional) populate the ``genes`` column.
    """
    def key(s: str) -> str:
        return re.sub(r"[^a-z0-9]+", "", s.lower())

    from .network import clean_pathway_name

    lookup = {key(clean_pathway_name(k)): k for k in gene_sets}
    lookup.update({key(k): k for k in gene_sets})
    qset = {normalize_gene(g) for g in query} if query is not None else None
    rows, missing = [], []
    for p in pathways:
        name = p if p in gene_sets else lookup.get(key(p))
        if name is None:
            missing.append(p)
            continue
        members = gene_sets[name]
        genes = [g for g in members if qset is None or g in qset]
        rows.append({"pathway": name, "fdr": np.nan, "pvalue": np.nan, "set_size": len(members),
                     "overlap": len(genes), "direction": "", "nes": np.nan, "genes": ";".join(genes)})
    if missing:
        import warnings
        warnings.warn(f"{len(missing)} pathway(s) not found in gene sets: {missing[:5]}{'...' if len(missing) > 5 else ''}", stacklevel=2)
    return pd.DataFrame(rows)
