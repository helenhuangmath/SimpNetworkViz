"""Single-cell marker tables -> per-cluster gene lists.

Reads the marker tables of Seurat (``FindAllMarkers``) and scanpy
(``rank_genes_groups``), or any long table with cluster, gene, log fold-change
and adjusted p-value columns, and returns ``{cluster: {gene: log2FC}}`` ready for
:func:`~simpnetworkviz.compare.compare_conditions`.
"""
from __future__ import annotations

import re
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd

from .genesets import normalize_gene

_COLS = {
    "cluster": ["cluster", "group", "groups", "cell_type", "celltype", "cluster_id", "ident", "leiden", "louvain", "comparison"],
    "gene": ["gene", "genes", "names", "name", "symbol", "gene_symbol", "feature", "features", "gene_name"],
    "logfc": ["avg_log2fc", "avg_logfc", "logfoldchanges", "log2fc", "logfc", "log2foldchange", "avg_diff", "lfc"],
    "padj": ["p_val_adj", "pvals_adj", "padj", "p.adjust", "fdr", "adj_p", "qval", "q_value", "p_val_adj_bonf"],
    "pval": ["p_val", "pvals", "pvalue", "p_value", "pval"],
    "score": ["scores", "score", "t", "stat", "z", "zscore", "auc"],
}


def _pick(cols: Dict[str, str], key: str) -> Optional[str]:
    for alias in _COLS[key]:
        if alias in cols:
            return cols[alias]
    return None


def read_marker_table(source, sep: Optional[str] = None) -> pd.DataFrame:
    """Read a marker table into columns ``cluster, gene, logfc, padj, score``.

    ``source`` is a path (CSV/TSV/XLSX) or a DataFrame. Seurat's ``FindAllMarkers``
    output (gene names in the row names or a ``gene`` column), scanpy's
    ``sc.get.rank_genes_groups_df(adata, None)`` and scanpy's wide export
    (``<group>_n``, ``<group>_l``, ``<group>_p`` columns) are recognized.
    """
    if isinstance(source, pd.DataFrame):
        raw = source.copy()
    else:
        path = str(source)
        if path.lower().endswith((".xlsx", ".xls")):
            raw = pd.read_excel(path)
        else:
            if sep is None:
                sep = "," if path.lower().endswith(".csv") else None
            raw = pd.read_csv(path, sep=sep, engine="python")
    cols = {str(c).strip().lower(): c for c in raw.columns}
    if _pick(cols, "cluster") is None and any(re.search(r"_(n|l|p|s)$", c) for c in cols):
        return _from_scanpy_wide(raw)

    cl, gene, lfc, padj, pval, score = (_pick(cols, k) for k in ["cluster", "gene", "logfc", "padj", "pval", "score"])
    if cl is None:
        raise ValueError("No cluster/group column found. Expected one of: " + ", ".join(_COLS["cluster"]))
    if gene is None:
        # Seurat: genes are row names, written as the first unnamed column.
        first = raw.columns[0]
        if str(first).lower().startswith("unnamed") or first == 0 or str(first).strip() == "":
            gene = first
        else:
            raise ValueError("No gene column found. Expected one of: " + ", ".join(_COLS["gene"]))
    out = pd.DataFrame({
        "cluster": raw[cl].astype(str).str.strip(),
        "gene": raw[gene].astype(str).map(lambda g: normalize_gene(re.sub(r"\.\d+$", "", g) if g.count(".") == 1 and g.rsplit(".", 1)[1].isdigit() else g)),
        "logfc": pd.to_numeric(raw[lfc], errors="coerce") if lfc else np.nan,
        "padj": pd.to_numeric(raw[padj], errors="coerce") if padj else (pd.to_numeric(raw[pval], errors="coerce") if pval else np.nan),
        "score": pd.to_numeric(raw[score], errors="coerce") if score else np.nan,
    })
    return out[out["gene"].str.len() > 0].reset_index(drop=True)


def _from_scanpy_wide(raw: pd.DataFrame) -> pd.DataFrame:
    """scanpy's ``rank_genes_groups`` wide export: ``<group>_n, _s, _l, _p``."""
    groups = sorted({re.sub(r"_(n|l|p|s|pa)$", "", str(c)) for c in raw.columns if re.search(r"_n$", str(c))})
    rows = []
    for g in groups:
        n = raw.get(f"{g}_n")
        if n is None:
            continue
        rows.append(pd.DataFrame({
            "cluster": g, "gene": n.astype(str).map(normalize_gene),
            "logfc": pd.to_numeric(raw.get(f"{g}_l", np.nan), errors="coerce"),
            "padj": pd.to_numeric(raw.get(f"{g}_pa", raw.get(f"{g}_p", np.nan)), errors="coerce"),
            "score": pd.to_numeric(raw.get(f"{g}_s", np.nan), errors="coerce"),
        }))
    return pd.concat(rows, ignore_index=True)


def markers_from_anndata(adata, key: str = "rank_genes_groups") -> pd.DataFrame:
    """Marker table from ``adata.uns[key]`` (after ``sc.tl.rank_genes_groups``)."""
    import scanpy as sc

    df = sc.get.rank_genes_groups_df(adata, None, key=key)
    return read_marker_table(df)


def markers_to_conditions(
    markers,
    padj_cutoff: Optional[float] = 0.05,
    min_logfc: Optional[float] = 0.25,
    top_n: Optional[int] = 200,
    direction: str = "up",
    clusters: Optional[Iterable[str]] = None,
    rank_by: str = "auto",
) -> Dict[str, Dict[str, Optional[float]]]:
    """``{cluster: {gene: log2FC}}`` from a marker table (path, DataFrame or AnnData).

    direction
        ``"up"`` keeps positive markers (the usual choice), ``"down"`` negative
        ones, ``"both"`` keeps both so up/down are tested separately.
    top_n
        keep the strongest ``top_n`` markers per cluster (ranked by ``score`` when
        present, else by |log2FC|); ``None`` keeps all that pass the filters.
    """
    if hasattr(markers, "uns"):
        df = markers_from_anndata(markers)
    elif isinstance(markers, pd.DataFrame) and {"cluster", "gene"} <= set(markers.columns) and "logfc" in markers.columns:
        df = markers.copy()
    else:
        df = read_marker_table(markers)
    if clusters is not None:
        wanted = {str(c) for c in clusters}
        df = df[df["cluster"].isin(wanted)]
        missing = wanted - set(df["cluster"])
        if missing:
            raise ValueError(f"Clusters not in the marker table: {sorted(missing)}. Available: {sorted(set(df['cluster']))}")
    if padj_cutoff is not None and df["padj"].notna().any():
        df = df[df["padj"].isna() | (df["padj"] <= padj_cutoff)]
    if df["logfc"].notna().any():
        if direction == "up":
            df = df[df["logfc"] > 0]
        elif direction == "down":
            df = df[df["logfc"] < 0]
        if min_logfc:
            df = df[df["logfc"].abs() >= min_logfc]
    if rank_by == "auto":
        rank_by = "score" if df["score"].notna().any() else "logfc"
    out: Dict[str, Dict[str, Optional[float]]] = {}
    for cl, grp in df.groupby("cluster", sort=False):
        grp = grp.assign(_r=grp[rank_by].abs().fillna(0)).sort_values("_r", ascending=False)
        if top_n:
            grp = grp.head(top_n)
        out[str(cl)] = {g: (None if pd.isna(v) else float(v)) for g, v in zip(grp["gene"], grp["logfc"]) if g}
    out = {k: v for k, v in out.items() if v}
    if not out:
        raise ValueError("No markers left after filtering. Relax padj_cutoff / min_logfc or check the columns.")
    return out
