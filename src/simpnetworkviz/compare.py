"""Enrichment of several conditions (contrasts, time points, cell clusters) side by side.

The result is one long table with a ``condition`` column. Passing it to
:func:`~simpnetworkviz.network.build_pathway_network` gives a network in which
every pathway node is a pie chart with one slice per condition, as in
EnrichmentMap with several data sets.
"""
from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Union

import numpy as np
import pandas as pd

from .enrichment import enrich, enrich_updown, read_enrichment, split_by_sign
from .genesets import GeneSets, load_gene_sets, normalize_gene

Condition = Union[Iterable[str], Dict[str, Optional[float]], pd.DataFrame, str]


def compare_conditions(
    conditions: Dict[str, Condition],
    gene_sets="demo",
    background: Optional[Iterable[str]] = None,
    method: str = "ora",
    split_updown: bool = True,
    **kwargs,
) -> pd.DataFrame:
    """Run enrichment for every condition and stack the results.

    Parameters
    ----------
    conditions : ``{name: genes}`` where genes is a list, a ``{gene: log2FC}``
        dict, a ranked Series/dict (``method="gsea"``), an enrichment table
        (DataFrame or path, used as is) or a path to a gene list.
    method : ``"ora"`` (hypergeometric test) or ``"gsea"`` (preranked GSEA;
        every condition then needs scores)
    kwargs : passed to :func:`enrich` or :func:`gsea_prerank`

    Returns the standard enrichment table plus a ``condition`` column.
    """
    gs = gene_sets if isinstance(gene_sets, dict) else load_gene_sets(gene_sets)
    parts = []
    for name, genes in conditions.items():
        if isinstance(genes, str):
            from .genesets import read_gene_list

            genes = read_gene_list(genes)
        if isinstance(genes, pd.DataFrame):
            table = genes if {"pathway", "fdr", "genes"} <= set(genes.columns) else read_enrichment(genes)
        elif method == "gsea":
            from .gsea import gsea_prerank

            table = gsea_prerank(genes, gs, **kwargs)
        else:
            table = _ora(genes, gs, background, split_updown, **kwargs)
        table = table.copy()
        table["condition"] = str(name)
        parts.append(table)
    if not parts:
        raise ValueError("No conditions given.")
    names = [str(n) for n in conditions]
    out = pd.concat([p for p in parts if len(p)] or parts[:1], ignore_index=True, sort=False)
    out.attrs["conditions"] = names  # kept even for conditions with no tested pathway
    return out


def _ora(genes, gs, background, split_updown, **kwargs):
    if isinstance(genes, pd.Series):
        genes = genes.to_dict()
    if isinstance(genes, dict):
        parts = split_by_sign(genes) if split_updown else None
        if parts and parts[0] and parts[1]:
            return enrich_updown(parts[0], parts[1], gs, background=background, **kwargs)
        if parts and (parts[0] or parts[1]):
            # All genes have the same sign (e.g. positive cluster markers): keep the direction.
            kwargs.setdefault("direction", "up" if parts[0] else "down")
        genes = list(genes)
    return enrich(genes, gs, background=background, **kwargs)


def read_gene_groups(path, sep: Optional[str] = None) -> Dict[str, Dict[str, Optional[float]]]:
    """Read a wide table with one gene group per column (header = group name).

    Columns may have different lengths (empty cells are ignored). Returns
    ``{group: {gene: None}}`` ready for :func:`compare_conditions`.
    """
    p = str(path)
    if p.lower().endswith((".xlsx", ".xls")):
        df = pd.read_excel(p, dtype=str)
    else:
        df = pd.read_csv(p, sep=sep or ("," if p.lower().endswith(".csv") else None), engine="python", dtype=str)
    groups = {}
    for col in df.columns:
        genes = [normalize_gene(g) for g in df[col].dropna().astype(str) if g.strip()]
        if genes:
            groups[str(col).strip()] = dict.fromkeys(genes)
    if not groups:
        raise ValueError(f"{path}: no genes found. Expected one column per group with the group name as header.")
    return groups


def condition_names(table: pd.DataFrame) -> List[str]:
    """Condition names in a comparison table (``table.attrs`` first, else order of appearance)."""
    if table.attrs.get("conditions"):
        return list(table.attrs["conditions"])
    if "condition" not in table.columns:
        return []
    return [c for c in pd.unique(table["condition"].astype(str)) if c and c != "nan"]


def is_comparison(table: pd.DataFrame) -> bool:
    return len(condition_names(table)) > 1


def enrichment_matrix(table: pd.DataFrame, value: str = "fdr", fdr_cutoff: Optional[float] = 0.05,
                      transform: Optional[str] = "neglog10") -> pd.DataFrame:
    """Wide pathway x condition matrix (for heatmaps or your own plots).

    ``value`` is ``"fdr"``, ``"nes"``, ``"fold_enrichment"`` or any column. With
    ``transform="neglog10"`` FDRs become -log10 values, signed by direction when
    a direction is known. Pathways not significant in any condition are dropped
    when ``fdr_cutoff`` is set.
    """
    df = table.copy()
    conds = condition_names(df)
    if fdr_cutoff is not None:
        keep = df.loc[df["fdr"] <= fdr_cutoff, "pathway"].unique()
        df = df[df["pathway"].isin(keep)]
    if value == "fdr" and transform == "neglog10":
        v = -np.log10(df["fdr"].clip(lower=1e-300))
        sign = df["direction"].map({"up": 1, "down": -1}).fillna(1) if "direction" in df else 1
        df["_v"] = v * sign + 0.0  # avoid -0.0
    else:
        df["_v"] = pd.to_numeric(df[value], errors="coerce")
    # A pathway tested both up and down in the same condition keeps the better entry.
    df = df.sort_values("fdr").drop_duplicates(["pathway", "condition"])
    mat = df.pivot(index="pathway", columns="condition", values="_v").reindex(columns=conds)
    order = mat.abs().max(axis=1).sort_values(ascending=False).index
    return mat.loc[order]


def merge_conditions(table: pd.DataFrame, fdr_cutoff: Optional[float] = 0.05,
                     conditions: Optional[List[str]] = None) -> pd.DataFrame:
    """Collapse a comparison table to one row per pathway.

    Adds ``fdr__<condition>``, ``direction__<condition>``, ``nes__<condition>``
    and ``genes__<condition>`` columns, ``conditions`` (names in which the
    pathway passes ``fdr_cutoff``) and ``n_conditions``. ``fdr`` is the best FDR,
    ``genes`` the union of hit genes and ``direction`` ``up``/``down`` when all
    conditions agree, else ``mixed``.
    """
    conds = list(conditions) if conditions else condition_names(table)
    if not conds:
        return table
    df = table.copy()
    df["condition"] = df["condition"].astype(str)
    df = df.sort_values("fdr", na_position="last").drop_duplicates(["pathway", "condition"])
    rows = []
    for pw, grp in df.groupby("pathway", sort=False):
        by = grp.set_index("condition")
        row = {"pathway": pw}
        sig, dirs, genes = [], set(), []
        for c in conds:
            if c in by.index:
                r = by.loc[c]
                row[f"fdr__{c}"] = r["fdr"]
                row[f"direction__{c}"] = r.get("direction", "")
                row[f"nes__{c}"] = r.get("nes", np.nan)
                row[f"genes__{c}"] = r.get("genes", "")
                if fdr_cutoff is None or (pd.notna(r["fdr"]) and r["fdr"] <= fdr_cutoff):
                    sig.append(c)
                    if r.get("direction") in ("up", "down"):
                        dirs.add(r["direction"])
                for g in str(r.get("genes", "") or "").split(";"):
                    if g and g not in genes:
                        genes.append(normalize_gene(g))
            else:
                row[f"fdr__{c}"] = np.nan
                row[f"direction__{c}"] = ""
                row[f"nes__{c}"] = np.nan
                row[f"genes__{c}"] = ""
        best = grp.iloc[0]
        for col in grp.columns:
            if col not in row and col != "condition":
                row[col] = best[col]
        row["fdr"] = grp["fdr"].min()
        row["genes"] = ";".join(genes)
        row["overlap"] = len(genes)
        row["direction"] = dirs.pop() if len(dirs) == 1 else ("mixed" if dirs else "")
        row["conditions"] = ";".join(sig)
        row["n_conditions"] = len(sig)
        rows.append(row)
    out = pd.DataFrame(rows)
    out.attrs["conditions"] = conds
    return out.sort_values(["fdr"], na_position="last").reset_index(drop=True)
