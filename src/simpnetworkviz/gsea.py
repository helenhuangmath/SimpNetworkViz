"""Preranked gene-set enrichment analysis (GSEA).

A compact implementation of the weighted Kolmogorov-Smirnov statistic of
Subramanian et al. (2005) with a gene-permutation null, as in ``fgsea`` and
``gseapy.prerank``. Results come back in the package's standard enrichment
table so they feed straight into :func:`~simpnetworkviz.network.build_pathway_network`.
"""
from __future__ import annotations

from typing import Dict, Iterable, Optional, Union

import numpy as np
import pandas as pd

from .enrichment import STANDARD_COLUMNS, bh_fdr
from .genesets import GeneSets, normalize_gene

Ranking = Union[Dict[str, float], pd.Series, pd.DataFrame]


def _as_ranking(ranking: Ranking) -> pd.Series:
    """``{gene: score}``, Series or (gene, score) frame -> Series sorted high to low."""
    if isinstance(ranking, pd.DataFrame):
        from .genesets import gene_values_from_frame

        ranking = gene_values_from_frame(ranking)
    if isinstance(ranking, dict):
        ranking = pd.Series({normalize_gene(k): v for k, v in ranking.items()}, dtype=float)
    s = pd.Series(ranking, dtype=float)
    s.index = [normalize_gene(str(g)) for g in s.index]
    s = s[~s.index.duplicated()].dropna()
    if s.empty:
        raise ValueError("The ranking has no numeric scores. GSEA needs a gene table with a statistic (log2FC, t, stat ...).")
    # Stable sort so ties keep input order, as GSEA does.
    return s.iloc[np.argsort(-s.values, kind="stable")]


def _es_from_hits(idx: np.ndarray, w: np.ndarray, N: int) -> np.ndarray:
    """Enrichment score for hit positions ``idx`` (sorted, shape ``(..., K)``).

    The running sum only changes at hits, so it is enough to look at the
    deviation just before and just after every hit.
    """
    K = idx.shape[-1]
    j = np.arange(K)
    phit = np.cumsum(w, axis=-1) / np.maximum(w.sum(axis=-1, keepdims=True), 1e-300)
    pmiss = (idx - j) / max(N - K, 1)
    after = phit - pmiss
    before = np.concatenate([np.zeros(idx.shape[:-1] + (1,)), phit[..., :-1]], axis=-1) - pmiss
    dev = np.concatenate([after, before], axis=-1)
    pick = np.argmax(np.abs(dev), axis=-1)
    return np.take_along_axis(dev, pick[..., None], axis=-1)[..., 0]


def gsea_prerank(
    ranking: Ranking,
    gene_sets: GeneSets,
    min_size: int = 15,
    max_size: int = 500,
    n_perm: int = 1000,
    weight: float = 1.0,
    seed: int = 0,
    fdr_method: str = "bh",
) -> pd.DataFrame:
    """Preranked GSEA.

    Parameters
    ----------
    ranking : ``{gene: score}``, a Series indexed by gene, or a table whose
        first column is the gene and which has a statistic column (``stat``,
        ``t``, ``log2FoldChange`` ...). Higher scores = up-regulated.
    gene_sets : ``{pathway: [genes]}``
    min_size, max_size : size limits after restricting a set to ranked genes
    n_perm : gene permutations for the null distribution (per set size)
    weight : exponent on the scores (1 = classic weighted GSEA, 0 = plain KS)
    fdr_method : ``"bh"`` (Benjamini-Hochberg over the empirical p-values, as
        ``fgsea`` reports) or ``"gsea"`` (the original NES-based FDR q-value)

    Returns the standard enrichment table. ``genes`` holds the leading-edge
    genes, ``overlap`` their number, ``direction`` the sign of the NES.
    """
    s = _as_ranking(ranking)
    genes = np.array(s.index)
    scores = np.abs(s.values.astype(float)) ** weight
    N = len(genes)
    pos = {g: i for i, g in enumerate(genes)}
    rng = np.random.default_rng(seed)

    sets = {}
    for name, members in gene_sets.items():
        idx = np.array(sorted({pos[g] for g in members if g in pos}), dtype=int)
        if min_size <= len(idx) <= max_size:
            sets[name] = idx
    if not sets:
        raise ValueError("No gene set overlaps the ranked genes within the size limits. Check gene symbols / --min-size.")

    # Null ES depend only on set size under gene permutation: share them per size.
    null: Dict[int, np.ndarray] = {}
    for K in sorted({len(v) for v in sets.values()}):
        perm = np.sort(np.stack([rng.choice(N, K, replace=False) for _ in range(n_perm)]), axis=1)
        null[K] = _es_from_hits(perm, scores[perm], N)

    rows = []
    for name, idx in sets.items():
        K = len(idx)
        es = float(_es_from_hits(idx, scores[idx], N))
        nl = null[K]
        same = nl[nl >= 0] if es >= 0 else nl[nl < 0]
        denom = float(np.abs(same).mean()) if len(same) else 1.0
        nes = es / denom if denom > 0 else 0.0
        extreme = int((same >= es).sum()) if es >= 0 else int((same <= es).sum())
        p = (extreme + 1) / (len(same) + 1)
        # Leading edge: hits up to the peak (positive ES) or from the peak on (negative ES).
        phit = np.cumsum(scores[idx]) / max(scores[idx].sum(), 1e-300)
        pmiss = (idx - np.arange(K)) / max(N - K, 1)
        dev_after, dev_before = phit - pmiss, np.concatenate([[0.0], phit[:-1]]) - pmiss
        if es >= 0:
            peak = int(np.argmax(dev_after))
            lead = genes[idx[: peak + 1]]
        else:
            peak = int(np.argmin(dev_before))
            lead = genes[idx[peak:]]
        rows.append({
            "pathway": name, "overlap": len(lead), "set_size": K, "query_size": N, "background_size": N,
            "fold_enrichment": np.nan, "pvalue": p, "direction": "up" if nes > 0 else "down",
            "genes": ";".join(lead), "es": es, "nes": nes, "nes_null": (nl, denom),
        })
    df = pd.DataFrame(rows)
    if fdr_method == "gsea":
        df["fdr"] = _gsea_fdr(df)
    else:
        df["fdr"] = bh_fdr(df["pvalue"].values)
    df = df.drop(columns="nes_null")
    cols = STANDARD_COLUMNS + ["es", "nes"]
    return df[cols].sort_values(["fdr", "pvalue"]).reset_index(drop=True)


def _gsea_fdr(df: pd.DataFrame) -> np.ndarray:
    """Original GSEA FDR q-value: compare observed NES to the pooled null NES."""
    null_nes = np.concatenate([nl / np.where(nl >= 0, d, np.abs(nl[nl < 0]).mean() if (nl < 0).any() else 1.0)
                               for nl, d in df["nes_null"]]) if len(df) else np.array([])
    obs = df["nes"].values
    q = np.ones(len(df))
    for i, n in enumerate(obs):
        if n >= 0:
            frac_null = (null_nes[null_nes >= 0] >= n).mean() if (null_nes >= 0).any() else 1.0
            frac_obs = (obs[obs >= 0] >= n).mean()
        else:
            frac_null = (null_nes[null_nes < 0] <= n).mean() if (null_nes < 0).any() else 1.0
            frac_obs = (obs[obs < 0] <= n).mean()
        q[i] = min(1.0, frac_null / max(frac_obs, 1e-12))
    return q


def read_ranking(path) -> pd.Series:
    """Read a ranked gene table (gene + statistic) into a sorted Series."""
    from .genesets import read_gene_list

    values = read_gene_list(path)
    if not any(v is not None for v in values.values()):
        raise ValueError(f"{path}: no numeric column found. GSEA needs a statistic such as log2FC, t or stat.")
    return _as_ranking({g: v for g, v in values.items() if v is not None})


def has_scores(genes: Optional[Iterable]) -> bool:
    """True when a ``{gene: value}`` dict carries numeric values."""
    return isinstance(genes, dict) and any(v is not None and v == v for v in genes.values())
