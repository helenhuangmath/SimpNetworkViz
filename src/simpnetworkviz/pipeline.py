"""One-call workflows used by the CLI and the GUI."""
from __future__ import annotations

import json
import os
from datetime import datetime
from typing import Dict, Iterable, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .compare import compare_conditions, enrichment_matrix, is_comparison
from .enrichment import enrich, enrich_updown, split_by_sign, table_from_pathways
from .genesets import GeneSets, load_gene_sets
from .network import Network, build_gene_network, build_library_network, build_pathway_network
from .plotting import plot_cluster_barplot


def analyze_genes(
    genes,
    gene_sets="demo",
    background: Optional[Iterable[str]] = None,
    split_updown: bool = True,
    method: str = "ora",
    **enrich_kwargs,
) -> pd.DataFrame:
    """Enrichment for a gene list or ``{gene: log2FC}`` dict.

    ``method="ora"``: hypergeometric test. With fold changes and
    ``split_updown=True``, up- and down-regulated genes are tested separately so
    the network can be colored red/blue.

    ``method="gsea"``: preranked GSEA on the ``{gene: score}`` values (all
    measured genes with their statistic, not only the significant ones).
    """
    gs = gene_sets if isinstance(gene_sets, dict) else load_gene_sets(gene_sets)
    if method == "gsea":
        from .gsea import gsea_prerank, has_scores

        if not has_scores(genes):
            raise ValueError("GSEA needs a ranked table: gene + statistic (log2FC, t, stat ...) for all measured genes.")
        return gsea_prerank({g: v for g, v in genes.items() if v is not None}, gs, **enrich_kwargs)
    if isinstance(genes, dict):
        parts = split_by_sign(genes) if split_updown else None
        if parts and parts[0] and parts[1]:
            return enrich_updown(parts[0], parts[1], gs, background=background, **enrich_kwargs)
        if parts and (parts[0] or parts[1]):
            enrich_kwargs.setdefault("direction", "up" if parts[0] else "down")
        genes = list(genes)
    return enrich(genes, gs, background=background, **enrich_kwargs)


def _merge_gene_values(conditions: Dict[str, object]) -> Optional[Dict[str, Optional[float]]]:
    """Union of the genes of all conditions; value = mean log2FC where available."""
    vals: Dict[str, list] = {}
    for genes in conditions.values():
        if isinstance(genes, pd.DataFrame):
            continue
        items = genes.items() if isinstance(genes, dict) else [(g, None) for g in genes]
        for g, v in items:
            vals.setdefault(g, [])
            if v is not None and v == v:
                vals[g].append(float(v))
    if not vals:
        return None
    return {g: (float(np.mean(v)) if v else None) for g, v in vals.items()}


def run_pipeline(
    genes=None,
    enrichment=None,
    pathways: Optional[Sequence[str]] = None,
    conditions: Optional[Dict[str, object]] = None,
    gene_sets="demo",
    outdir: str = "snv_output",
    background: Optional[Iterable[str]] = None,
    method: str = "ora",
    fdr_cutoff: float = 0.05,
    max_pathways: int = 100,
    similarity_cutoff: float = 0.375,
    similarity_metric: str = "combined",
    use_full_gene_sets: bool = True,
    cluster_method: str = "mcl",
    cluster_kwargs: Optional[dict] = None,
    library_network: bool = False,
    gene_network_mode: str = "shared",
    color_by: str = "auto",
    formats: Sequence[str] = ("pdf", "svg", "html", "cyjs", "graphml", "csv"),
    title: Optional[str] = None,
    js_path: Optional[str] = None,
    seed: int = 0,
    dpi: int = 300,
    split_updown: bool = True,
    gene_network_kwargs: Optional[dict] = None,
    enrich_kwargs: Optional[dict] = None,
    pathway_map=None,
) -> Dict[str, object]:
    """Gene list / enrichment table / pathway list / conditions -> networks + figures on disk.

    Exactly one of ``genes`` (list or ``{gene: log2FC}``), ``enrichment`` (table
    or path), ``pathways`` (names in ``gene_sets``) or ``conditions``
    (``{name: genes}``, see :func:`~simpnetworkviz.compare.compare_conditions`)
    is the starting point; ``genes`` may also be given with the others to color
    gene nodes.

    ``method``: ``"ora"`` or ``"gsea"`` (``genes``/``conditions`` must then carry
    scores for all measured genes). ``library_network=True`` draws every pathway
    of the library and highlights the enriched ones. ``cluster_kwargs`` passes
    ``categories``, ``themes``, ``louvain_resolution``, ``mcl_inflation`` or
    ``hierarchical_threshold`` to :func:`build_pathway_network`. ``pathway_map``
    (a map table or the CSV path from :func:`~simpnetworkviz.pathmap.build_pathway_map`)
    adds ``pathway_map.pdf`` with the enriched pathways highlighted on the fixed
    library map, plus ``pathway_map_hits.csv`` and ``pathway_map_regions.csv``.
    """
    os.makedirs(outdir, exist_ok=True)
    gs: Optional[GeneSets] = None
    if gene_sets is not None:
        gs = gene_sets if isinstance(gene_sets, dict) else load_gene_sets(gene_sets)
    gene_values = genes if isinstance(genes, dict) else ({g: None for g in genes} if genes else None)
    enrich_kwargs = dict(enrich_kwargs or {})

    if conditions:
        table = compare_conditions(conditions, gs, background=background, method=method, split_updown=split_updown, **enrich_kwargs)
        table.to_csv(os.path.join(outdir, "enrichment.csv"), index=False)
        if is_comparison(table):
            enrichment_matrix(table, fdr_cutoff=fdr_cutoff).to_csv(os.path.join(outdir, "enrichment_matrix.csv"))
        if gene_values is None:
            gene_values = _merge_gene_values(conditions)
    elif enrichment is not None:
        table = enrichment
    elif pathways is not None:
        if gs is None:
            raise ValueError("A pathway list needs gene sets (library name or GMT) to look up the pathway genes.")
        table = table_from_pathways(pathways, gs, query=list(gene_values) if gene_values else None)
        fdr_cutoff = None
    elif gene_values:
        table = analyze_genes(gene_values, gs, background=background, split_updown=split_updown, method=method, **enrich_kwargs)
        table.to_csv(os.path.join(outdir, "enrichment.csv"), index=False)
        if method == "gsea":
            # The gene network should show leading-edge genes, not every ranked gene.
            lead = {g for s in table.loc[table["fdr"] <= (fdr_cutoff or 1), "genes"] for g in str(s).split(";") if g}
            gene_values = {g: v for g, v in gene_values.items() if g in lead}
    else:
        raise ValueError("Provide genes, an enrichment table, a list of pathways or conditions.")

    net_kwargs = dict(similarity_cutoff=similarity_cutoff, similarity_metric=similarity_metric,
                      use_full_gene_sets=use_full_gene_sets, cluster_method=cluster_method, seed=seed, **(cluster_kwargs or {}))
    if library_network:
        if gs is None:
            raise ValueError("library_network=True needs gene sets (library name or GMT).")
        pnet = build_library_network(table, gs, fdr_cutoff=fdr_cutoff, max_pathways=max(max_pathways, 400), **net_kwargs)
    else:
        pnet = build_pathway_network(table, gs, fdr_cutoff=fdr_cutoff, max_pathways=max_pathways, **net_kwargs)
    files = {"pathway_network": pnet.save(outdir, "pathway_network", formats, title=title, color_by=color_by, js_path=js_path, dpi=dpi)}
    fig = plot_cluster_barplot(pnet, None, dpi=dpi)
    if fig is not None:
        for f in ("pdf", "svg", "png"):
            if f in formats:
                fig.savefig(os.path.join(outdir, f"cluster_barplot.{f}"), dpi=dpi, facecolor="white")
        plt.close(fig)

    result = {"enrichment": table, "pathway_network": pnet, "files": files}
    if pathway_map is not None:
        from .pathmap import load_map, map_hits, plot_pathway_map, region_summary

        pmap = pathway_map if isinstance(pathway_map, pd.DataFrame) else load_map(pathway_map)
        if not isinstance(table, pd.DataFrame):
            from .enrichment import read_enrichment

            table = read_enrichment(table)
        hits = map_hits(pmap, table, fdr_cutoff)
        hits.to_csv(os.path.join(outdir, "pathway_map_hits.csv"), index=False)
        region_summary(pmap, hits).to_csv(os.path.join(outdir, "pathway_map_regions.csv"), index=False)
        fig = plot_pathway_map(pmap, table, None, fdr_cutoff=fdr_cutoff, title=title, dpi=dpi)
        files["pathway_map"] = {}
        for f in ("pdf", "svg", "png"):
            if f in formats:
                p = os.path.join(outdir, f"pathway_map.{f}")
                fig.savefig(p, dpi=dpi, facecolor="white", bbox_inches="tight")
                files["pathway_map"][f] = p
        plt.close(fig)
        result["pathway_map_hits"] = hits
    if gene_network_mode:
        kwargs = dict(gene_network_kwargs or {})
        if library_network:
            # Only the highlighted (enriched) pathways contribute genes.
            kwargs.setdefault("pathways", [n for n, d in pnet.graph.nodes(data=True) if d.get("highlight")])
        gnet = build_gene_network(pnet, gs, mode=gene_network_mode, gene_values=gene_values, seed=seed, **kwargs)
        files["gene_network"] = gnet.save(outdir, "gene_network", formats, color_by="auto", js_path=js_path, dpi=dpi)
        result["gene_network"] = gnet
        if gene_network_mode != "cnet":
            kwargs.pop("max_genes", None)
            kwargs.pop("string_species", None)
            cnet = build_gene_network(pnet, gs, mode="cnet", gene_values=gene_values, seed=seed, **kwargs)
            files["gene_pathway_network"] = cnet.save(outdir, "gene_pathway_network", formats, color_by="auto", js_path=js_path, dpi=dpi)
            result["gene_pathway_network"] = cnet
    _write_params(outdir, locals())
    return result


def _write_params(outdir, params):
    """Record the settings of a run (``run_params.json``) for reproducibility."""
    from . import __version__

    skip = {"gs", "gene_values", "table", "pnet", "gnet", "cnet", "fig", "files", "result", "kwargs", "net_kwargs", "lead"}
    out = {"simpnetworkviz": __version__, "date": datetime.now().isoformat(timespec="seconds")}
    for k, v in params.items():
        if k in skip or k.startswith("_"):
            continue
        if isinstance(v, dict):
            v = {str(a): (f"<{len(b)} genes>" if isinstance(b, (dict, list, set)) else str(b)) for a, b in v.items()} if k == "conditions" else {str(a): str(b) for a, b in v.items()}
        elif isinstance(v, (list, tuple, set)):
            v = f"<{len(v)} items>" if k in ("genes", "pathways", "background") else list(v)
        elif isinstance(v, pd.DataFrame):
            v = f"<table with {len(v)} rows>"
        elif v is not None and not isinstance(v, (str, int, float, bool)):
            v = str(v)
        out[k] = v
    with open(os.path.join(outdir, "run_params.json"), "w") as fh:
        json.dump(out, fh, indent=2)
