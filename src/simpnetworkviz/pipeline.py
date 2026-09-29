"""One-call workflows used by the CLI and the GUI."""
from __future__ import annotations

import os
from typing import Dict, Iterable, Optional, Sequence

import matplotlib.pyplot as plt
import pandas as pd

from .enrichment import enrich, enrich_updown, split_by_sign, table_from_pathways
from .genesets import GeneSets, load_gene_sets
from .network import Network, build_gene_network, build_pathway_network
from .plotting import plot_cluster_barplot


def analyze_genes(
    genes,
    gene_sets="demo",
    background: Optional[Iterable[str]] = None,
    split_updown: bool = True,
    **enrich_kwargs,
) -> pd.DataFrame:
    """Enrichment for a gene list or ``{gene: log2FC}`` dict.

    With fold changes and ``split_updown=True``, up- and down-regulated genes are
    tested separately so the network can be colored red/blue.
    """
    gs = gene_sets if isinstance(gene_sets, dict) else load_gene_sets(gene_sets)
    if isinstance(genes, dict):
        parts = split_by_sign(genes) if split_updown else None
        if parts and parts[0] and parts[1]:
            return enrich_updown(parts[0], parts[1], gs, background=background, **enrich_kwargs)
        genes = list(genes)
    return enrich(genes, gs, background=background, **enrich_kwargs)


def run_pipeline(
    genes=None,
    enrichment=None,
    pathways: Optional[Sequence[str]] = None,
    gene_sets="demo",
    outdir: str = "snv_output",
    background: Optional[Iterable[str]] = None,
    fdr_cutoff: float = 0.05,
    max_pathways: int = 100,
    similarity_cutoff: float = 0.375,
    similarity_metric: str = "combined",
    cluster_method: str = "mcl",
    gene_network_mode: str = "shared",
    color_by: str = "auto",
    formats: Sequence[str] = ("png", "pdf", "html", "cyjs", "graphml", "csv"),
    title: Optional[str] = None,
    js_path: Optional[str] = None,
    seed: int = 0,
    dpi: int = 300,
    split_updown: bool = True,
    gene_network_kwargs: Optional[dict] = None,
) -> Dict[str, object]:
    """Gene list / enrichment table / pathway list -> networks + figures on disk.

    Exactly one of ``genes`` (list or ``{gene: log2FC}``), ``enrichment`` (table
    or path) or ``pathways`` (names in ``gene_sets``) is the starting point;
    ``genes`` may also be given with the other two to color gene nodes.
    """
    os.makedirs(outdir, exist_ok=True)
    gs: Optional[GeneSets] = None
    if gene_sets is not None:
        gs = gene_sets if isinstance(gene_sets, dict) else load_gene_sets(gene_sets)
    gene_values = genes if isinstance(genes, dict) else ({g: None for g in genes} if genes else None)

    if enrichment is not None:
        table = enrichment
    elif pathways is not None:
        if gs is None:
            raise ValueError("A pathway list needs gene sets (library name or GMT) to look up the pathway genes.")
        table = table_from_pathways(pathways, gs, query=list(gene_values) if gene_values else None)
        fdr_cutoff = None
    elif gene_values:
        table = analyze_genes(gene_values, gs, background=background, split_updown=split_updown)
        table.to_csv(os.path.join(outdir, "enrichment.csv"), index=False)
    else:
        raise ValueError("Provide genes, an enrichment table or a list of pathways.")

    pnet = build_pathway_network(
        table, gs, fdr_cutoff=fdr_cutoff, max_pathways=max_pathways, similarity_cutoff=similarity_cutoff,
        similarity_metric=similarity_metric, cluster_method=cluster_method, seed=seed,
    )
    files = {"pathway_network": pnet.save(outdir, "pathway_network", formats, title=title, color_by=color_by, js_path=js_path, dpi=dpi)}
    fig = plot_cluster_barplot(pnet, os.path.join(outdir, "cluster_barplot.png") if "png" in formats else None, dpi=dpi)
    if fig is not None:
        if "pdf" in formats:
            fig.savefig(os.path.join(outdir, "cluster_barplot.pdf"))
        plt.close(fig)

    result = {"enrichment": table, "pathway_network": pnet, "files": files}
    if gene_network_mode:
        kwargs = dict(gene_network_kwargs or {})
        gnet = build_gene_network(pnet, gs, mode=gene_network_mode, gene_values=gene_values, seed=seed, **kwargs)
        files["gene_network"] = gnet.save(outdir, "gene_network", formats, color_by="auto", js_path=js_path, dpi=dpi)
        result["gene_network"] = gnet
        if gene_network_mode != "cnet":
            cnet = build_gene_network(pnet, gs, mode="cnet", gene_values=gene_values, seed=seed)
            files["gene_pathway_network"] = cnet.save(outdir, "gene_pathway_network", formats, color_by="auto", js_path=js_path, dpi=dpi)
            result["gene_pathway_network"] = cnet
    return result
