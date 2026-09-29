"""SimpNetworkViz: pathway enrichment -> EnrichmentMap-style pathway & gene networks.

Quick start::

    import simpnetworkviz as snv

    gene_sets = snv.load_gene_sets("MSigDB_Hallmark_2020")   # or a .gmt path, or "demo"
    table = snv.enrich(my_genes, gene_sets)                    # ORA, BH-corrected
    pnet = snv.build_pathway_network(table, gene_sets)          # cluster + annotate
    pnet.plot("pathway_network.pdf")
    pnet.to_html("pathway_network.html")                        # interactive
    gnet = snv.build_gene_network(pnet, gene_sets, mode="shared")
    gnet.plot("gene_network.pdf")
"""
from .enrichment import bh_fdr, enrich, enrich_updown, read_enrichment, read_pathway_list, table_from_pathways
from .export import save_all, write_cyjs, write_graphml, write_html
from .genesets import POPULAR_LIBRARIES, fetch_enrichr_library, load_gene_sets, read_gene_list, read_gmt, write_gmt
from .network import Network, annotate_cluster, build_gene_network, build_pathway_network, clean_pathway_name
from .pipeline import analyze_genes, run_pipeline
from .plotting import plot_cluster_barplot, plot_network, plot_pathway_gene_network

__version__ = "0.1.0"

__all__ = [
    "Network", "POPULAR_LIBRARIES", "analyze_genes", "annotate_cluster", "bh_fdr", "build_gene_network",
    "build_pathway_network", "clean_pathway_name", "enrich", "enrich_updown", "fetch_enrichr_library",
    "load_gene_sets", "plot_cluster_barplot", "plot_network", "plot_pathway_gene_network", "read_enrichment", "read_gene_list", "read_gmt",
    "read_pathway_list", "run_pipeline", "save_all", "table_from_pathways", "write_cyjs", "write_gmt",
    "write_graphml", "write_html",
]
