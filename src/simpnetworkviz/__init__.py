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

Several conditions (contrasts, time points, single-cell clusters)::

    table = snv.compare_conditions({"KO": ko_genes, "WT": wt_genes}, gene_sets)
    pnet = snv.build_pathway_network(table, gene_sets)          # pie-chart nodes

Preranked GSEA::

    table = snv.gsea_prerank(ranked_scores, gene_sets)
"""
from .compare import compare_conditions, enrichment_matrix, merge_conditions, read_gene_groups
from .enrichment import bh_fdr, enrich, enrich_updown, read_enrichment, read_pathway_list, table_from_pathways
from .export import save_all, write_cyjs, write_graphml, write_html
from .genesets import (MSIGDB_COLLECTIONS, POPULAR_LIBRARIES, SPECIES, fetch_enrichr_library, fetch_msigdb, load_gene_sets,
                       normalize_species, read_gene_list, read_gmt, species_taxon, write_gmt)
from .gsea import gsea_prerank, read_ranking
from .network import (CLUSTER_METHODS, Network, annotate_cluster, build_gene_network, build_library_network,
                      build_pathway_network, clean_pathway_name, group_by_category, group_by_keywords, pathway_category)
from .pathmap import (build_category_maps, build_pathway_map, load_map, map_hits, plot_category_maps, plot_pathway_map,
                      region_summary, save_category_maps, save_map, split_by_category)
from .maphtml import write_map_html
from .pipeline import analyze_genes, run_pipeline
from .plotting import plot_cluster_barplot, plot_network, plot_pathway_gene_network
from .singlecell import markers_to_conditions, read_marker_table

__version__ = "1.0.0"

__all__ = [
    "CLUSTER_METHODS", "MSIGDB_COLLECTIONS", "Network", "POPULAR_LIBRARIES", "SPECIES", "analyze_genes", "annotate_cluster", "bh_fdr", "build_gene_network",
    "build_category_maps", "build_library_network", "build_pathway_map", "build_pathway_network", "clean_pathway_name", "compare_conditions", "enrich", "enrich_updown",
    "enrichment_matrix", "fetch_enrichr_library", "fetch_msigdb", "group_by_category", "group_by_keywords", "gsea_prerank", "load_gene_sets", "load_map", "map_hits",
    "markers_to_conditions", "merge_conditions", "normalize_species", "pathway_category", "plot_category_maps", "plot_cluster_barplot", "plot_network", "plot_pathway_map",
    "plot_pathway_gene_network", "read_enrichment", "read_gene_list", "read_gmt", "read_marker_table", "read_gene_groups", "read_pathway_list",
    "read_ranking", "region_summary", "run_pipeline", "save_all", "save_category_maps", "save_map", "species_taxon", "split_by_category", "table_from_pathways", "write_cyjs", "write_gmt", "write_graphml", "write_html", "write_map_html",
]
