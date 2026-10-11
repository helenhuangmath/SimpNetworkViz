"""Command line interface: ``simpnetworkviz`` (alias ``snv``)."""
from __future__ import annotations

import argparse
import os
import sys

from . import __version__
from .enrichment import read_enrichment, read_pathway_list
from .genesets import POPULAR_LIBRARIES, demo_path, load_gene_sets, read_gene_list
from .network import CLUSTER_METHODS, read_category_file, read_theme_file
from .pipeline import analyze_genes, run_pipeline

EPILOG = """examples:
  snv demo                                              # try it with bundled data
  snv run genes.txt --library MSigDB_Hallmark_2020      # genes -> enrichment -> networks
  snv run degs.csv --library KEGG_2021_Human,Reactome_2022 --background measured.txt
  snv run deseq2_all_genes.csv --method gsea --library MSigDB_Hallmark_2020   # preranked GSEA
  snv run --enrichment clusterprofiler_result.csv --gmt c2.cp.v2024.gmt
  snv run --pathways my_pathways.txt --library GO_Biological_Process_2023
  snv run genes.txt --library MSigDB_Hallmark_2020 --whole-library            # all 50 hallmarks, hits highlighted
  snv compare KO=ko_degs.csv WT=wt_degs.csv --library Reactome_2022           # pie-chart nodes
  snv markers seurat_markers.csv --library MSigDB_Hallmark_2020 --whole-library
  snv map --gmt m5.go.bp.v2025.1.Mm.symbols.gmt -o gobp_map.csv       # once per library
  snv run genes.txt --gmt m5.go.bp.v2025.1.Mm.symbols.gmt --map gobp_map.csv   # hits on the fixed map
  snv enrich genes.txt --gmt h.all.v2024.1.Hs.symbols.gmt -o enrichment.csv
  snv libraries
"""


def _add_network_args(p):
    p.add_argument("-o", "--outdir", default="snv_output", help="output folder (default: snv_output)")
    p.add_argument("--fdr", type=float, default=0.05, help="FDR cutoff for pathways (default 0.05)")
    p.add_argument("--max-pathways", type=int, default=100, help="max pathways in the network (default 100)")
    p.add_argument("--similarity", choices=["combined", "jaccard", "overlap"], default="combined",
                   help="gene-set similarity metric (default combined = EnrichmentMap default)")
    p.add_argument("--similarity-genes", choices=["full", "hits"], default="full",
                   help="compare full gene sets (EnrichmentMap default) or only the hit genes")
    p.add_argument("--cutoff", type=float, default=0.375, help="similarity cutoff for edges (default 0.375)")
    p.add_argument("--cluster", choices=list(CLUSTER_METHODS), default="mcl",
                   help="grouping: mcl (AutoAnnotate default), louvain, greedy, label_propagation, hierarchical, "
                        "components, category (by source/--categories), keyword (--themes), none")
    p.add_argument("--categories", help="TSV 'pathway<TAB>category' for --cluster category (default: name prefix)")
    p.add_argument("--themes", help="TSV 'theme<TAB>keyword1,keyword2' for --cluster keyword")
    p.add_argument("--resolution", type=float, default=1.0, help="louvain/greedy resolution (higher = more clusters)")
    p.add_argument("--inflation", type=float, default=2.0, help="MCL inflation (higher = more clusters)")
    p.add_argument("--linkage-threshold", type=float, default=None, help="hierarchical: max average distance, 1-similarity (default 1 - cutoff)")
    p.add_argument("--whole-library", action="store_true",
                   help="draw every pathway of the library and highlight the enriched ones")
    p.add_argument("--map", metavar="MAP.csv", default=None,
                   help="library map from 'snv map': also draw the enriched pathways on it (pathway_map.pdf)")
    p.add_argument("--color-by", default="auto",
                   help="auto | cluster | direction | significance | nes | condition (pathway network coloring)")
    p.add_argument("--gene-network", choices=["shared", "string", "cnet", "none"], default="shared",
                   help="gene network type: shared pathways (offline), STRING PPI (internet), cnet, or none")
    p.add_argument("--max-genes", type=int, default=150, help="max genes in gene network (default 150)")
    p.add_argument("--species", type=int, default=9606, help="NCBI taxon for STRING (9606 human, 10090 mouse)")
    p.add_argument("--formats", default="pdf,svg,html,cyjs,graphml,csv", help="comma separated output formats (add png for raster figures)")
    p.add_argument("--title", default=None, help="figure title")
    p.add_argument("--cytoscape-js", default=None, help="local cytoscape.min.js to embed (offline HTML)")
    p.add_argument("--seed", type=int, default=0, help="random seed for layouts")
    p.add_argument("--dpi", type=int, default=300, help="PNG resolution in dots per inch (default 300)")


def _add_enrich_args(p):
    p.add_argument("--library", default=None, help="Enrichr library name(s), comma separated, or 'demo'")
    p.add_argument("--gmt", action="append", help="GMT gene-set file (repeatable)")
    p.add_argument("--background", help="background/universe gene list (recommended: all measured genes)")
    p.add_argument("--method", choices=["ora", "gsea"], default="ora",
                   help="ora = hypergeometric test on a gene list; gsea = preranked GSEA on gene + statistic for all genes")
    p.add_argument("--no-split", action="store_true", help="do not split up/down genes by log2FC sign")
    p.add_argument("--min-size", type=int, default=None, help="min gene-set size (default 5 ora, 15 gsea)")
    p.add_argument("--max-size", type=int, default=None, help="max gene-set size (default 1000 ora, 500 gsea)")
    p.add_argument("--permutations", type=int, default=1000, help="gsea: gene permutations (default 1000)")


def _gene_sets_arg(args):
    if getattr(args, "gmt", None):
        return ",".join(args.gmt)
    return args.library


def _enrich_kwargs(args):
    kw = {}
    if args.min_size is not None:
        kw["min_size"] = args.min_size
    if args.max_size is not None:
        kw["max_size"] = args.max_size
    if args.method == "gsea":
        kw["n_perm"] = args.permutations
    return kw


def _pipeline_kwargs(args):
    cluster_kwargs = {"mcl_inflation": args.inflation, "louvain_resolution": args.resolution,
                      "hierarchical_threshold": args.linkage_threshold}
    if args.categories:
        cluster_kwargs["categories"] = read_category_file(args.categories)
    if args.themes:
        cluster_kwargs["themes"] = read_theme_file(args.themes)
    if args.cluster == "keyword" and not args.themes:
        sys.exit("error: --cluster keyword needs --themes FILE (theme<TAB>keyword1,keyword2 per line)")
    return dict(
        outdir=args.outdir, fdr_cutoff=args.fdr, max_pathways=args.max_pathways, similarity_cutoff=args.cutoff,
        similarity_metric=args.similarity, use_full_gene_sets=args.similarity_genes == "full", cluster_method=args.cluster,
        cluster_kwargs=cluster_kwargs, library_network=args.whole_library,
        gene_network_mode=None if args.gene_network == "none" else args.gene_network, color_by=args.color_by,
        formats=[f.strip() for f in args.formats.split(",") if f.strip()], title=args.title, js_path=args.cytoscape_js,
        seed=args.seed, dpi=args.dpi, gene_network_kwargs={"max_genes": args.max_genes, "string_species": args.species},
        pathway_map=args.map,
    )


def _require_gene_sets(args, needed=True):
    gene_sets = _gene_sets_arg(args)
    if gene_sets is None and needed:
        sys.exit("error: choose gene sets with --library NAME or --gmt FILE (or 'snv libraries' to list)")
    return gene_sets


def cmd_run(args):
    genes = read_gene_list(args.genes) if args.genes else None
    background = list(read_gene_list(args.background)) if args.background else None
    enrichment = read_enrichment(args.enrichment) if args.enrichment else None
    pathways = read_pathway_list(args.pathways) if args.pathways else None
    if not any([genes, enrichment is not None, pathways]):
        sys.exit("error: give a gene list, --enrichment TABLE or --pathways FILE (see snv run -h)")
    gene_sets = _require_gene_sets(args, needed=enrichment is None or args.whole_library)
    res = run_pipeline(
        genes=genes, enrichment=enrichment, pathways=pathways, gene_sets=gene_sets, background=background,
        method=args.method, split_updown=not args.no_split, enrich_kwargs=_enrich_kwargs(args), **_pipeline_kwargs(args),
    )
    _report(res, args.outdir)


def _parse_conditions(items):
    conds = {}
    for item in items:
        if "=" in item:
            name, path = item.split("=", 1)
        else:
            path = item
            name = os.path.splitext(os.path.basename(path))[0]
        if not os.path.exists(path):
            sys.exit(f"error: file not found: {path}")
        conds[name.strip()] = read_gene_list(path)
    return conds


def cmd_compare(args):
    if args.groups:
        from .compare import read_gene_groups

        conditions = read_gene_groups(args.groups)
        if args.conditions:
            conditions.update(_parse_conditions(args.conditions))
    else:
        conditions = _parse_conditions(args.conditions)
    if len(conditions) < 2:
        sys.exit("error: give at least two conditions as NAME=FILE, or --groups TABLE with one column per group")
    print(f"{len(conditions)} groups: " + ", ".join(f"{k} ({len(v)} genes)" for k, v in conditions.items()))
    background = list(read_gene_list(args.background)) if args.background else None
    res = run_pipeline(conditions=conditions, gene_sets=_require_gene_sets(args), background=background, method=args.method,
                       split_updown=not args.no_split, enrich_kwargs=_enrich_kwargs(args), **_pipeline_kwargs(args))
    _report(res, args.outdir)


def cmd_markers(args):
    from .singlecell import markers_to_conditions

    conditions = markers_to_conditions(
        args.markers, padj_cutoff=args.padj, min_logfc=args.min_logfc, top_n=args.top_genes, direction=args.direction,
        clusters=[c.strip() for c in args.clusters.split(",")] if args.clusters else None,
    )
    print(f"{len(conditions)} clusters: " + ", ".join(f"{k} ({len(v)} markers)" for k, v in conditions.items()))
    background = list(read_gene_list(args.background)) if args.background else None
    res = run_pipeline(conditions=conditions, gene_sets=_require_gene_sets(args), background=background, method="ora",
                       split_updown=args.direction == "both", enrich_kwargs=_enrich_kwargs(args), **_pipeline_kwargs(args))
    _report(res, args.outdir)


def _report(res, outdir):
    pnet = res["pathway_network"]
    ct = pnet.cluster_table()
    n_hl = sum(1 for _, d in pnet.graph.nodes(data=True) if d.get("highlight"))
    extra = f", {n_hl} enriched (highlighted)" if any("highlight" in d for _, d in pnet.graph.nodes(data=True)) else ""
    print(f"Pathway network: {pnet.graph.number_of_nodes()} pathways, {pnet.graph.number_of_edges()} edges, {len(ct)} clusters{extra}")
    if pnet.conditions:
        print(f"Conditions (pie slices): {', '.join(pnet.conditions)}")
    for _, r in ct.head(15).iterrows():
        print(f"  [{r['cluster']:>2}] {r['label']}  ({r['n_nodes']} pathways)")
    if len(ct) > 15:
        print(f"  ... {len(ct) - 15} more (see pathway_network_clusters.csv)")
    if "gene_network" in res:
        g = res["gene_network"].graph
        print(f"Gene network: {g.number_of_nodes()} genes, {g.number_of_edges()} edges")
    print(f"\nResults written to {os.path.abspath(outdir)}/")
    print("  open pathway_network.html in a browser for the interactive view;")
    print("  import *.cyjs into Cytoscape (File > Import > Network from File) with *_style.xml.")


def cmd_enrich(args):
    genes = read_gene_list(args.genes)
    background = list(read_gene_list(args.background)) if args.background else None
    table = analyze_genes(genes, _require_gene_sets(args), background=background, split_updown=not args.no_split,
                          method=args.method, **_enrich_kwargs(args))
    table.to_csv(args.output, index=False)
    sig = (table["fdr"] <= 0.05).sum()
    print(f"{len(table)} pathways tested with overlap, {sig} with FDR <= 0.05 -> {args.output}")


def cmd_demo(args):
    genes = read_gene_list(demo_path("demo_genes.tsv"))
    res = run_pipeline(genes=genes, gene_sets="demo", outdir=args.outdir, similarity_cutoff=0.3,
                       js_path=args.cytoscape_js, dpi=args.dpi)
    _report(res, args.outdir)
    if args.compare:
        up = {g: v for g, v in genes.items() if v is not None and v > 0}
        down = {g: v for g, v in genes.items() if v is not None and v < 0}
        out = os.path.join(args.outdir, "compare")
        res = run_pipeline(conditions={"Up genes": up, "Down genes": down, "All genes": genes}, gene_sets="demo",
                           outdir=out, similarity_cutoff=0.3, library_network=True,
                           js_path=args.cytoscape_js, dpi=args.dpi)
        _report(res, out)


def cmd_map(args):
    from .pathmap import build_pathway_map, plot_pathway_map, save_map

    gs = load_gene_sets(_require_gene_sets(args))
    print(f"{len(gs)} gene sets; building the map ({args.method}) ...")
    pmap = build_pathway_map(gs, method=args.method, min_size=args.min_size, max_size=args.max_size,
                             n_neighbors=args.neighbors, cluster_resolution=args.resolution, seed=args.seed)
    save_map(pmap, args.output)
    stem = os.path.splitext(args.output)[0]
    fig = plot_pathway_map(pmap, None, f"{stem}.pdf", dpi=args.dpi)
    import matplotlib.pyplot as plt

    plt.close(fig)
    regions = pmap[pmap["cluster"] > 0].groupby("cluster").agg(n=("pathway", "size"), label=("cluster_label", "first"))
    print(f"{len(pmap)} pathways placed with {pmap.attrs['method'].upper()}, {len(regions)} labelled regions -> {args.output}, {stem}.pdf")
    for c, r in regions.head(20).iterrows():
        print(f"  [{c:>2}] {r['label']}  ({r['n']} pathways)")
    if len(regions) > 20:
        print(f"  ... {len(regions) - 20} more")
    print("Use it with:  snv run genes.txt --library ... --map " + args.output)


def cmd_libraries(args):
    print("Popular Enrichr libraries (any Enrichr library name works with --library):\n")
    for k, v in POPULAR_LIBRARIES.items():
        print(f"  {k:<30} {v}")
    print("\nFull list: https://maayanlab.cloud/Enrichr/#libraries")
    print("MSigDB GMT files: https://www.gsea-msigdb.org/gsea/downloads.jsp (use --gmt FILE)")
    print("Offline test set: --library demo")


def cmd_gui(args):
    try:
        import streamlit  # noqa: F401
    except ImportError:
        sys.exit("The GUI needs streamlit:  pip install 'simpnetworkviz[gui]'")
    import subprocess

    app = os.path.join(os.path.dirname(__file__), "app.py")
    sys.exit(subprocess.call([sys.executable, "-m", "streamlit", "run", app]))


def build_parser():
    parser = argparse.ArgumentParser(
        prog="snv", description="Pathway enrichment + EnrichmentMap-style pathway and gene network plots.",
        epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", version=f"simpnetworkviz {__version__}")
    sub = parser.add_subparsers(dest="command")

    p = sub.add_parser("run", help="gene list / enrichment table / pathway list -> networks",
                       epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("genes", nargs="?", help="gene list (txt, one per line) or table with gene + log2FC/statistic columns")
    p.add_argument("--enrichment", help="existing enrichment table (this tool, Enrichr, clusterProfiler, g:Profiler, GSEA)")
    p.add_argument("--pathways", help="text file with pathway names (one per line) to draw directly")
    _add_enrich_args(p)
    _add_network_args(p)
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("compare", help="several conditions (NAME=FILE ...) -> one network with pie-chart nodes",
                       epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("conditions", nargs="*", metavar="NAME=FILE", help="gene lists / DE tables, one per condition")
    p.add_argument("--groups", metavar="TABLE", help="wide CSV/TSV/XLSX with one gene group per column (header = group name)")
    _add_enrich_args(p)
    _add_network_args(p)
    p.set_defaults(func=cmd_compare)

    p = sub.add_parser("markers", help="single-cell marker table (Seurat FindAllMarkers / scanpy) -> per-cluster networks",
                       epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("markers", help="marker table with cluster, gene, log2FC and adjusted p-value columns")
    p.add_argument("--padj", type=float, default=0.05, help="adjusted p-value cutoff for markers (default 0.05)")
    p.add_argument("--min-logfc", type=float, default=0.25, help="min |log2FC| for markers (default 0.25)")
    p.add_argument("--top-genes", type=int, default=200, help="max markers per cluster (default 200)")
    p.add_argument("--direction", choices=["up", "down", "both"], default="up", help="marker sign to keep (default up)")
    p.add_argument("--clusters", help="comma separated cluster names to include (default all)")
    _add_enrich_args(p)
    _add_network_args(p)
    p.set_defaults(func=cmd_markers)

    p = sub.add_parser("map", help="build a fixed 2-D map of a whole gene-set library (reuse it with --map)")
    p.add_argument("--library", default=None, help="Enrichr library name")
    p.add_argument("--gmt", action="append", help="GMT gene-set file (repeatable)")
    p.add_argument("-o", "--output", default="pathway_map.csv", help="map table (a PDF with the same stem is written too)")
    p.add_argument("--method", choices=["auto", "umap", "tsne"], default="auto", help="embedding (auto = UMAP if installed, else t-SNE)")
    p.add_argument("--min-size", type=int, default=10)
    p.add_argument("--max-size", type=int, default=500)
    p.add_argument("--neighbors", type=int, default=15, help="neighbours for the embedding and the regions (default 15)")
    p.add_argument("--resolution", type=float, default=1.0, help="region resolution (higher = more, smaller regions)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--dpi", type=int, default=300)
    p.set_defaults(func=cmd_map)

    p = sub.add_parser("enrich", help="only run enrichment (ORA or GSEA), no network")
    p.add_argument("genes")
    _add_enrich_args(p)
    p.add_argument("-o", "--output", default="enrichment.csv")
    p.set_defaults(func=cmd_enrich)

    p = sub.add_parser("demo", help="run the full workflow on bundled demo data (offline)")
    p.add_argument("-o", "--outdir", default="demo_output")
    p.add_argument("--compare", action="store_true", help="also run a three-condition comparison into demo_output/compare")
    p.add_argument("--cytoscape-js", default=None)
    p.add_argument("--dpi", type=int, default=300, help="PNG resolution in dots per inch (default 300)")
    p.set_defaults(func=cmd_demo)

    sub.add_parser("libraries", help="list popular gene-set libraries").set_defaults(func=cmd_libraries)
    sub.add_parser("gui", help="launch the point-and-click web app (needs streamlit)").set_defaults(func=cmd_gui)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 0
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
