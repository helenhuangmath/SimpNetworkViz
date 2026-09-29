"""Command line interface: ``simpnetworkviz`` (alias ``snv``)."""
from __future__ import annotations

import argparse
import os
import sys

from . import __version__
from .enrichment import read_enrichment, read_pathway_list
from .genesets import POPULAR_LIBRARIES, demo_path, read_gene_list
from .pipeline import analyze_genes, run_pipeline

EPILOG = """examples:
  snv demo                                              # try it with bundled data
  snv run genes.txt --library MSigDB_Hallmark_2020      # genes -> enrichment -> networks
  snv run degs.csv --library KEGG_2021_Human,Reactome_2022 --background measured.txt
  snv run --enrichment clusterprofiler_result.csv --gmt c2.cp.v2024.gmt
  snv run --pathways my_pathways.txt --library GO_Biological_Process_2023
  snv enrich genes.txt --gmt h.all.v2024.1.Hs.symbols.gmt -o enrichment.csv
  snv libraries
"""


def _add_network_args(p):
    p.add_argument("-o", "--outdir", default="snv_output", help="output folder (default: snv_output)")
    p.add_argument("--fdr", type=float, default=0.05, help="FDR cutoff for pathways (default 0.05)")
    p.add_argument("--max-pathways", type=int, default=100, help="max pathways in the network (default 100)")
    p.add_argument("--similarity", choices=["combined", "jaccard", "overlap"], default="combined",
                   help="gene-set similarity metric (default combined = EnrichmentMap default)")
    p.add_argument("--cutoff", type=float, default=0.375, help="similarity cutoff for edges (default 0.375)")
    p.add_argument("--cluster", choices=["mcl", "louvain", "components", "none"], default="mcl",
                   help="clustering algorithm (default mcl, as in AutoAnnotate)")
    p.add_argument("--color-by", default="auto",
                   help="auto | cluster | direction | significance | nes (pathway network coloring)")
    p.add_argument("--gene-network", choices=["shared", "string", "cnet", "none"], default="shared",
                   help="gene network type: shared pathways (offline), STRING PPI (internet), cnet, or none")
    p.add_argument("--max-genes", type=int, default=150, help="max genes in gene network (default 150)")
    p.add_argument("--species", type=int, default=9606, help="NCBI taxon for STRING (9606 human, 10090 mouse)")
    p.add_argument("--formats", default="png,pdf,svg,html,cyjs,graphml,csv", help="comma separated output formats")
    p.add_argument("--title", default=None, help="figure title")
    p.add_argument("--cytoscape-js", default=None, help="local cytoscape.min.js to embed (offline HTML)")
    p.add_argument("--seed", type=int, default=0, help="random seed for layouts")
    p.add_argument("--dpi", type=int, default=300, help="PNG resolution in dots per inch (default 300)")


def _gene_sets_arg(args):
    if getattr(args, "gmt", None):
        return ",".join(args.gmt)
    return args.library


def cmd_run(args):
    genes = read_gene_list(args.genes) if args.genes else None
    background = list(read_gene_list(args.background)) if args.background else None
    enrichment = read_enrichment(args.enrichment) if args.enrichment else None
    pathways = read_pathway_list(args.pathways) if args.pathways else None
    if not any([genes, enrichment is not None, pathways]):
        sys.exit("error: give a gene list, --enrichment TABLE or --pathways FILE (see snv run -h)")
    gene_sets = _gene_sets_arg(args)
    if gene_sets is None and enrichment is None:
        sys.exit("error: choose gene sets with --library NAME or --gmt FILE (or 'snv libraries' to list)")
    res = run_pipeline(
        genes=genes, enrichment=enrichment, pathways=pathways, gene_sets=gene_sets, outdir=args.outdir,
        background=background, fdr_cutoff=args.fdr, max_pathways=args.max_pathways,
        similarity_cutoff=args.cutoff, similarity_metric=args.similarity, cluster_method=args.cluster,
        gene_network_mode=None if args.gene_network == "none" else args.gene_network, color_by=args.color_by,
        formats=[f.strip() for f in args.formats.split(",") if f.strip()], title=args.title,
        js_path=args.cytoscape_js, seed=args.seed, dpi=args.dpi, split_updown=not args.no_split,
        gene_network_kwargs={"max_genes": args.max_genes, "string_species": args.species},
    )
    _report(res, args.outdir)


def _report(res, outdir):
    pnet = res["pathway_network"]
    ct = pnet.cluster_table()
    print(f"Pathway network: {pnet.graph.number_of_nodes()} pathways, {pnet.graph.number_of_edges()} edges, {len(ct)} clusters")
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
    table = analyze_genes(genes, _gene_sets_arg(args), background=background, split_updown=not args.no_split,
                          min_size=args.min_size, max_size=args.max_size)
    table.to_csv(args.output, index=False)
    sig = (table["fdr"] <= 0.05).sum()
    print(f"{len(table)} pathways tested with overlap, {sig} with FDR <= 0.05 -> {args.output}")


def cmd_demo(args):
    genes = read_gene_list(demo_path("demo_genes.tsv"))
    res = run_pipeline(genes=genes, gene_sets="demo", outdir=args.outdir, similarity_cutoff=0.3,
                       title="Demo pathway network", js_path=args.cytoscape_js, dpi=args.dpi)
    _report(res, args.outdir)


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
    p.add_argument("genes", nargs="?", help="gene list (txt, one per line) or table with gene + log2FC columns")
    p.add_argument("--enrichment", help="existing enrichment table (this tool, Enrichr, clusterProfiler, g:Profiler, GSEA)")
    p.add_argument("--pathways", help="text file with pathway names (one per line) to draw directly")
    p.add_argument("--library", default=None, help="Enrichr library name(s), comma separated, or 'demo'")
    p.add_argument("--gmt", action="append", help="GMT gene-set file (repeatable)")
    p.add_argument("--background", help="background/universe gene list (recommended: all measured genes)")
    p.add_argument("--no-split", action="store_true", help="do not split up/down genes by log2FC sign")
    _add_network_args(p)
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("enrich", help="only run enrichment (over-representation analysis)")
    p.add_argument("genes")
    p.add_argument("--library", default=None)
    p.add_argument("--gmt", action="append")
    p.add_argument("--background")
    p.add_argument("--min-size", type=int, default=5)
    p.add_argument("--max-size", type=int, default=1000)
    p.add_argument("--no-split", action="store_true")
    p.add_argument("-o", "--output", default="enrichment.csv")
    p.set_defaults(func=cmd_enrich)

    p = sub.add_parser("demo", help="run the full workflow on bundled demo data (offline)")
    p.add_argument("-o", "--outdir", default="demo_output")
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
    if getattr(args, "library", None) is None and not getattr(args, "gmt", None) and args.command == "enrich":
        sys.exit("error: choose gene sets with --library NAME or --gmt FILE")
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
