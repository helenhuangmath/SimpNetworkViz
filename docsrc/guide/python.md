# Using the Python API

The command line is a thin layer over a Python API. Use it to rename clusters, choose subsets, or place the figures in your own layout.

```python
import simpnetworkviz as snv

genes = snv.read_gene_list("degs.csv")                  # {gene: log2FC} or plain list
gene_sets = snv.load_gene_sets("MSigDB_Hallmark_2020")  # Enrichr name, .gmt path, or "demo"
measured = list(snv.read_gene_list("measured.txt"))

table = snv.analyze_genes(genes, gene_sets, background=measured)   # enrichment table
pnet = snv.build_pathway_network(table, gene_sets)                 # similarity + MCL + labels
pnet.cluster_table()                                               # inspect the clusters
pnet.group_labels[3] = "Interferon signalling"                     # rename one
pnet.plot("pathway_network.pdf", color_by="direction")
pnet.to_html("pathway_network.html")

gnet = snv.build_gene_network(pnet, gene_sets, mode="shared", clusters=[1, 3], gene_values=genes)
gnet.plot("gene_network.pdf")
snv.plot_pathway_gene_network(pnet, gnet, "panel.pdf")           # a | b figure with both
```

Comparisons, GSEA and single cell:

```python
table = snv.compare_conditions({"KO": ko_genes, "WT": wt_genes}, gene_sets)   # pie-chart nodes
snv.enrichment_matrix(table)                                                   # pathway x condition
snv.build_pathway_network(table, gene_sets).plot("compare.pdf")
snv.build_library_network(table, gene_sets).plot("library.pdf")                # whole library, hits highlighted

gsea = snv.gsea_prerank(ranked_series, gene_sets)                              # NES, leading edge
snv.build_pathway_network(gsea, gene_sets, fdr_cutoff=0.25).plot("gsea.pdf", color_by="nes")

per_cluster = snv.markers_to_conditions("seurat_markers.csv")                  # or an AnnData object
table = snv.compare_conditions(per_cluster, gene_sets)

snv.build_pathway_network(table, gene_sets, cluster_method="keyword",
                          themes={"Proliferation": ["cycle", "E2F"], "Immune": ["interferon"]})
```

`run_pipeline(...)` does everything `snv run` does and returns the tables and networks. `examples/quickstart.py` and `examples/compare_and_singlecell.py` are complete, runnable walkthroughs.

Every function is listed in the [API reference](../api.md).
