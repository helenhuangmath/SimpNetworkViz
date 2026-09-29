"""SimpNetworkViz quick start (runs offline with the bundled demo data).

Swap "demo" for an Enrichr library name such as "MSigDB_Hallmark_2020" or
"GO_Biological_Process_2023", or for the path to an MSigDB .gmt file.
"""
import simpnetworkviz as snv

# 1. Genes (with optional log2 fold changes) and gene sets
genes = snv.read_gene_list(snv.genesets.demo_path("demo_genes.tsv"))   # {gene: log2FC}
gene_sets = snv.load_gene_sets("demo")

# 2. Pathway enrichment (up and down genes tested separately -> red/blue nodes)
table = snv.analyze_genes(genes, gene_sets)
table.to_csv("enrichment.csv", index=False)
print(table.head())

# 3. Pathway network: similarity edges, MCL clusters, automatic cluster labels
pnet = snv.build_pathway_network(table, gene_sets, fdr_cutoff=0.05, similarity_cutoff=0.3)
print(pnet.cluster_table()[["cluster", "label", "n_nodes"]])
pnet.plot("pathway_network.pdf", title="Pathway network")
pnet.to_html("pathway_network.html")        # interactive, open in a browser
pnet.to_cytoscape("pathway_network.cyjs")    # Cytoscape desktop (+ _style.xml)

# Rename clusters by hand if you like, then re-plot
pnet.group_labels[1] = "ECM remodelling"
pnet.plot("pathway_network_renamed.png", color_by="cluster")

# 4. Gene networks for all clusters (or pick some: clusters=[1, 2])
gnet = snv.build_gene_network(pnet, gene_sets, mode="shared", gene_values=genes)
gnet.plot("gene_network.pdf")
cnet = snv.build_gene_network(pnet, gene_sets, mode="cnet", gene_values=genes)
cnet.plot("gene_pathway_network.pdf")
# mode="string" pulls protein-protein interactions from STRING (needs internet)
