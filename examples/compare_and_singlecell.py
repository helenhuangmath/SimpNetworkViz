"""Compare conditions, run GSEA, and go from single-cell markers to networks.

Runs offline on the bundled demo data. Outputs go to ``example_output/``.
"""
import os

import numpy as np
import pandas as pd

import simpnetworkviz as snv
from simpnetworkviz.genesets import demo_path

out = "example_output"
os.makedirs(out, exist_ok=True)
gene_sets = snv.load_gene_sets("demo")               # use "MSigDB_Hallmark_2020" etc. for real data
genes = snv.read_gene_list(demo_path("demo_genes.tsv"))  # {gene: log2FC}

# 1. Several conditions in one network (pie-chart nodes) -------------------------
conditions = {
    "Up genes": {g: v for g, v in genes.items() if v > 0},
    "Down genes": {g: v for g, v in genes.items() if v < 0},
    "All genes": genes,
}
table = snv.compare_conditions(conditions, gene_sets)          # long table with a 'condition' column
snv.enrichment_matrix(table).to_csv(f"{out}/pathway_by_condition.csv")  # wide signed -log10 FDR matrix
pnet = snv.build_pathway_network(table, gene_sets, similarity_cutoff=0.3)
pnet.plot(f"{out}/compare_pathways.pdf", color_by="condition")  # slice color = condition
pnet.plot(f"{out}/compare_direction.pdf", color_by="direction")  # slice color = up/down in that condition
pnet.to_html(f"{out}/compare_pathways.html")

# 2. Whole library with the enriched pathways highlighted -------------------------
lib = snv.build_library_network(table, gene_sets, similarity_cutoff=0.3)
lib.plot(f"{out}/library_highlighted.pdf")

# 3. Preranked GSEA ----------------------------------------------------------------
# A real analysis passes ALL measured genes with their statistic (DESeq2 'stat', limma 't', ...).
rng = np.random.default_rng(0)
ranking = pd.Series(rng.normal(size=3000), index=[f"GENE{i}" for i in range(3000)])
ranking = pd.concat([ranking, pd.Series({g: v for g, v in genes.items()})])
gsea = snv.gsea_prerank(ranking, gene_sets, min_size=10)
gsea.to_csv(f"{out}/gsea.csv", index=False)
gnet = snv.build_pathway_network(gsea, gene_sets, fdr_cutoff=0.25, similarity_cutoff=0.3)
gnet.plot(f"{out}/gsea_network.pdf", color_by="nes")

# 4. Single-cell markers -> per-cluster enrichment -> one network ------------------
# Any Seurat FindAllMarkers / scanpy rank_genes_groups table works; here we fake one.
markers = pd.DataFrame({
    "gene": list(genes), "cluster": ["Cycling" if v > 0 else "Fibroblast-like" for v in genes.values()],
    "avg_log2FC": [abs(v) for v in genes.values()], "p_val_adj": 1e-6,
})
per_cluster = snv.markers_to_conditions(markers, top_n=100)      # {cluster: {gene: log2FC}}
sc_table = snv.compare_conditions(per_cluster, gene_sets)
sc_net = snv.build_library_network(sc_table, gene_sets, similarity_cutoff=0.3)
sc_net.plot(f"{out}/singlecell_library.pdf")

# 5. Other ways to group pathways ---------------------------------------------------
single = snv.analyze_genes(genes, gene_sets)
for method in ["louvain", "hierarchical", "category"]:
    net = snv.build_pathway_network(single, gene_sets, similarity_cutoff=0.3, cluster_method=method)
    net.plot(f"{out}/grouping_{method}.pdf")
themes = {"Proliferation": ["cycle", "mitotic", "E2F", "replication"], "Immune": ["interferon"],
          "Matrix": ["matrix", "collagen", "mesenchymal", "adhesion"]}
net = snv.build_pathway_network(single, gene_sets, similarity_cutoff=0.3, cluster_method="keyword", themes=themes)
net.plot(f"{out}/grouping_keyword.pdf")
print("done ->", out)
