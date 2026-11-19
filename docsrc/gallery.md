# Gallery of outputs

All figures come from `snv demo`, which runs offline on 88 genes with log2 fold changes and 32 shortened gene sets. The package writes vector PDF/SVG files; the previews below are PNG renderings of them. Click a figure for full size.

Pathway network, colored by direction (default) or by cluster:

| `--color-by direction` | `--color-by cluster` |
|---|---|
| [![](../docs/images/pathway_network.png)](../docs/images/pathway_network.png) | [![](../docs/images/pathway_network_by_cluster.png)](../docs/images/pathway_network_by_cluster.png) |

Gene network and gene–pathway ("cnet") network:

| Gene network | Gene–pathway network |
|---|---|
| [![](../docs/images/gene_network.png)](../docs/images/gene_network.png) | [![](../docs/images/gene_pathway_network.png)](../docs/images/gene_pathway_network.png) |

Pathway and gene network side by side (`plot_pathway_gene_network` in Python), sharing cluster colors and labels:

[![](../docs/images/pathway_gene_panel.png)](../docs/images/pathway_gene_panel.png)

Several conditions in one network ([Comparing several conditions](commands/compare.md)). Each pathway is a pie chart with one slice per condition; a white slice means not enriched there. Left: slice color = condition, drawn on the whole library so non-enriched pathways appear grey. Right: slices colored red/blue by direction.

| Condition colors, whole library | Direction colors |
|---|---|
| [![](../docs/images/compare_library.png)](../docs/images/compare_library.png) | [![](../docs/images/compare_direction.png)](../docs/images/compare_direction.png) |

GSEA ([Enrichment test: ORA or GSEA](commands/run.md#enrichment-test-ora-or-gsea)) colored by normalized enrichment score, and grouping by source database ([Grouping pathways](guide/grouping.md)):

| GSEA, colored by NES | Grouped by source |
|---|---|
| [![](../docs/images/gsea_network.png)](../docs/images/gsea_network.png) | [![](../docs/images/grouped_by_source.png)](../docs/images/grouped_by_source.png) |

Cluster bar plot, ranking the themes by the FDR of their best pathway:

[![](../docs/images/cluster_barplot.png)](../docs/images/cluster_barplot.png)

The pathway map ([The whole-library view and the pathway map](guide/library_map.md)): all 5,489 GO Biological Process terms embedded once from gene-set similarity, the enriched terms of a Cl13 vs Arm comparison (from `examples/Cl13_vs_Arm`) drawn as pies on top. The same 35 labelled regions appear in both; only the 2-D layout differs between t-SNE and UMAP.

| t-SNE (`--method tsne`) | UMAP (`--method umap`) |
|---|---|
| [![](../docs/images/pathway_map_tsne.png)](../docs/images/pathway_map_tsne.png) | [![](../docs/images/pathway_map_umap.png)](../docs/images/pathway_map_umap.png) |

The same two maps without hits (the output of `snv map`), and with the six Cl13 subgroups as pie slices:

| | t-SNE | UMAP |
|---|---|---|
| Library only | [![](../docs/images/pathway_map_tsne_library.png)](../docs/images/pathway_map_tsne_library.png) | [![](../docs/images/pathway_map_umap_library.png)](../docs/images/pathway_map_umap_library.png) |
| Six Cl13 subgroups | [![](../docs/images/pathway_map_tsne_within_Cl13.png)](../docs/images/pathway_map_tsne_within_Cl13.png) | [![](../docs/images/pathway_map_umap_within_Cl13.png)](../docs/images/pathway_map_umap_within_Cl13.png) |

The same test gene sets against the **MSigDB C7 ImmuneSigDB** collection (4,854 immunologic signatures, human v7.4 symbols; the maps are from `examples/msigdb_maps/c7.all`). ImmuneSigDB sets are up/down pairs of 200 genes from different GEO comparisons and rarely share genes, so the comparison network stays sparse (80 most significant sets, 19 edges) while the map shows where the hits fall among the T-cell, Treg, macrophage and dendritic-cell signature regions. Because T-cell gene lists hit almost half of this collection at FDR 0.05, the map overlays use FDR ≤ 1e-8 (300 sets).

| Cl13 vs Arm network, C7 | Cl13 vs Arm, t-SNE map, C7 | UMAP |
|---|---|---|
| [![](../docs/images/pathway_network_c7.png)](../docs/images/pathway_network_c7.png) | [![](../docs/images/pathway_map_tsne_c7.png)](../docs/images/pathway_map_tsne_c7.png) | [![](../docs/images/pathway_map_umap_c7.png)](../docs/images/pathway_map_umap_c7.png) |

| | t-SNE | UMAP |
|---|---|---|
| C7 library only | [![](../docs/images/pathway_map_tsne_c7_library.png)](../docs/images/pathway_map_tsne_c7_library.png) | [![](../docs/images/pathway_map_umap_c7_library.png)](../docs/images/pathway_map_umap_c7_library.png) |
| Six Cl13 subgroups on C7 | [![](../docs/images/pathway_map_tsne_c7_within_Cl13.png)](../docs/images/pathway_map_tsne_c7_within_Cl13.png) | [![](../docs/images/pathway_map_umap_c7_within_Cl13.png)](../docs/images/pathway_map_umap_c7_within_Cl13.png) |

**GO, KEGG, Reactome and WikiPathways** (MSigDB C5 and C2 collections). Per-database maps are in `examples/msigdb_maps/c5.all` and `c2.cp`; the overall map below embeds all six databases together (8,851 sets, 36 regions, `examples/msigdb_maps/overall`), points colored by database, with the Cl13 vs Arm hits on top (FDR ≤ 0.05):

| | t-SNE | UMAP |
|---|---|---|
| Overall library, by database | [![](../docs/images/pathway_map_tsne_overall_library.png)](../docs/images/pathway_map_tsne_overall_library.png) | [![](../docs/images/pathway_map_umap_overall_library.png)](../docs/images/pathway_map_umap_overall_library.png) |
| Cl13 vs Arm on the overall map | [![](../docs/images/pathway_map_tsne_overall.png)](../docs/images/pathway_map_tsne_overall.png) | [![](../docs/images/pathway_map_umap_overall.png)](../docs/images/pathway_map_umap_overall.png) |

**Every collection at once, human v7.4 collections.** `examples/Cl13_vs_Arm/run_msigdb_all.py --species human` reruns the same comparison against all sixteen human collection maps (summary table in [Step 6 of the Cl13 vs Arm tutorial](tutorial/cl13_vs_arm.md#step-6-place-every-analysis-on-the-pathway-map)). Each panel is one collection's map with the Cl13 vs Arm hits as pies. The mouse version, which is the right one for these mouse genes, follows below.

| t-SNE maps | UMAP maps |
|---|---|
| [![](../docs/images/msigdb_all_human_tsne.png)](../docs/images/msigdb_all_human_tsne.png) | [![](../docs/images/msigdb_all_human_umap.png)](../docs/images/msigdb_all_human_umap.png) |

**The same with the mouse collections**, which are the right ones for these mouse genes: `run_msigdb_all.py --species mouse` downloads MSigDB's native mouse collections (Hallmark, GO, Reactome, WikiPathways, BioCarta, cell types, TF targets; release 2025.1) through the package, uses their maps in `examples/msigdb_maps/mouse`, and adds the human C7 ImmuneSigDB map as the one cross-species exception. Mouse sets recover more genes and give more hits (GO BP 364 Cl13 and 22 Arm sets instead of 290 and 18; the full table is in the tutorial):

| t-SNE maps, mouse | UMAP maps, mouse |
|---|---|
| [![](../docs/images/msigdb_all_mouse_tsne.png)](../docs/images/msigdb_all_mouse_tsne.png) | [![](../docs/images/msigdb_all_mouse_umap.png)](../docs/images/msigdb_all_mouse_umap.png) |
