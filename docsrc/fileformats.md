# File formats

## Input formats

**Gene list** (`snv run FILE`): text with one gene per line, or a CSV/TSV with the gene in the first column. A numeric column named `log2FoldChange`, `logFC`, `avg_log2FC`, `log2FC`, `stat`, `t`, `score` or `value` is used as the fold change (the first numeric column otherwise). Header optional.

**Ranked table for GSEA** (`snv run FILE --method gsea`): the same format, but with *all* measured genes and a statistic for each (DESeq2 `stat`, limma `t`, or log2FC). Higher = up-regulated.

**Enrichment table** (`snv run --enrichment FILE`): results from Enrichr (`Term`, `Adjusted P-value`, `Genes`), clusterProfiler (`Description`, `p.adjust`, `geneID`), g:Profiler (`term_name`, `p_value`, `intersections`), GSEA desktop (`NAME`, `NES`, `FDR q-val`), Metascape-like exports, or this tool's own `enrichment.csv`. Columns are detected by name. Give the matching gene sets (`--library` or `--gmt`) so pathway overlaps can be computed from the full gene sets.

**Pathway list** (`snv run --pathways FILE`): one pathway name per line. Matching ignores case and punctuation, so `G2M checkpoint` finds `HALLMARK_G2M_CHECKPOINT`.

**Several gene lists** (`snv compare NAME=FILE ...`): one file per condition, each in the gene-list format; or one wide table with a column per group (`snv compare --groups FILE`).

**Marker table** (`snv markers FILE`): Seurat `FindAllMarkers` output (`cluster`, `gene`, `avg_log2FC`, `p_val_adj`; gene names may be in the row names), scanpy `sc.get.rank_genes_groups_df(adata, None)` (`group`, `names`, `logfoldchanges`, `pvals_adj`), scanpy's wide export (`<group>_n`, `<group>_l`, `<group>_p`), or any table with cluster, gene, log fold-change and adjusted p-value columns.

**Gene sets**: Enrichr library names (several at once, comma-separated: `--library KEGG_2021_Human,Reactome_2022`), local `.gmt` files (`--gmt`, repeatable), or `demo`.

## Output files

| File | What it is |
|---|---|
| `pathway_network.pdf` / `.svg` | The figure, as vector graphics with editable text. Add `--formats png,pdf,svg,...` if you also want a PNG. |
| `pathway_network.html` | Interactive viewer: open in a browser |
| `pathway_map.html` | Interactive map (with `--map`): hover pathways and regions, lasso-select groups |
| `gene_network.*`, `gene_pathway_network.*` | Gene-level figures |
| `cluster_barplot.pdf` / `.svg` | Cluster themes ranked by significance (a companion panel) |
| `enrichment.csv` | Enrichment table: pathway, overlap, set size, fold enrichment, p-value, FDR, direction, genes |
| `pathway_network_clusters.csv` | Cluster number, automatic label, size, best FDR, member pathways |
| `pathway_network_nodes.csv`, `_edges.csv` | Everything behind the figure, for your own plots |
| `pathway_network.cyjs` + `_style.xml` | For fine-tuning in Cytoscape ([Fine-tuning in Cytoscape](guide/figures.md#fine-tuning-in-cytoscape)) |
| `pathway_network.graphml` | For igraph, Gephi and other network tools |
| `run_params.json` | Every setting of this run, the package version and the date |

Comparisons ([`snv compare`](commands/compare.md), [`snv markers`](commands/markers.md)) also write `enrichment_matrix.csv`, a pathway × condition table of signed −log10 FDR. Runs with `--map` also write `pathway_map.pdf`, `pathway_map.html`, `pathway_map_hits.csv` and `pathway_map_regions.csv`.

## Notes

* The bundled `demo` gene sets are shortened versions of real pathways, meant only for testing and demos. Use full libraries for real analyses.
* Figures use Arial at regular weight (Liberation Sans, with identical metrics, where Arial is not installed). PDF text is embedded as editable type (fonttype 42); SVG text stays text.
