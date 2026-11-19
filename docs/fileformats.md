# File formats

This page describes the files SimpNetworkViz reads and the tables it writes.

## Input files

### Gene list

Either a plain text file with one gene symbol per line:

```text
MKI67
TOP2A
CDK1
```

or a CSV/TSV table whose first column is the gene. If the table has a numeric column, it is read as the fold change:

```text
gene	log2FC
MCM2	1.34
MCM3	3.12
COL1A1	-2.40
```

* A column named `log2FC`, `logFC`, `log2FoldChange` or `avg_log2FC` is used as the fold change, so results from DESeq2, edgeR, limma and Seurat can be used after filtering to the significant genes. If there is no column with one of these names, the first numeric column is used.
* The first row is treated as a header when its first cell is a word such as `gene` or `symbol`.
* With fold changes, genes with positive and negative values are tested separately, and gene nodes are colored by fold change.
* Gene symbols are matched without regard to case.
* Lines starting with `#` are ignored.

:::{note}
Pass only the genes you consider differentially expressed. SimpNetworkViz does not filter a full results table by p-value; it treats every gene in the file as a hit. Give the full set of measured genes as the background instead.
:::

### Background gene list

The same format as a gene list; only the gene column is used. It should hold every gene that could have been a hit, such as all genes that passed the expression filter in an RNA-seq experiment.

### Gene sets (GMT)

A tab-separated file with one gene set per line: the name, a description, then the genes.

```text
HALLMARK_G2M_CHECKPOINT	description	AURKA	AURKB	BUB1	CCNA2	...
HALLMARK_E2F_TARGETS	description	AURKA	BRCA1	CDK1	MCM2	...
```

This is the format [MSigDB](https://www.gsea-msigdb.org/gsea/downloads.jsp) distributes. Choose the file with gene symbols, not Entrez IDs, unless your gene list uses Entrez IDs too.

### Enrichment table

A CSV, TSV or Excel table with one pathway per row. The columns are recognized by name:

| Source | Pathway column | Significance column | Genes column |
|---|---|---|---|
| SimpNetworkViz | `pathway` | `fdr` | `genes` |
| Enrichr | `Term` | `Adjusted P-value` | `Genes` |
| clusterProfiler | `Description` | `p.adjust` | `geneID` |
| g:Profiler | `term_name` | `p_value` | `intersections` |
| GSEA | `NAME` | `FDR q-val` | none |

* A `NES` column, as in GSEA results, is read as the normalized enrichment score and can be used to color the nodes.
* A simple two-column table of pathway and FDR also works.
* Pass the gene sets the table was computed from with `--gmt` or `--library`. They are required when the table has no genes column, and recommended otherwise.

### Pathway list

A text file with one pathway name per line:

```text
G2M checkpoint
E2F targets
Interferon gamma response
```

Names are matched against the gene sets ignoring case and punctuation, so `G2M checkpoint` matches `HALLMARK_G2M_CHECKPOINT`.

## Output tables

### enrichment.csv

Written by [`snv enrich`](commands/enrich.md), and by [`snv run`](commands/run.md) when it starts from a gene list.

| Column | Meaning |
|---|---|
| `pathway` | Gene-set name |
| `overlap` | Number of your genes in the gene set |
| `set_size` | Size of the gene set, within the background |
| `query_size` | Number of your genes, within the background |
| `background_size` | Size of the background |
| `fold_enrichment` | Observed overlap divided by the overlap expected by chance |
| `pvalue` | Hypergeometric p-value |
| `fdr` | Benjamini–Hochberg adjusted p-value |
| `direction` | `up` or `down`, when genes were split by fold change |
| `genes` | The overlapping genes, separated by `;` |

### PREFIX_nodes.csv

One row per node of a network.

| Column | Meaning |
|---|---|
| `id`, `name` | The node identifier: the pathway name or gene symbol |
| `label` | The cleaned-up name drawn in the figure |
| `type` | `pathway` or `gene` |
| `size_value` | The value that sets the node size |
| `cluster` | Cluster number; 0 means the node is in no cluster |
| `cluster_label` | The cluster's label |
| `x`, `y` | Layout coordinates |

Pathway nodes also have:

| Column | Meaning |
|---|---|
| `fdr`, `pvalue`, `nes`, `direction` | Carried over from the enrichment table |
| `set_size`, `overlap` | Gene-set size and number of your genes in it |
| `genes` | Your genes in the pathway, separated by `;` |

Gene nodes also have:

| Column | Meaning |
|---|---|
| `logfc` | The fold change from your gene list, if given |
| `n_pathways` | Number of drawn pathways the gene belongs to |
| `pathways` | Those pathways, separated by `; ` |

### PREFIX_edges.csv

| Column | Network | Meaning |
|---|---|---|
| `source`, `target` | all | The two nodes |
| `weight` | all | Edge weight: gene-set similarity in a pathway network, number of shared pathways in a `shared` gene network |
| `shared` | pathway | Number of genes the two pathways share |
| `shared_pathways` | gene (`shared`) | Number of pathways the two genes share |
| `similarity` | gene (`shared`) | Fraction of their pathways the two genes share (Jaccard) |

### pathway_network_clusters.csv

| Column | Meaning |
|---|---|
| `cluster` | Cluster number, ordered from the most significant cluster |
| `label` | The automatic label |
| `n_nodes` | Number of pathways in the cluster |
| `best_fdr` | The smallest FDR among its pathways |
| `members` | The member pathways, separated by `; ` |

## Network files

| File | Use |
|---|---|
| `PREFIX.html` | Interactive viewer for any browser |
| `PREFIX.cyjs` + `PREFIX_style.xml` | Cytoscape desktop; see [Fine-tune in Cytoscape](tutorial/figures.md#fine-tune-in-cytoscape) |
| `PREFIX.graphml` | Gephi, igraph, networkx and other network tools |
