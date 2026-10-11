# SimpNetworkViz: pathway enrichment and EnrichmentMap-style network figures from a gene list

**Author:** Hua Huang
**Package version:** 0.2.0
**Last updated:** 2026-10 

## Abstract

Pathway enrichment analysis of transcriptomic data typically yields tens to hundreds of overlapping gene sets, and the enrichment map, a network in which gene sets sharing genes are connected and related sets are grouped into annotated clusters, has become the standard way to summarize such results. Producing one, however, requires an interactive multi-step procedure in Cytoscape with the EnrichmentMap and AutoAnnotate applications. Here we present SimpNetworkViz, a Python package that generates enrichment maps from a gene list in a single command. Given differentially expressed genes, a ranked gene table, an enrichment table from another tool, or a list of pathway names, SimpNetworkViz performs over-representation analysis or preranked gene-set enrichment analysis, constructs the gene-set similarity network, identifies and labels clusters by Markov clustering, computes an overlap-free layout, and derives a matching gene-level network whose modules carry the same annotation. Multiple conditions, including the marker genes of single-cell clusters, are compared in one network through pie-chart nodes, and an entire gene-set library can be rendered with the enriched sets highlighted to place results in the context of the full pathway landscape. Outputs comprise vector figures (PDF, SVG), an interactive HTML viewer, tabular results and Cytoscape session files. Default parameters reproduce those of EnrichmentMap and AutoAnnotate, and the statistical procedures are verified against independent implementations. SimpNetworkViz is available at <https://github.com/helenhuangmath/SimpNetworkViz> under the MIT license.

**Contents**

- [1 Standard workflow](#1-standard-workflow)
  - [1.1 Quick start](#11-quick-start)
  - [1.2 How to get help](#12-how-to-get-help)
  - [1.3 Input data](#13-input-data)
  - [1.4 Running the analysis](#14-running-the-analysis)
  - [1.5 Exploring and exporting results](#15-exploring-and-exporting-results)
  - [1.6 Gallery of outputs](#16-gallery-of-outputs)
- [2 How to use in detail](#2-how-to-use-in-detail)
  - [2.1 Installation](#21-installation)
  - [2.2 Input formats](#22-input-formats)
  - [2.3 Enrichment test: ORA or GSEA](#23-enrichment-test-ora-or-gsea)
  - [2.4 Comparing several conditions](#24-comparing-several-conditions)
  - [2.5 Single-cell marker genes](#25-single-cell-marker-genes)
  - [2.6 The whole-library view](#26-the-whole-library-view)
  - [2.7 Grouping pathways](#27-grouping-pathways)
  - [2.8 Customizing the figures](#28-customizing-the-figures)
  - [2.9 Fine-tuning in Cytoscape](#29-fine-tuning-in-cytoscape)
  - [2.10 The point-and-click app](#210-the-point-and-click-app)
  - [2.11 Using the Python API](#211-using-the-python-api)
- [3 Theory behind SimpNetworkViz](#3-theory-behind-simpnetworkviz)
- [4 Reference](#4-reference)

---

## 1 Standard workflow

### 1.1 Quick start

Here we show the most basic steps for going from a list of differentially expressed genes to a pathway network figure. Install the package once:

```bash
git clone https://github.com/helenhuangmath/SimpNetworkViz
cd SimpNetworkViz
pip install -e .
```

Then, with a DE table `degs.csv` (gene symbol in the first column, a log2 fold-change column) and, ideally, a file `measured.txt` listing all genes measured in the experiment:

```bash
snv run degs.csv --library MSigDB_Hallmark_2020 --background measured.txt -o results
```

That is it. `results/pathway_network.pdf` is the figure, `results/pathway_network.html` the interactive version, `results/enrichment.csv` the enrichment table. To try the package without your own data, run the bundled example:

```bash
snv demo
```

The rest of section 1 explains each step; section 2 covers every option.

### 1.2 How to get help

* `snv -h` lists the commands, `snv run -h` every option of a command.
* Error messages say what to change (for example, relax `--fdr` when no pathway passes the cutoff).
* Bugs and feature requests: open an issue at <https://github.com/helenhuangmath/SimpNetworkViz/issues>.

If `snv` is not found after installing, use `python -m simpnetworkviz.cli` wherever this guide says `snv`.

### 1.3 Input data

SimpNetworkViz starts from one of four things. The standard workflow uses the first.

**A gene list.** A text file with one gene symbol per line, or (better) a CSV/TSV table with the gene in the first column and a log2 fold-change column, for example a DESeq2, edgeR or limma result filtered to the significant genes:

```
gene,log2FoldChange
MCM2,1.34
COL1A1,-2.80
IFIT1,2.95
...
```

With fold changes, up- and down-regulated genes are tested separately and the figure is colored red/blue. Column names such as `log2FoldChange`, `logFC` and `avg_log2FC` are recognized automatically. Gene symbols are matched regardless of case.

**A background.** A text file listing every gene that was measured (for RNA-seq, every gene that passed your expression filter). This is the universe for the statistical test. Without it, every gene in the library is used, which overstates significance. Pass it with `--background`.

**A gene-set library.** Pathways to test against. `snv libraries` lists the common choices:

| Library | When to use |
|---|---|
| `MSigDB_Hallmark_2020` | 50 broad, well-curated processes. Start here. |
| `KEGG_2021_Human`, `Reactome_2022`, `WikiPathways_2023_Human` | Signalling and metabolic pathways |
| `GO_Biological_Process_2023` | Fine-grained processes; large and redundant, so the clustering helps most |
| `KEGG_2019_Mouse`, `WikiPathways_2019_Mouse` | Mouse pathways |

Any [Enrichr library name](https://maayanlab.cloud/Enrichr/#libraries) works; it is downloaded once and cached in `~/.cache/simpnetworkviz`. Local `.gmt` files, e.g. from [MSigDB](https://www.gsea-msigdb.org/gsea/downloads.jsp), work with `--gmt file.gmt`.

The other starting points, a ranked table for GSEA, an enrichment table from another tool, a pathway list, and marker tables from single-cell data, are described in section 2.2.

### 1.4 Running the analysis

```bash
snv run degs.csv --library MSigDB_Hallmark_2020 --background measured.txt -o results
```

The command runs the enrichment, keeps pathways with FDR ≤ 0.05 (at most 100), connects pathways whose gene sets overlap, clusters and labels them, lays out the figure, and builds the gene networks. The terminal reports what it found:

```
Pathway network: 17 pathways, 16 edges, 5 clusters
  [ 1] Epithelial mesenchymal transition  (5 pathways)
  [ 2] Mitotic checkpoint  (3 pathways)
  [ 3] E2F targets  (3 pathways)
  [ 4] Interferon response  (3 pathways)
  [ 5] Oxidative phosphorylation  (2 pathways)
Gene network: 88 genes, 893 edges

Results written to .../results/
```

The three settings you will adjust most often:

```bash
--fdr 0.1             # relax the FDR cutoff if few pathways pass (default 0.05)
--max-pathways 60     # fewer nodes for a cleaner figure (default 100)
--cutoff 0.5          # stricter similarity for edges = smaller, tighter clusters (default 0.375)
```

A run overwrites its output folder, so use `-o` to keep versions apart.

### 1.5 Exploring and exporting results

The output folder contains:

| File | What it is |
|---|---|
| `pathway_network.pdf` / `.svg` | The figure, as vector graphics with editable text. Add `--formats png,pdf,svg,...` if you also want a PNG. |
| `pathway_network.html` | Interactive viewer: open in a browser |
| `gene_network.*`, `gene_pathway_network.*` | Gene-level figures |
| `cluster_barplot.pdf` / `.svg` | Cluster themes ranked by significance (a companion panel) |
| `enrichment.csv` | Enrichment table: pathway, overlap, set size, fold enrichment, p-value, FDR, direction, genes |
| `pathway_network_clusters.csv` | Cluster number, automatic label, size, best FDR, member pathways |
| `pathway_network_nodes.csv`, `_edges.csv` | Everything behind the figure, for your own plots |
| `pathway_network.cyjs` + `_style.xml` | For fine-tuning in Cytoscape (section 2.9) |
| `pathway_network.graphml` | For igraph, Gephi and other network tools |
| `run_params.json` | Every setting of this run, the package version and the date |

A sensible order to look at them:

1. `enrichment.csv`: do the top pathways make sense?
2. `pathway_network.pdf`: each circle is a pathway, each shaded outline a cluster of related pathways with an automatic label.
3. `pathway_network_clusters.csv`: which pathways make up each cluster.
4. `pathway_network.html`: click a node for its genes and FDR, click a cluster for its members, drag things around, search for a gene, export a PNG.
5. If a cluster label is awkward, rename it (section 2.8) or try another grouping (section 2.7).

**How to read the pathway network.** Node = pathway, size = number of genes in the gene set. Color: red = up-regulated, blue = down-regulated, darker = more significant. Edge = the two pathways share genes; width = similarity. Shaded outline + label = a cluster of related pathways. A node without an outline belongs to no cluster.

**How to read the gene network.** Node = one of your genes that hits a pathway in the network; size = number of pathways it belongs to; color = your log2FC. Edge = the two genes share a pathway (or interact in STRING, with `--gene-network string`). Outlines and labels are the same as in the pathway network, so each gene module matches its pathway theme.

### 1.6 Gallery of outputs

All figures come from `snv demo`, which runs offline on 88 genes with log2 fold changes and 32 shortened gene sets. The package writes vector PDF/SVG files; the previews below are PNG renderings of them. Click a figure for full size.

Pathway network, colored by direction (default) or by cluster:

| `--color-by direction` | `--color-by cluster` |
|---|---|
| [![](docs/images/pathway_network.png)](docs/images/pathway_network.png) | [![](docs/images/pathway_network_by_cluster.png)](docs/images/pathway_network_by_cluster.png) |

Gene network and gene–pathway ("cnet") network:

| Gene network | Gene–pathway network |
|---|---|
| [![](docs/images/gene_network.png)](docs/images/gene_network.png) | [![](docs/images/gene_pathway_network.png)](docs/images/gene_pathway_network.png) |

Pathway and gene network side by side (`plot_pathway_gene_network` in Python), sharing cluster colors and labels:

[![](docs/images/pathway_gene_panel.png)](docs/images/pathway_gene_panel.png)

Several conditions in one network (section 2.4). Each pathway is a pie chart with one slice per condition; a white slice means not enriched there. Left: slice color = condition, drawn on the whole library so non-enriched pathways appear grey. Right: slices colored red/blue by direction.

| Condition colors, whole library | Direction colors |
|---|---|
| [![](docs/images/compare_library.png)](docs/images/compare_library.png) | [![](docs/images/compare_direction.png)](docs/images/compare_direction.png) |

GSEA (section 2.3) colored by normalized enrichment score, and grouping by source database (section 2.7):

| GSEA, colored by NES | Grouped by source |
|---|---|
| [![](docs/images/gsea_network.png)](docs/images/gsea_network.png) | [![](docs/images/grouped_by_source.png)](docs/images/grouped_by_source.png) |

Cluster bar plot, ranking the themes by the FDR of their best pathway:

[![](docs/images/cluster_barplot.png)](docs/images/cluster_barplot.png)

The pathway map (section 2.6): all 5,489 GO Biological Process terms embedded once from gene-set similarity, the enriched terms of a Cl13 vs Arm comparison (from `examples/Cl13_vs_Arm`) drawn as pies on top. The same 35 labelled regions appear in both; only the 2-D layout differs between t-SNE and UMAP.

| t-SNE (`--method tsne`) | UMAP (`--method umap`) |
|---|---|
| [![](docs/images/pathway_map_tsne.png)](docs/images/pathway_map_tsne.png) | [![](docs/images/pathway_map_umap.png)](docs/images/pathway_map_umap.png) |

The same two maps without hits (the output of `snv map`), and with the six Cl13 subgroups as pie slices:

| | t-SNE | UMAP |
|---|---|---|
| Library only | [![](docs/images/pathway_map_tsne_library.png)](docs/images/pathway_map_tsne_library.png) | [![](docs/images/pathway_map_umap_library.png)](docs/images/pathway_map_umap_library.png) |
| Six Cl13 subgroups | [![](docs/images/pathway_map_tsne_within_Cl13.png)](docs/images/pathway_map_tsne_within_Cl13.png) | [![](docs/images/pathway_map_umap_within_Cl13.png)](docs/images/pathway_map_umap_within_Cl13.png) |

---

## 2 How to use in detail

### 2.1 Installation

Python 3.9 or newer is required.

```bash
git clone https://github.com/helenhuangmath/SimpNetworkViz
cd SimpNetworkViz
pip install -e .             # command line + Python API
pip install -e ".[gui]"      # optional: the point-and-click app (section 2.10)
pip install -e ".[map]"      # optional: UMAP/t-SNE for the pathway map (section 2.6)
pip install -e ".[test]"     # optional: run the tests with  pytest -q
snv --version
```

Dependencies are numpy, pandas, scipy, networkx, matplotlib and requests; they are installed automatically. Internet is needed only to download a library the first time and for STRING interactions. The demo runs offline.

### 2.2 Input formats

**Gene list** (`snv run FILE`): text with one gene per line, or a CSV/TSV with the gene in the first column. A numeric column named `log2FoldChange`, `logFC`, `avg_log2FC`, `log2FC`, `stat`, `t`, `score` or `value` is used as the fold change (the first numeric column otherwise). Header optional.

**Ranked table for GSEA** (`snv run FILE --method gsea`): the same format, but with *all* measured genes and a statistic for each (DESeq2 `stat`, limma `t`, or log2FC). Higher = up-regulated.

**Enrichment table** (`snv run --enrichment FILE`): results from Enrichr (`Term`, `Adjusted P-value`, `Genes`), clusterProfiler (`Description`, `p.adjust`, `geneID`), g:Profiler (`term_name`, `p_value`, `intersections`), GSEA desktop (`NAME`, `NES`, `FDR q-val`), Metascape-like exports, or this tool's own `enrichment.csv`. Columns are detected by name. Give the matching gene sets (`--library` or `--gmt`) so pathway overlaps can be computed from the full gene sets.

**Pathway list** (`snv run --pathways FILE`): one pathway name per line. Matching ignores case and punctuation, so `G2M checkpoint` finds `HALLMARK_G2M_CHECKPOINT`.

**Several gene lists** (`snv compare NAME=FILE ...`): one file per condition, each in the gene-list format; or one wide table with a column per group (`snv compare --groups FILE`).

**Marker table** (`snv markers FILE`): Seurat `FindAllMarkers` output (`cluster`, `gene`, `avg_log2FC`, `p_val_adj`; gene names may be in the row names), scanpy `sc.get.rank_genes_groups_df(adata, None)` (`group`, `names`, `logfoldchanges`, `pvals_adj`), scanpy's wide export (`<group>_n`, `<group>_l`, `<group>_p`), or any table with cluster, gene, log fold-change and adjusted p-value columns.

**Gene sets**: Enrichr library names (several at once, comma-separated: `--library KEGG_2021_Human,Reactome_2022`), local `.gmt` files (`--gmt`, repeatable), or `demo`.

### 2.3 Enrichment test: ORA or GSEA

**Over-representation analysis (default).** A hypergeometric test asks whether your gene list contains more genes of a pathway than expected by chance, given the background. P-values are corrected with Benjamini–Hochberg across all gene sets of the library. With fold changes, up and down genes are tested separately (`--no-split` to test them together). Gene sets with fewer than 5 or more than 1000 genes are skipped (`--min-size`, `--max-size`).

```bash
snv run degs.csv --library MSigDB_Hallmark_2020 --background measured.txt
```

**Preranked GSEA.** If you would rather not choose a fold-change or p-value cutoff, give the full ranked table:

```bash
snv run deseq2_all_genes.csv --method gsea --library MSigDB_Hallmark_2020 --fdr 0.25
```

Gene sets of 15–500 genes are tested with the weighted Kolmogorov–Smirnov statistic and 1000 gene permutations (`--permutations`). The table gains `es` and `nes` columns, the figure is colored by NES, and the gene network shows the leading-edge genes. FDR 0.25 is the customary GSEA cutoff.

Enrichment only, no figures: `snv enrich degs.csv --library ... -o enrichment.csv`.

### 2.4 Comparing several conditions

Knock-out vs. wild type, several time points, treatments or tissues: give each a name and a file.

```bash
snv compare KO=ko_degs.csv WT=wt_degs.csv Day7=day7_degs.csv --library Reactome_2022
```

If the groups are columns of one table (header = group name, columns may differ in length, as in a table of gene modules), use `--groups`:

```bash
snv compare --groups gene_modules.csv --library MSigDB_Hallmark_2020
```

Each condition is tested on its own. A pathway enters the network when it passes the FDR cutoff in at least one condition, and the node becomes a **pie chart** with one slice per condition: slice 1 at 12 o'clock, then clockwise, in the order you gave. A colored slice means enriched in that condition (darker = more significant); a white slice means not enriched. `--color-by direction` colors the slices red/blue by direction instead. The legend shows the slice order.

Additional outputs: `enrichment.csv` has a `condition` column; `enrichment_matrix.csv` is a pathway × condition table of signed −log10 FDR (for heatmaps); `pathway_network_nodes.csv` has `fdr__<condition>` and `direction__<condition>` columns; the HTML viewer lists the per-condition results when you click a node. In the gene networks, genes are colored by their mean log2FC across conditions, and the viewer shows which conditions each gene is a hit in.

`--method gsea` works here too when every file is a full ranked table.

A worked example on real data, six gene subgroups of exhausted (LCMV clone 13) versus memory (LCMV Armstrong) CD8 T cells, is in [examples/Cl13_vs_Arm/TUTORIAL.md](examples/Cl13_vs_Arm/TUTORIAL.md).

### 2.5 Single-cell marker genes

Export the markers from Seurat (`write.csv(FindAllMarkers(obj), "markers.csv")`) or scanpy (`sc.get.rank_genes_groups_df(adata, None).to_csv("markers.csv")`), then:

```bash
snv markers markers.csv --library MSigDB_Hallmark_2020 --whole-library
```

Each cell cluster becomes one condition (one pie slice), so the result is the comparison network of section 2.4. By default positive markers with adjusted p ≤ 0.05 and log2FC ≥ 0.25 are kept, at most 200 per cluster, ranked by the test score:

```bash
--padj 0.01 --min-logfc 0.5 --top-genes 100     # stricter markers
--direction both                                # keep negative markers too; up/down tested separately
--clusters "CD8 T,NK,B cells"                   # only these clusters
```

Because marker lists are short and the universe is the whole library, Hallmark, KEGG or Reactome give more readable results than GO.

### 2.6 The whole-library view

`--whole-library` draws **every pathway of the library** as a node, clusters the library itself, and highlights the pathways enriched in your data; the rest are grey, and clusters without any enriched pathway get a grey outline. This shows which part of the pathway landscape your genes, or each single-cell cluster, fall into, and makes figures from different experiments directly comparable because the layout describes the library, not the hits.

```bash
snv run degs.csv --library MSigDB_Hallmark_2020 --whole-library
snv markers markers.csv --library KEGG_2021_Human --whole-library
```

Libraries with more than 400 gene sets are trimmed to the enriched pathways plus the largest others; GO is too large to draw whole as a network. For large libraries use the pathway map below instead.

**The pathway map: a fixed coordinate system for a whole library.** `snv map` embeds every gene set of a library in two dimensions (UMAP, or t-SNE when `umap-learn` is not installed) from gene-set similarity, divides the map into labelled regions, and saves it as a CSV. Every later analysis can be drawn on that same map, so results from different experiments, conditions or cell clusters are directly comparable: the enriched pathways are colored points (or pies, for comparisons) on a grey background of the whole library, and `pathway_map_regions.csv` counts how many hits fall into each region.

```bash
pip install "simpnetworkviz[map]"                                   # umap-learn + scikit-learn, once
snv map --gmt m5.go.bp.v2025.1.Mm.symbols.gmt -o gobp_map.csv       # once per library (writes gobp_map.pdf too)
snv run genes.txt --gmt m5.go.bp.v2025.1.Mm.symbols.gmt --map gobp_map.csv
snv compare KO=ko.csv WT=wt.csv --gmt m5.go.bp.v2025.1.Mm.symbols.gmt --map gobp_map.csv
```

Each run adds `pathway_map.pdf`, `pathway_map_hits.csv` (every enriched pathway with its coordinates and region) and `pathway_map_regions.csv` (hits per region and condition). The map is built once and reused; keep the CSV with the library it was built from. `--method tsne|umap`, `--neighbors` and `--resolution` (more or fewer regions) tune the map. The regions do not depend on the embedding (they are Louvain communities in the 50-dimensional similarity space), so t-SNE and UMAP maps of the same library share their region labels and counts; UMAP keeps related regions closer together and leaves more white space, t-SNE spreads the regions more evenly (compare the two in the gallery, section 1.6). Both maps of the example are included under `examples/Cl13_vs_Arm`. In Python: `build_pathway_map`, `save_map`, `load_map`, `plot_pathway_map`, `map_hits`, `region_summary`.

### 2.7 Grouping pathways

The default grouping is Markov clustering (MCL) on the gene-overlap network, as in AutoAnnotate, with labels built from the words the member pathways share. Alternatives, chosen with `--cluster`:

| Method | What it does | Tuning |
|---|---|---|
| `mcl` (default) | Markov clustering on the similarity graph | `--inflation` (higher = more, smaller clusters; default 2) |
| `louvain`, `greedy` | Community detection | `--resolution` (higher = more clusters; default 1) |
| `label_propagation` | Fast community detection | |
| `hierarchical` | Average-linkage clustering on gene-set distance; also groups pathways related through intermediates | `--linkage-threshold` (default 1 − cutoff) |
| `components` | Connected components of the graph | |
| `category` | By source database: Hallmark, KEGG, Reactome, GO ..., read from the name, or from `--categories FILE` (`pathway<TAB>category`, e.g. MSigDB sub-collections) | |
| `keyword` | By your own themes from `--themes FILE` | |
| `none` | No grouping | |

A themes file has one theme per line with its keywords:

```
Proliferation	cycle, mitotic, E2F, replication, checkpoint
Immune	interferon, cytokine, antigen
Matrix	matrix, collagen, adhesion
```

The first matching theme wins, so list specific themes before general ones. Pathways matching nothing stay unclustered.

Two related switches: `--similarity jaccard|overlap|combined` changes the similarity measure, and `--similarity-genes hits` computes it from your hit genes only instead of the full gene sets (useful when you care about which of *your* genes the pathways share).

### 2.8 Customizing the figures

* **Colors.** `--color-by direction` (red/blue, shade = FDR), `cluster` (one pastel color per cluster), `significance` (single scale), `nes` (GSEA), `condition` (comparisons). Gene networks are colored by log2FC when available.
* **Formats, size and text.** Figures are written as PDF and SVG (vector, editable text); add `png` to `--formats` for a raster copy, with `--dpi` setting its resolution (300 by default). Cluster labels 14 pt, node labels 9 pt, legend 11–13 pt, all Arial (Liberation Sans where Arial is missing); the layout leaves room for them. `--title "..."` adds a title. In Python: `plot(cluster_label_size=16, node_label_size=8, max_width=12)`.
* **Node labels.** Shown on pathway networks with ≤ 30 nodes; `plot(node_labels="top")` labels only the best pathway per cluster, `"all"` or `"none"` as named.
* **Cluster outlines.** Each cluster has its own color so pathway and gene figures match; `plot(hull_style="neutral")` gives grey outlines.
* **Rename a cluster.** In Python, `pnet.group_labels[3] = "Interferon signalling"` then `pnet.plot(...)`; or edit the text in the PDF/SVG.
* **Gene network.** `--gene-network shared` (default, offline), `string` (STRING interactions; `--species 10090` for mouse), `cnet` (gene–pathway only) or `none`. `--max-genes` limits the genes shown (default 150).
* **Reproducible layout.** The layout is deterministic for a given `--seed` (default 0).
* **Offline HTML.** `--cytoscape-js path/to/cytoscape.min.js` embeds the JavaScript library in the viewer.

### 2.9 Fine-tuning in Cytoscape

1. *File ▸ Import ▸ Network from File*: choose `pathway_network.cyjs`. Positions are kept.
2. *File ▸ Import ▸ Styles from File*: choose `pathway_network_style.xml`, then apply the **SimpNetworkViz** style.
3. Optionally run AutoAnnotate on the existing `cluster` column.

For comparisons, the `fdr__<condition>` columns are exported so charts can be added in Cytoscape.

### 2.10 The point-and-click app

```bash
pip install -e ".[gui]"
snv gui
```

A page opens in the browser. Choose the input type in the sidebar (gene list, compare conditions, single-cell markers, enrichment table, pathway list), paste or upload the data, choose the library and the options, and click **Build networks**. Tabs show the figures, the interactive viewer and the tables; the Download tab gives everything as a zip.

### 2.11 Using the Python API

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

---

## 3 Theory behind SimpNetworkViz

**Enrichment.** Over-representation uses the hypergeometric test: for a pathway with *K* genes in a universe of *N*, and a list of *n* genes with *k* in the pathway, the p-value is P(X ≥ k). Every gene set of the library that falls within the size limits counts toward the Benjamini–Hochberg correction. Preranked GSEA follows Subramanian et al. (2005): the weighted Kolmogorov–Smirnov running sum (weight 1) gives the enrichment score; the null distribution comes from gene permutations, shared across sets of the same size; NES = ES / mean |null ES| of the same sign; empirical p-values are corrected with BH as in `fgsea` (`fdr_method="gsea"` in Python gives the original NES-based q-value).

**Similarity.** Two pathways are compared by their full gene sets with the combined coefficient of EnrichmentMap, ½·Jaccard + ½·overlap coefficient, and connected when it reaches 0.375. The overlap coefficient lets a small pathway attach to a large one that contains it; Jaccard tempers that for very unequal sizes.

**Clustering and labels.** Markov clustering (inflation 2) on the weighted similarity graph, as in AutoAnnotate's default. Each cluster is labelled with the informative words shared by at least two of its pathways, ordered as they appear in the names, or with its most significant pathway when it has fewer than three members. Clusters with a single pathway are left unlabelled.

**Layout.** Each cluster is laid out separately (Kamada–Kawai, then nodes pushed apart so that circles and labels never overlap) and rotated so its long axis is horizontal. Clusters are then packed in rows, each with room for its label; unclustered nodes go into a grid at the end. This mirrors AutoAnnotate's "layout network to prevent cluster overlap".

**Gene networks.** Genes are those of your list that hit the selected pathways (the leading edge for GSEA; the full gene sets when you start from a pathway list without genes). Each gene is assigned to the pathway cluster it belongs to most, so gene modules carry the same labels and colors as the pathway network. Edges are shared pathways (weight = number shared) or STRING interactions (score ≥ 0.4).

**Comparisons.** Each condition is tested separately. The merged node keeps the best FDR, the union of hit genes, and per-condition FDR, direction and genes. The direction of the node is `up` or `down` when all conditions agree and `mixed` otherwise.

**Whole library.** All pathways of the library are added before similarity and clustering, so the clusters describe the library; enriched pathways are highlighted.

**Validation.** `tests/test_validation.py` checks every ORA p-value against a one-sided Fisher exact test (scipy) and against `gseapy.enrich`; BH FDR against `statsmodels`/R `p.adjust`; GSEA enrichment scores and leading-edge genes against `gseapy.prerank` (identical); and the similarity formula against the EnrichmentMap definition. Live tests download `MSigDB_Hallmark_2020` from Enrichr and query STRING for known interactions, and are skipped when offline. Run them with `pip install -e ".[test]"` and `pytest -q`.

---

## 4 Reference

### Commands

| Command | Purpose |
|---|---|
| `snv run FILE` | gene list / ranked table / `--enrichment` / `--pathways` → networks |
| `snv compare NAME=FILE ...` / `--groups TABLE` | several conditions → one network with pie nodes |
| `snv markers FILE` | single-cell marker table → one slice per cluster |
| `snv map` | build a fixed 2-D map of a whole library, to reuse with `--map` |
| `snv enrich FILE` | enrichment table only, no figures |
| `snv demo [--compare]` | the bundled example |
| `snv libraries` | list common libraries |
| `snv gui` | the point-and-click app |

### Options of `snv run`, `snv compare` and `snv markers`

| Option | Default | Meaning |
|---|---|---|
| `--library NAME[,NAME]`, `--gmt FILE` | | gene sets |
| `--background FILE` | all library genes | universe for the test |
| `--method ora\|gsea` | ora | hypergeometric test or preranked GSEA |
| `--no-split` | off | do not test up/down genes separately |
| `--min-size`, `--max-size` | 5 / 1000 (ORA), 15 / 500 (GSEA) | gene-set size limits |
| `--permutations` | 1000 | GSEA gene permutations |
| `--fdr` | 0.05 | FDR cutoff for pathways in the network |
| `--max-pathways` | 100 | most significant pathways kept |
| `--similarity combined\|jaccard\|overlap` | combined | gene-set similarity for edges |
| `--similarity-genes full\|hits` | full | full gene sets or only hit genes |
| `--cutoff` | 0.375 | minimum similarity for an edge |
| `--cluster METHOD` | mcl | grouping method (section 2.7) |
| `--inflation`, `--resolution`, `--linkage-threshold` | 2 / 1 / 1 − cutoff | tuning of mcl / louvain, greedy / hierarchical |
| `--categories FILE`, `--themes FILE` | | for `--cluster category` / `keyword` |
| `--whole-library` | off | draw all pathways, highlight enriched |
| `--map MAP.csv` | | also draw the hits on a library map from `snv map` |
| `--color-by` | auto | cluster, direction, significance, nes, condition |
| `--gene-network shared\|string\|cnet\|none` | shared | gene network type |
| `--max-genes` | 150 | genes in the gene network |
| `--species` | 9606 | NCBI taxon for STRING (10090 = mouse) |
| `--formats` | pdf,svg,html,cyjs,graphml,csv | outputs to write (add `png` for raster figures) |
| `--title`, `--dpi`, `--seed`, `-o` | | title, PNG resolution if requested (300), layout seed (0), output folder |
| `snv markers` only: `--padj`, `--min-logfc`, `--top-genes`, `--direction`, `--clusters` | 0.05 / 0.25 / 200 / up / all | marker filters |

### Output files

`pathway_network.{pdf,svg,html,cyjs,graphml}`, `pathway_network_style.xml`, `pathway_network_{nodes,edges,clusters}.csv`, `gene_network.*`, `gene_pathway_network.*`, `cluster_barplot.{pdf,svg}`, `enrichment.csv`, `enrichment_matrix.csv` (comparisons), `run_params.json`.

### Notes

* The bundled `demo` gene sets are shortened versions of real pathways, meant only for testing and demos. Use full libraries for real analyses.
* Figures use Arial at regular weight (Liberation Sans, with identical metrics, where Arial is not installed). PDF text is embedded as editable type (fonttype 42); SVG text stays text.

License: MIT
