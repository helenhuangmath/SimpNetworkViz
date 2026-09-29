# SimpNetworkViz

SimpNetworkViz takes a gene list or a list of pathways and draws **EnrichmentMap / AutoAnnotate-style network figures**, the Cytoscape look you see in many papers (for example Fig. 6a of [Nat Cell Biol 2024](https://www.nature.com/articles/s41556-024-01570-0)). You don't need to install Cytoscape or click through it by hand.

* **Gene list → pathway enrichment.** Runs a hypergeometric test with Benjamini–Hochberg FDR. If you supply log2FC, up- and down-regulated genes are tested separately.
* **Pathways → clustered pathway network.** Nodes are pathways. Edges join pathways that share genes. Clusters come from MCL (as in AutoAnnotate) and get automatic word labels, with a shaded outline around each cluster.
* **Pathways → gene networks.** Three kinds:
  * genes linked by the pathways they share
  * STRING protein–protein interactions
  * a gene–pathway "cnet" plot
* **Output formats:**
  * publication-ready PNG/PDF/SVG (text stays editable in Illustrator)
  * an **interactive HTML** viewer (drag, search, click for gene lists)
  * `.cyjs` and a style file for **Cytoscape desktop**
  * GraphML and CSV tables

## Example results (bundled demo data)

Everything below comes from `snv demo`, which runs offline on the bundled demo data: 88 genes with log2FC (53 up, 35 down) and 32 demo gene sets. Figures are 300 dpi PNGs; click one to see it at full size.

### 1. Enrichment

Up- and down-regulated genes are tested separately (hypergeometric test, BH FDR). 20 of the 32 gene sets overlap the genes, and **17 of those 20 have FDR ≤ 0.05.** The top of `enrichment.csv`:

| Pathway | Direction | Overlap / set size | Fold enrichment | FDR |
|---|---|---|---|---|
| HALLMARK_G2M_CHECKPOINT | up | 23 / 30 | 7.0 | 1.8e-17 |
| HALLMARK_E2F_TARGETS | up | 24 / 34 | 6.5 | 3.2e-17 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | down | 18 / 28 | 8.9 | 5.2e-15 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | down | 16 / 28 | 7.9 | 2.3e-12 |
| REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION | down | 16 / 28 | 7.9 | 2.3e-12 |
| REACTOME_RESPIRATORY_ELECTRON_TRANSPORT | down | 12 / 25 | 6.7 | 6.0e-08 |
| GOBP_MITOTIC_SISTER_CHROMATID_SEGREGATION | up | 13 / 25 | 4.8 | 1.1e-06 |
| REACTOME_INTERFERON_SIGNALING | up | 14 / 29 | 4.4 | 1.1e-06 |
| KEGG_CELL_CYCLE | up | 14 / 30 | 4.3 | 1.6e-06 |
| REACTOME_MITOTIC_SPINDLE_CHECKPOINT | up | 11 / 20 | 5.0 | 3.6e-06 |

Not significant: HALLMARK_DNA_REPAIR (FDR 0.83), KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION (1.0), GOBP_DOUBLE_STRAND_BREAK_REPAIR (1.0).

### 2. Pathway network

**17 pathways, 16 similarity edges, 5 MCL clusters.** KEGG_CITRATE_CYCLE_TCA_CYCLE shares too few genes with the others to join a cluster, so it stays a lone node. Cluster labels are generated automatically (`pathway_network_clusters.csv`):

| # | Cluster label | Pathways | Best FDR | Members |
|---|---|---|---|---|
| 1 | Epithelial mesenchymal transition | 5 | 2.3e-12 | EMT, ECM organization, collagen formation, ECM–receptor interaction, focal adhesion (all down) |
| 2 | Mitotic checkpoint | 3 | 1.8e-17 | G2M checkpoint, mitotic sister chromatid segregation, mitotic spindle checkpoint (all up) |
| 3 | E2F targets | 3 | 3.2e-17 | E2F targets, cell cycle, DNA replication (all up) |
| 4 | Interferon response | 3 | 1.1e-06 | Interferon signaling, IFN-γ response, IFN-α response (all up) |
| 5 | Oxidative phosphorylation | 2 | 5.2e-15 | Oxidative phosphorylation, respiratory electron transport (all down) |

The data show proliferation and interferon programs going up, while ECM/EMT and mitochondrial respiration go down.

| Colored by direction (red = up, blue = down, shade = −log10 FDR) | Colored by cluster |
|---|---|
| [![](docs/images/pathway_network.png)](docs/images/pathway_network.png) | [![](docs/images/pathway_network_by_cluster.png)](docs/images/pathway_network_by_cluster.png) |

Cluster themes ranked by their best pathway (`cluster_barplot.png`):

[![](docs/images/cluster_barplot.png)](docs/images/cluster_barplot.png)

### 3. Gene networks

* **Gene network (shared pathways): 88 genes, 893 edges.** Each gene is assigned to the pathway cluster it belongs to most, which gives 17–18 genes per cluster. Nodes are colored by log2FC.
* **Gene–pathway (cnet) network: 103 nodes (genes + pathways), 209 gene–pathway edges.**

| Gene network (grouped by pathway cluster, colored by log2FC) | Gene–pathway network |
|---|---|
| [![](docs/images/gene_network.png)](docs/images/gene_network.png) | [![](docs/images/gene_pathway_network.png)](docs/images/gene_pathway_network.png) |

The pathway network and the gene network are saved as **separate figures**. Clusters keep the same outline color and label in both, so each gene module matches its pathway theme.

## Install

```bash
git clone https://github.com/helenhuangmath/SimpNetworkViz
cd SimpNetworkViz
pip install -e .            # command line + Python API
pip install -e ".[gui]"     # optional point-and-click web app
```

## Three ways to use it

### 1. Point-and-click app

```bash
snv gui
```

A browser page opens. Paste your genes, pick a library, and click **Build networks**. You can then download every figure and table as one zip.

### 2. Command line

```bash
snv demo                                                   # offline example -> demo_output/

# gene list -> enrichment -> pathway + gene networks
snv run genes.txt --library MSigDB_Hallmark_2020
snv run deseq2_results.csv --library GO_Biological_Process_2023 --background all_measured_genes.txt

# you already have enrichment results (Enrichr, clusterProfiler, g:Profiler, GSEA, ...)
snv run --enrichment clusterProfiler_enrichGO.csv --gmt c5.go.bp.v2024.1.Hs.symbols.gmt

# you have a list of pathways to draw
snv run --pathways my_pathways.txt --library Reactome_2022 --gene-network cnet

snv libraries                                              # list popular gene-set libraries
snv run -h                                                 # all options
```

Every run writes this to `snv_output/` (change it with `-o`):

| File | What it is |
|---|---|
| `pathway_network.png/.pdf/.svg` | Clustered pathway network figure |
| `pathway_network.html` | Interactive viewer: open it in a browser |
| `pathway_network.cyjs` + `_style.xml` | Import into Cytoscape desktop |
| `pathway_network_clusters.csv` | Cluster labels, sizes, best FDR and member pathways |
| `gene_network.*`, `gene_pathway_network.*` | Gene-level networks |
| `cluster_barplot.png` | Cluster themes ranked by significance (a companion panel) |
| `enrichment.csv` | Enrichment table (when you start from genes) |

### 3. Python

```python
import simpnetworkviz as snv

genes = snv.read_gene_list("degs.csv")                    # {gene: log2FC} or plain list
gene_sets = snv.load_gene_sets("MSigDB_Hallmark_2020")     # Enrichr name, .gmt path, or "demo"

table = snv.analyze_genes(genes, gene_sets)                # enrichment (up/down split)
pnet = snv.build_pathway_network(table, gene_sets)         # similarity + MCL + labels
pnet.plot("pathway_network.pdf", color_by="direction")     # or "cluster", "significance", "nes"
pnet.to_html("pathway_network.html")

pnet.group_labels[3] = "Interferon signalling"             # rename a cluster by hand
gnet = snv.build_gene_network(pnet, gene_sets, mode="shared", clusters=[1, 3], gene_values=genes)
gnet.plot("gene_network.pdf")
# optional: both networks side by side in one a|b figure
# snv.plot_pathway_gene_network(pnet, gnet, "pathway_gene_network.pdf")
```

`examples/quickstart.py` shows the full walkthrough.

## Inputs

* **Gene list.** A text file with one gene per line, or a CSV/TSV whose first column is the gene. If the file has a log2FC column (`log2FoldChange`, `logFC`, `avg_log2FC`, …), up- and down-regulated genes are tested separately and gene nodes are colored by fold change. Gene symbols are matched without regard to case.
* **Gene sets.**
  * any [Enrichr library](https://maayanlab.cloud/Enrichr/#libraries) by name (downloaded once, then cached in `~/.cache/simpnetworkviz`)
  * a local `.gmt` file, such as one from [MSigDB](https://www.gsea-msigdb.org/gsea/downloads.jsp)
  * several libraries at once, comma-separated: `--library KEGG_2021_Human,Reactome_2022`
* **Enrichment tables.** Columns are detected automatically. Supported formats:
  * Enrichr (`Term`, `Adjusted P-value`, `Genes`)
  * clusterProfiler (`Description`, `p.adjust`, `geneID`)
  * g:Profiler (`term_name`, `p_value`, `intersections`)
  * GSEA (`NAME`, `NES`, `FDR q-val`)
  * this tool's own output
* **Pathway list.** One name per line. Matching ignores case and punctuation, so `G2M checkpoint` matches `HALLMARK_G2M_CHECKPOINT`.

**Background genes:** for RNA-seq, pass all expressed or measured genes with `--background`. Without it, the universe is every gene in the library.

## How the network is built

The defaults follow EnrichmentMap and AutoAnnotate:

1. **Filter.** Keep pathways with FDR ≤ 0.05 (`--fdr`), up to 100 of them (`--max-pathways`).
2. **Similarity.** Compare the full gene sets using combined similarity = ½·Jaccard + ½·overlap coefficient. Draw an edge when similarity is ≥ 0.375 (`--cutoff`, `--similarity`).
3. **Clusters.** Cluster with Markov clustering, MCL (`--cluster mcl|louvain|components`). Each cluster is labeled with its most frequent informative words, or with its most significant pathway when the cluster is small.
4. **Layout.** Lay out each cluster on its own, then pack the clusters so they don't overlap and each label has room. This is AutoAnnotate's "layout to prevent cluster overlap".
5. **Styling.**
   * Node size = gene-set size.
   * Node color = direction and −log10 FDR, or cluster.
   * Edge width = similarity.

**Gene networks.** Genes come from the pathway hits (your genes). If you give only a pathway list and no genes, the full gene sets are used. Each gene is assigned to the pathway cluster it belongs to most, so gene modules carry the same labels as the pathway network. Use `--gene-network string --species 10090` for mouse STRING interactions.

## Tips for publication figures

* Default text sizes are cluster labels 14 pt, node labels 9 pt, legend 11–13 pt and title 18 pt. The layout leaves room for text of these sizes. For a one-off change, pass `cluster_label_size=` or `node_label_size=` to `plot()`. To change the defaults, edit the constants at the top of `layout.py` and `plotting.py`.
* PNGs are saved at 300 dpi. For sharper images use `--dpi 600` on the command line, or `plot(..., dpi=600)` in Python. PDF and SVG are vector files and look sharp at any zoom.
* Cluster outlines use the cluster colors by default so separate pathway and gene figures match; use `plot(hull_style="neutral")` for grey outlines.
* All figure text is Arial, regular weight (Liberation Sans, which has identical metrics, is used where Arial is not installed). The PDF/SVG text stays editable. Font sizes are chosen for a figure about 16 inches wide or smaller; set `max_width=` in `plot()` to change that.
* `--color-by cluster` gives the pastel "one color per theme" look. `direction` gives the classic red/blue EnrichmentMap look.
* To fine-tune in Cytoscape:
  1. Import the `.cyjs` file (*File ▸ Import ▸ Network from File*).
  2. Import `_style.xml` (*File ▸ Import ▸ Styles from File*) and apply the **SimpNetworkViz** style.
  3. Optionally run AutoAnnotate using the existing `cluster` column.
* `--cytoscape-js path/to/cytoscape.min.js` embeds the JavaScript library, so the HTML file also works offline.

## Development

```bash
pip install -e ".[test]"
pytest -q
```

The bundled `demo` gene sets (`simpnetworkviz/data/demo_genesets.gmt`) are **shortened versions of real pathways, meant only for testing and demos**. Use full libraries for real analyses.

License: MIT
