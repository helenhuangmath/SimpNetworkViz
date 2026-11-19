# Standard workflow

This page goes from a list of differentially expressed genes to a pathway network figure, one step at a time. Each subcommand and option is described in [Subcommands](commands/index.md) and the [User guide](guide/index.md).

## Quick start

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


## How to get help

* `snv -h` lists the commands, `snv run -h` every option of a command.
* Error messages say what to change (for example, relax `--fdr` when no pathway passes the cutoff).
* Bugs and feature requests: open an issue at <https://github.com/helenhuangmath/SimpNetworkViz/issues>.

If `snv` is not found after installing, use `python -m simpnetworkviz.cli` wherever this guide says `snv`.

## Input data

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

**Species.** Enrichr libraries are mostly human. For mouse data, or to use MSigDB directly, choose an MSigDB collection with `--msigdb NAME --species mouse|human` (downloaded once from the Broad release server, release 2025.1): `Hallmark`, `GO_BP`, `GO_CC`, `GO_MF`, `Reactome`, `WikiPathways`, `BioCarta`, `TFT`, `CellType`, `All` exist for both species; `KEGG`, `PID`, `HPO`, `Oncogenic`, `ImmuneSigDB` and `Vaccine` for human only (`snv libraries` prints the table). Mouse collections are MSigDB's native mouse sets, not orthology-mapped human ones. `--species` also sets the STRING taxon for `--gene-network string`.

```bash
snv run degs.csv --msigdb Hallmark --msigdb GO_BP --species mouse --background measured.txt -o results
```

The other starting points, a ranked table for GSEA, an enrichment table from another tool, a pathway list, and marker tables from single-cell data, are described in [File formats](fileformats.md).

## Running the analysis

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

## Exploring and exporting results

The output folder contains:

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

A sensible order to look at them:

1. `enrichment.csv`: do the top pathways make sense?
2. `pathway_network.pdf`: each circle is a pathway, each shaded outline a cluster of related pathways with an automatic label.
3. `pathway_network_clusters.csv`: which pathways make up each cluster.
4. `pathway_network.html`: click a node for its genes and FDR, click a cluster for its members, drag things around, search for a gene, export a PNG.
5. If a cluster label is awkward, rename it ([Customizing the figures](guide/figures.md)) or try another grouping ([Grouping pathways](guide/grouping.md)).

**How to read the pathway network.** Node = pathway, size = number of genes in the gene set. Color: red = up-regulated, blue = down-regulated, darker = more significant. Edge = the two pathways share genes; width = similarity. Shaded outline + label = a cluster of related pathways. A node without an outline belongs to no cluster.

**How to read the gene network.** Node = one of your genes that hits a pathway in the network; size = number of pathways it belongs to; color = your log2FC. Edge = the two genes share a pathway (or interact in STRING, with `--gene-network string`). Outlines and labels are the same as in the pathway network, so each gene module matches its pathway theme.
