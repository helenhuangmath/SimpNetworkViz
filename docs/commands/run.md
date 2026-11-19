# run

## Overview

`snv run` is the main command. It takes one of three starting points, builds a clustered pathway network and a gene network, and writes the figures and tables to an output folder.

| You have | You pass | What `run` does first |
|---|---|---|
| A gene list, with or without log2FC | `GENES` | Runs the enrichment test |
| An enrichment table from another tool | `--enrichment TABLE` | Reads the table and skips the test |
| A list of pathway names | `--pathways FILE` | Looks the pathways up in the gene sets and draws all of them |

## Examples of command-line usage

1. Gene list, one Enrichr library:

   ```bash
   snv run genes.txt --library MSigDB_Hallmark_2020
   ```

2. Differential expression table with a log2FC column, two libraries and a background:

   ```bash
   snv run degs.csv --library KEGG_2021_Human,Reactome_2022 --background measured_genes.txt
   ```

3. A local MSigDB GMT file in place of an Enrichr library:

   ```bash
   snv run degs.csv --gmt h.all.v2024.1.Hs.symbols.gmt
   ```

4. An enrichment table from clusterProfiler, with the GMT it was computed from:

   ```bash
   snv run --enrichment clusterProfiler_enrichGO.csv --gmt c5.go.bp.v2024.1.Hs.symbols.gmt
   ```

5. A list of pathway names, with a gene–pathway plot as the gene network:

   ```bash
   snv run --pathways my_pathways.txt --library Reactome_2022 --gene-network cnet
   ```

6. A stricter network, colored by cluster, PDF only:

   ```bash
   snv run degs.csv --library GO_Biological_Process_2023 --fdr 0.01 --max-pathways 60 \
       --color-by cluster --formats pdf,csv -o go_bp_network
   ```

7. Mouse genes with STRING interactions as the gene network:

   ```bash
   snv run degs_mouse.csv --library KEGG_2019_Mouse --gene-network string --species 10090
   ```

## Essential command-line options

### Input

* `GENES`

  A gene list: a text file with one gene per line, or a CSV/TSV whose first column is the gene. If the table has a numeric column such as `log2FoldChange`, up- and down-regulated genes are tested separately and gene nodes are colored by fold change. See [File formats](../fileformats.md#gene-list).

  `GENES` can also be given together with `--enrichment` or `--pathways`. In that case no test is run, and the genes are used to choose and color the nodes of the gene network.

* `--enrichment TABLE`

  An existing enrichment table. Tables from this tool, Enrichr, clusterProfiler, g:Profiler and GSEA are recognized by their column names. See [File formats](../fileformats.md#enrichment-table).

* `--pathways FILE`

  A text file with one pathway name per line. Matching ignores case and punctuation, so `G2M checkpoint` matches `HALLMARK_G2M_CHECKPOINT`. All listed pathways are drawn; the FDR cutoff is not applied.

* `--library NAME`

  One or more [Enrichr library](https://maayanlab.cloud/Enrichr/#libraries) names, comma-separated, or `demo` for the bundled test set. A library is downloaded on first use and cached. Run [`snv libraries`](libraries.md) for a short list.

* `--gmt FILE`

  A local GMT gene-set file. Repeat the option to use several files. Files can be combined with `--library` or `--msigdb`.

  Gene sets are required for a gene list and for a pathway list. With `--enrichment` they are optional but recommended: with them, pathway similarity is computed from the full gene sets, as EnrichmentMap does; without them, only the genes listed in the table are used.

* `--background FILE`

  The gene universe for the enrichment test, in the same format as a gene list. For RNA-seq, pass all expressed or measured genes. Without it, the universe is every gene in the gene sets, which usually makes p-values too small.

* `--no-split`

  Test all genes together even when fold changes are present. By default, genes with positive and negative log2FC are tested separately.

### Options controlling the pathway network

* `--fdr FLOAT`

  Keep pathways with FDR at or below this value. Default: 0.05.

* `--max-pathways INT`

  Draw at most this many pathways, taking the most significant first. Default: 100.

* `--similarity {combined,jaccard,overlap}`

  How two gene sets are compared. `combined` is ½·Jaccard + ½·overlap coefficient, the EnrichmentMap default. Default: `combined`.

* `--cutoff FLOAT`

  Draw an edge when similarity is at or above this value. Lower it to connect more pathways, raise it to split large clusters. Default: 0.375.

* `--cluster {mcl,louvain,components,none}`

  The clustering method. `mcl` is Markov clustering, as in AutoAnnotate. `components` makes each connected group one cluster. `none` draws no clusters. Default: `mcl`.

* `--color-by {auto,cluster,direction,significance,nes}`

  How pathway nodes are colored. `direction` is red for up and blue for down, shaded by −log10 FDR. `cluster` gives one color per cluster. `significance` shades by −log10 FDR. `nes` uses the normalized enrichment score from a GSEA table. `auto` picks one based on the columns in the table. Default: `auto`.

### Options controlling the gene network

* `--gene-network {shared,string,cnet,none}`

  * `shared`: genes are linked when they belong to the same pathways. Runs offline.
  * `string`: protein–protein interactions from STRING. Needs internet and sends the gene names to STRING.
  * `cnet`: a gene–pathway plot, where genes are linked to the pathways they belong to.
  * `none`: no gene network.

  Default: `shared`. With `shared` or `string`, a `cnet` plot is also written as `gene_pathway_network.*`.

* `--max-genes INT`

  Draw at most this many genes. Default: 150.

* `--species SPECIES`

  Species for MSigDB and STRING. Use `human` (the default) or `mouse`. The taxon IDs 9606 and 10090 are also accepted.

### Output and other options

* `-o, --outdir DIR`

  The output folder. It is created if missing, and existing files with the same names are overwritten. Default: `snv_output`.

* `--formats LIST`

  Comma-separated output formats. Default: `png,pdf,svg,html,cyjs,graphml,csv`.

* `--title TEXT`

  A title for the pathway network figure.

* `--dpi INT`

  PNG resolution. Default: 300.

* `--cytoscape-js FILE`

  A local copy of `cytoscape.min.js` to embed in the HTML viewer, so the HTML works offline. Without it, the viewer loads the library from a CDN.

* `--seed INT`

  The random seed for layouts. Default: 0. See [Why do two runs give slightly different figures?](../qa.md#why-do-two-runs-give-slightly-different-figures)

## Output files

Each network is written with a prefix: `pathway_network`, `gene_network` or `gene_pathway_network`.

1. `PREFIX.png`, `PREFIX.pdf`, `PREFIX.svg`

   The figure. PDF and SVG are vector files, and their text stays editable.

2. `PREFIX.html`

   An interactive viewer. Open it in a browser to drag nodes, search, and click a node for its gene list.

3. `PREFIX.cyjs` and `PREFIX_style.xml`

   The network and its style for Cytoscape desktop.

4. `PREFIX.graphml`

   The network in GraphML, for other network tools.

5. `PREFIX_nodes.csv`, `PREFIX_edges.csv`

   The node and edge tables.

6. `pathway_network_clusters.csv`

   One row per cluster: its label, size, best FDR and member pathways.

7. `cluster_barplot.png`, `cluster_barplot.pdf`

   Cluster themes ranked by their most significant pathway.

8. `enrichment.csv`

   The enrichment table. Written only when you start from a gene list.

The columns of each table are described in [File formats](../fileformats.md#output-tables).

The command also prints a summary:

```text
Pathway network: 17 pathways, 16 edges, 5 clusters
  [ 1] Extracellular matrix organization  (5 pathways)
  [ 2] Mitotic checkpoint  (3 pathways)
  [ 3] E2F targets  (3 pathways)
  [ 4] Interferon response  (3 pathways)
  [ 5] Oxidative phosphorylation  (2 pathways)
Gene network: 88 genes, 893 edges
```
