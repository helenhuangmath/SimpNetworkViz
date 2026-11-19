# Common Q & A

## Installing and running

### The shell says `snv: command not found`

The environment you installed into is not active. Activate it, or run the tool as a module:

```bash
python -m simpnetworkviz.cli demo
```

### `snv gui` says it needs streamlit

The web app is an optional extra. Install it with:

```bash
pip install -e ".[gui]"
```

### Which steps need the internet?

Most of the package runs locally. These steps contact outside services:

| Step | Service | What is sent |
|---|---|---|
| Loading a gene-set library by name | Enrichr | The library name only. Your genes are not sent. |
| `--gene-network string` | STRING | The genes in the network |
| `snv.enrichr()` in Python | Enrichr | Your gene list |
| Opening the interactive HTML | A CDN for Cytoscape.js | Nothing from your data |

To work fully offline, use local `.gmt` files, keep the default `shared` gene network and embed Cytoscape.js. The default enrichment test is always computed on your machine.

### A library fails to download

Check the internet connection and the spelling of the name against the [Enrichr library page](https://maayanlab.cloud/Enrichr/#libraries). On a cluster without internet access, download a `.gmt` file on another machine and pass it with `--gmt`.

## Input

### Can I use mouse genes?

Yes. Gene symbols are matched without regard to case, so mouse symbols such as `Mki67` match human gene sets where the ortholog has the same name. Genes whose mouse and human symbols differ are lost. To avoid that, use a mouse library (`KEGG_2019_Mouse`, `WikiPathways_2019_Mouse`) or a mouse GMT file. For STRING, add `--species 10090`.

### Can I use Ensembl or Entrez IDs?

Only if the gene sets use the same kind of ID. SimpNetworkViz does not convert identifiers. Enrichr libraries use gene symbols, so convert your IDs to symbols first, or use a GMT file with matching IDs.

### Can I give it a full DESeq2 table?

Not as it is. Every gene in the file is treated as a hit, so filter the table to your significant genes first. See [Step 1 of the tutorial](tutorial/your_data.md#step-1-prepare-the-gene-list).

### What should the background be?

Every gene that could have been a hit: for RNA-seq, the genes that were expressed and tested. Without `--background`, all genes in the library are used, which usually makes p-values too small.

### Does it run GSEA?

No. The built-in test is over-representation analysis on a gene list. To draw a network from a ranked-list GSEA, run GSEA elsewhere and pass its result table with `--enrichment`.

## Results

### No pathways are drawn

Nothing passed the FDR cutoff. In order of likelihood:

1. The gene symbols do not match the library: wrong species, or IDs where symbols are expected.
2. The gene list is small, or the library is very large, so nothing survives correction. Try a smaller library such as `MSigDB_Hallmark_2020`.
3. The signal is weak. `--fdr 0.1` relaxes the cutoff; say so when you report the figure.

`snv enrich` writes the full table, so you can see how close the top pathways are.

### Most nodes are single and unconnected

The pathways share few genes at the current cutoff. Lower it, for example `--cutoff 0.3`. Libraries with a hierarchy, such as GO, give more connected networks than Hallmark, which was built to be non-redundant.

### Everything is one giant cluster

Raise the cutoff, for example `--cutoff 0.5`, or draw fewer pathways with `--max-pathways`.

### The figure is too crowded

Reduce `--max-pathways` or `--max-genes`, label fewer nodes with `node_labels="top"` in Python, or open the HTML viewer or the `.cyjs` file and rearrange by hand.

### A cluster label is not what I would call it

The label is built from words in the pathway names. Rename it in Python with `pnet.group_labels[3] = "My label"` and plot again. See [Rename the clusters](tutorial/figures.md#rename-the-clusters).

### Why do the FDR values differ from Enrichr or clusterProfiler?

The tools differ in the background, in the gene-set size limits and in which gene sets count toward the correction. SimpNetworkViz also tests up- and down-regulated genes separately when fold changes are given. To draw the network for another tool's numbers exactly, pass that tool's table with `--enrichment`.

### Why do two runs give slightly different figures?

The pathways, edges and cluster memberships are the same on every run. Two things can differ between runs:

* the label of a cluster whose top pathways have the same FDR
* the exact positions of nodes

`--seed` fixes the random part of the layout, but Python also randomizes the order of some internal collections each time it starts. To get identical files on every run, fix that as well:

```bash
PYTHONHASHSEED=0 snv run degs.csv --library MSigDB_Hallmark_2020
```

For a Python script, set the variable when you start Python: `PYTHONHASHSEED=0 python my_script.py`.

## Output

### How do I make the HTML viewer work offline?

Download `cytoscape.min.js` from the [Cytoscape.js releases](https://github.com/cytoscape/cytoscape.js/tree/master/dist) and embed it:

```bash
snv run degs.csv --library MSigDB_Hallmark_2020 --cytoscape-js cytoscape.min.js
```

### Why is there no SVG from `snv demo`?

`snv demo` writes PNG, PDF, HTML, Cytoscape, GraphML and CSV files. `snv run` also writes SVG by default. To get SVG for the demo data, run:

```bash
snv run "$(python -c 'import simpnetworkviz as s; print(s.genesets.demo_path("demo_genes.tsv"))')" \
    --library demo --cutoff 0.3 -o demo_svg
```

### The font in my PDF is not Arial

Arial is not installed on that machine, as on many Linux servers. Liberation Sans is used in its place. It has the same metrics, so the layout is the same, and you can switch the font to Arial when you edit the PDF.

### How do I cite it?

See [Citation](index.md#citation).
