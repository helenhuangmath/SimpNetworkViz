# FAQ

## Installing and running

### The shell says `snv: command not found`

Use `python -m simpnetworkviz.cli` wherever this documentation says `snv`, or activate the environment you installed into.

### Which steps need the internet?

Internet is needed only to download a gene-set library the first time it is used, and for STRING interactions (`--gene-network string`). The demo runs offline. The interactive HTML files load a JavaScript library from a CDN when opened: Cytoscape.js for the network viewer and plotly.js for the pathway map.

### How do I make the network viewer work offline?

Pass a local copy of the library, and it is embedded in the HTML file:

```bash
snv run degs.csv --library MSigDB_Hallmark_2020 --cytoscape-js path/to/cytoscape.min.js
```

## Input

### Which gene-set library should I use?

Start with Hallmark: 50 broad, well-curated processes. KEGG, Reactome and WikiPathways cover signalling and metabolic pathways. GO biological process is fine-grained, large and redundant, so the clustering helps most there. [`snv libraries`](commands/libraries.md) lists the choices.

### My genes are mouse symbols. What should I use?

Use MSigDB's native mouse collections with `--msigdb NAME --species mouse`, for example `--msigdb Hallmark --msigdb GO_BP --species mouse`. Enrichr libraries are mostly human. Human libraries still match mouse symbols, because symbols are compared regardless of case, but genes are lost: in the [Cl13 vs Arm tutorial](tutorial/cl13_vs_arm.md), 22 % and 34 % of the two gene lists were lost against human collections.

### What should the background be?

Every gene that was measured: for RNA-seq, every gene that passed your expression filter. Without `--background`, every gene in the library is used, which overstates significance.

### Should I use ORA or GSEA?

Over-representation analysis, the default, needs a list of significant genes. Preranked GSEA (`--method gsea`) takes all measured genes with a statistic and needs no fold-change or p-value cutoff. See [Enrichment test: ORA or GSEA](commands/run.md#enrichment-test-ora-or-gsea).

### Can I use enrichment results from another tool?

Yes. Pass the table with `snv run --enrichment FILE`, and the matching gene sets with `--library` or `--gmt` so that pathway overlaps can be computed from the full gene sets. The recognized tables are listed in [File formats](fileformats.md).

## Results

### No pathway passes the cutoff

Relax the FDR cutoff with `--fdr 0.1`, and check that the gene symbols match the species of the library. For small gene groups, look at `enrichment.csv` rather than the network.

### The figure has too many pathways, or the clusters are too large

`--max-pathways 60` gives fewer nodes. `--cutoff 0.5` asks for stricter similarity, which gives smaller, tighter clusters. Other ways to group the pathways are in [Grouping pathways](guide/grouping.md).

### A cluster label is awkward

Labels are automatic, built from the words the member pathways share. Rename a cluster in Python with `pnet.group_labels[3] = "Interferon signalling"` and plot again, or edit the text in the PDF or SVG.

### Small, odd GO terms appear

Terms such as photoreceptor development appear when three or four broadly expressed genes coincide. `--min-size 20` removes most of them, and the cluster size and `best_fdr` columns make them easy to filter.

### How do I read a pie-chart node?

Each slice is one condition, starting at 12 o'clock and going clockwise in the order you gave. A colored slice means the pathway is enriched in that condition, and a white slice means it is not. The shade is significance, not effect size; fold enrichment is in `enrichment.csv`.

### GO is too large to draw as a whole library

Libraries with more than 400 gene sets are trimmed in the whole-library network. Use the [pathway map](guide/library_map.md) for large libraries.

## Output

### How do I get PNG figures?

Figures are written as PDF and SVG by default. Add `png` to `--formats`, for example `--formats png,pdf,html,csv`, and set the resolution with `--dpi`.

### The font in my PDF is not Arial

Where Arial is not installed, Liberation Sans is used. It has identical metrics, so the layout is the same. PDF text is embedded as editable type, so you can change the font afterwards.

### How do I cite SimpNetworkViz?

See [Citation](documentation.md#citation).
