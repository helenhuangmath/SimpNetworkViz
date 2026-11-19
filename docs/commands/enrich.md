# enrich

## Overview

`snv enrich` runs only the enrichment test and writes the result table. Use it when you want to look at the table, or filter it, before drawing a network with [`snv run --enrichment`](run.md).

The test is over-representation analysis: a hypergeometric test for each gene set, with Benjamini–Hochberg FDR. When the gene list has fold changes, up- and down-regulated genes are tested separately and the `direction` column says which.

## Examples of command-line usage

1. A gene list against an Enrichr library:

   ```bash
   snv enrich genes.txt --library GO_Biological_Process_2023
   ```

2. A differential expression table against a local GMT, with a background:

   ```bash
   snv enrich degs.csv --gmt h.all.v2024.1.Hs.symbols.gmt --background measured_genes.txt -o hallmark.csv
   ```

3. Then draw the network from the table:

   ```bash
   snv run --enrichment hallmark.csv --gmt h.all.v2024.1.Hs.symbols.gmt
   ```

## Command-line options

* `GENES`

  Required. A gene list, in the same format as for [`run`](run.md#input).

* `--library NAME`, `--gmt FILE`

  The gene sets. One of the two is required. They work as in [`run`](run.md#input).

* `--background FILE`

  The gene universe. Without it, the universe is every gene in the gene sets.

* `--min-size INT`, `--max-size INT`

  Test only gene sets with this many genes, counted after restricting to the background. Defaults: 5 and 1000.

* `--no-split`

  Test all genes together even when fold changes are present.

* `-o, --output FILE`

  The output table. Default: `enrichment.csv`.

## Output

A CSV table with one row per tested pathway, sorted by significance. Pathways with fewer than two of your genes are not reported. The columns are described in [File formats](../fileformats.md#enrichmentcsv).

The command prints a one-line summary:

```text
20 pathways tested with overlap, 17 with FDR <= 0.05 -> enrichment.csv
```
