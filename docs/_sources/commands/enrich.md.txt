# enrich

## Overview

`snv enrich` runs only the enrichment test, over-representation analysis or preranked GSEA, and writes the table. No figures are made. Use it to look at the table, or to filter it, before drawing a network with `snv run --enrichment`.

## Examples of command-line usage

```bash
snv enrich genes.txt --gmt h.all.v2024.1.Hs.symbols.gmt -o enrichment.csv
snv enrich deseq2_all_genes.csv --method gsea --library MSigDB_Hallmark_2020 -o gsea.csv
snv run --enrichment enrichment.csv --gmt h.all.v2024.1.Hs.symbols.gmt
```

The two tests are described under [`run`](run.md#enrichment-test-ora-or-gsea).

## Command-line options

```{eval-rst}
.. argparse::
   :module: _cli_for_docs
   :func: build_parser
   :prog: snv
   :path: enrich
   :noepilog:
```
