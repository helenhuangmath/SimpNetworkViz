# run

## Overview

`snv run` is the main command. It takes one of four starting points, builds a clustered pathway network and the matching gene networks, and writes the figures and tables to an output folder.

| You have | You pass | What `run` does first |
|---|---|---|
| A gene list, with or without log2FC | `FILE` | Over-representation analysis |
| A ranked table of all measured genes | `FILE --method gsea` | Preranked GSEA |
| An enrichment table from another tool | `--enrichment FILE` | Reads the table and skips the test |
| A list of pathway names | `--pathways FILE` | Looks the pathways up in the gene sets and draws them |

The input formats are described in [File formats](../fileformats.md).

## Examples of command-line usage

1. A differential expression table, one Enrichr library and a background:

   ```bash
   snv run degs.csv --library MSigDB_Hallmark_2020 --background measured.txt -o results
   ```

2. Two libraries at once:

   ```bash
   snv run degs.csv --library KEGG_2021_Human,Reactome_2022 --background measured.txt
   ```

3. Mouse genes against MSigDB mouse collections:

   ```bash
   snv run degs.csv --msigdb Hallmark --msigdb GO_BP --species mouse --background measured.txt -o results
   ```

4. Preranked GSEA on all measured genes:

   ```bash
   snv run deseq2_all_genes.csv --method gsea --library MSigDB_Hallmark_2020 --fdr 0.25
   ```

5. An enrichment table from clusterProfiler, with the gene sets it was computed from:

   ```bash
   snv run --enrichment clusterprofiler_result.csv --gmt c2.cp.v2024.gmt
   ```

6. A list of pathway names:

   ```bash
   snv run --pathways my_pathways.txt --library GO_Biological_Process_2023
   ```

7. The whole library with the enriched pathways highlighted:

   ```bash
   snv run genes.txt --library MSigDB_Hallmark_2020 --whole-library
   ```

8. The hits drawn on a fixed library map built with [`snv map`](map.md):

   ```bash
   snv run genes.txt --gmt m5.go.bp.v2025.1.Mm.symbols.gmt --map gobp_map.csv
   ```

## Enrichment test: ORA or GSEA

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

## The settings you will adjust most often

```bash
```
```
```
```bash
--fdr 0.1             # relax the FDR cutoff if few pathways pass (default 0.05)
--max-pathways 60     # fewer nodes for a cleaner figure (default 100)
--cutoff 0.5          # stricter similarity for edges = smaller, tighter clusters (default 0.375)
```

A run overwrites its output folder, so use `-o` to keep versions apart.

## Command-line options

```{eval-rst}
.. argparse::
   :module: _cli_for_docs
   :func: build_parser
   :prog: snv
   :path: run
   :noepilog:
```

## Output files

`pathway_network.{pdf,svg,html,cyjs,graphml}`, `pathway_network_style.xml`, `pathway_network_{nodes,edges,clusters}.csv`, `gene_network.*`, `gene_pathway_network.*`, `cluster_barplot.{pdf,svg}`, `enrichment.csv`, `enrichment_matrix.csv` (comparisons), `run_params.json`.

Each file is described in [File formats](../fileformats.md#output-files).
