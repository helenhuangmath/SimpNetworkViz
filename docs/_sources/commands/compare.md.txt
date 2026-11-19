# compare

## Overview

`snv compare` tests several conditions separately and draws them in one network, where each pathway is a pie chart with one slice per condition.

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

A worked example on real data, six gene subgroups of exhausted (LCMV clone 13) versus memory (LCMV Armstrong) CD8 T cells, is in the [Cl13 vs Arm tutorial](../tutorial/cl13_vs_arm.md).

## Command-line options

```{eval-rst}
.. argparse::
   :module: _cli_for_docs
   :func: build_parser
   :prog: snv
   :path: compare
   :noepilog:
```

## Output files

The same files as [`snv run`](run.md#output-files), plus `enrichment_matrix.csv`.
