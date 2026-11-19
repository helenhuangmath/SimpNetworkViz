# markers

## Overview

`snv markers` reads a single-cell marker table, makes each cell cluster one condition, and draws the comparison network.

Export the markers from Seurat (`write.csv(FindAllMarkers(obj), "markers.csv")`) or scanpy (`sc.get.rank_genes_groups_df(adata, None).to_csv("markers.csv")`), then:

```bash
snv markers markers.csv --library MSigDB_Hallmark_2020 --whole-library
```

Each cell cluster becomes one condition (one pie slice), so the result is the comparison network of [Comparing several conditions](compare.md). By default positive markers with adjusted p ≤ 0.05 and log2FC ≥ 0.25 are kept, at most 200 per cluster, ranked by the test score:

```bash
--padj 0.01 --min-logfc 0.5 --top-genes 100     # stricter markers
--direction both                                # keep negative markers too; up/down tested separately
--clusters "CD8 T,NK,B cells"                   # only these clusters
```

Because marker lists are short and the universe is the whole library, Hallmark, KEGG or Reactome give more readable results than GO.

The accepted marker-table layouts are listed in [File formats](../fileformats.md#input-formats).

## Command-line options

```{eval-rst}
.. argparse::
   :module: _cli_for_docs
   :func: build_parser
   :prog: snv
   :path: markers
   :noepilog:
```

## Output files

The same files as [`snv compare`](compare.md#output-files).
