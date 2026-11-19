# libraries

## Overview

`snv libraries` prints a short list of widely used gene-set libraries. It takes no options and does not use the internet.

```bash
snv libraries
```

```text
Popular Enrichr libraries (any Enrichr library name works with --library):

  MSigDB_Hallmark_2020           MSigDB Hallmark gene sets (50 broad processes)
  KEGG_2021_Human                KEGG pathways (human)
  Reactome_2022                  Reactome pathways (human)
  WikiPathways_2023_Human        WikiPathways (human)
  GO_Biological_Process_2023     Gene Ontology biological process
  GO_Molecular_Function_2023     Gene Ontology molecular function
  GO_Cellular_Component_2023     Gene Ontology cellular component
  KEGG_2019_Mouse                KEGG pathways (mouse)
  WikiPathways_2019_Mouse        WikiPathways (mouse)
  BioPlanet_2019                 NCATS BioPlanet pathways
```

## Using a library

Pass a name to `--library` in [`run`](run.md) or [`enrich`](enrich.md). Several names can be joined with commas:

```bash
snv run genes.txt --library KEGG_2021_Human,Reactome_2022
```

The list is only a selection. Any name from the [Enrichr library page](https://maayanlab.cloud/Enrichr/#libraries) works. A library is downloaded the first time it is used and then cached in `~/.cache/simpnetworkviz`.

For gene sets that are not in Enrichr, or to pin an exact release, download a GMT file (for example from [MSigDB](https://www.gsea-msigdb.org/gsea/downloads.jsp)) and pass it with `--gmt`.

`--library demo` selects the bundled test set, which works offline.
