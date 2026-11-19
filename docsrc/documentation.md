# SimpNetworkViz documentation

SimpNetworkViz is an open-source Python package that turns gene lists, ranked scores or existing enrichment results into annotated pathway and gene networks, the EnrichmentMap / AutoAnnotate figures usually made by hand in Cytoscape. It also compares several conditions in one network and places their hits on a map of an entire gene-set library.

::::{grid} 1 2 2 3
:gutter: 3

:::{grid-item-card} Install
:link: install
:link-type: doc

Requirements, installation from GitHub and the optional extras.
:::

:::{grid-item-card} Get started
:link: tutorial/index
:link-type: doc

Step-by-step tutorials with code, output and figures.
:::

:::{grid-item-card} Subcommands
:link: commands/index
:link-type: doc

Every `snv` subcommand with its options.
:::

:::{grid-item-card} User guide
:link: guide/index
:link-type: doc

The pathway map, grouping, figures and the Python API.
:::

:::{grid-item-card} Gallery
:link: gallery
:link-type: doc

Example outputs of every kind.
:::

:::{grid-item-card} Reference
:link: api
:link-type: doc

All functions of the Python API.
:::
::::

## Quick start

Install from GitHub and run the bundled demo. It generates an enrichment table, network figures and an interactive HTML viewer.

```bash
git clone https://github.com/helenhuangmath/SimpNetworkViz.git
cd SimpNetworkViz
python -m pip install -e .
snv demo
```

Then follow the [guided tutorial](tutorial/demo.md), which explains the same analysis step by step.

## Usage

Example for a differential expression table with a log2 fold-change column:

```bash
snv run degs.csv --library MSigDB_Hallmark_2020 --background measured.txt -o results
```

Example for preranked GSEA on all measured genes:

```bash
snv run deseq2_all_genes.csv --method gsea --library MSigDB_Hallmark_2020 --fdr 0.25
```

Example for comparing several conditions in one network:

```bash
snv compare KO=ko_degs.csv WT=wt_degs.csv Day7=day7_degs.csv --library Reactome_2022
```

Example for single-cell marker genes:

```bash
snv markers markers.csv --library MSigDB_Hallmark_2020 --whole-library
```

Example for mouse genes on a fixed map of all GO biological process terms:

```bash
snv map --msigdb GO_BP --species mouse -o gobp_map.csv
snv run degs.csv --msigdb GO_BP --species mouse --map gobp_map.csv
```

These are the subcommands:

| Subcommand | Description |
|---|---|
| [`run`](commands/run.md) | Gene list, ranked table, enrichment table or pathway list in; pathway and gene networks out. |
| [`compare`](commands/compare.md) | Several conditions in one network, with pie-chart nodes. |
| [`markers`](commands/markers.md) | Single-cell marker table in; one pie slice per cell cluster. |
| [`map`](commands/map.md) | Build a fixed 2-D map of a whole gene-set library, to reuse with `--map`. |
| [`enrich`](commands/enrich.md) | Run only the enrichment test and write the table. |
| [`demo`](commands/demo.md) | Run the full workflow on the bundled demo data, offline. |
| [`libraries`](commands/libraries.md) | List the MSigDB collections and common Enrichr libraries. |
| [`gui`](commands/gui.md) | Start the point-and-click app. |

## Contribute

Questions, bug reports and pull requests are welcome. Open an [issue](https://github.com/helenhuangmath/SimpNetworkViz/issues) on GitHub. For a bug, include the command you ran, the output of `snv --version`, and a small input file that reproduces it if you can share one.

## Author

Hua Huang

## Citation

There is no paper for SimpNetworkViz yet. If you use it, please cite the repository and state the version you used:

> Huang H. SimpNetworkViz: pathway enrichment and EnrichmentMap-style network figures from a gene list. <https://github.com/helenhuangmath/SimpNetworkViz>

Please also cite the methods it builds on, listed in [Theory behind SimpNetworkViz](theory.md#references), and the gene-set collection you used.

## Other useful links

* [Cytoscape](https://cytoscape.org/), with the [EnrichmentMap](https://apps.cytoscape.org/apps/enrichmentmap) and [AutoAnnotate](https://apps.cytoscape.org/apps/autoannotate) apps
* [MSigDB gene-set collections](https://www.gsea-msigdb.org/gsea/msigdb/)
* [Enrichr gene-set libraries](https://maayanlab.cloud/Enrichr/#libraries)
* [STRING](https://string-db.org/)

```{toctree}
:maxdepth: 2
:hidden:

install
tutorial/index
commands/index
guide/index
gallery
workflow
fileformats
theory
qa
api
```
