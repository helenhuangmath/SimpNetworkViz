# Get started

These tutorials go through each kind of analysis step by step, with the code, its output and the figures. All but one run offline on data that comes with the package or the repository.

New to SimpNetworkViz? Start with the guided tutorial.

## Introductory tutorials

::::{grid} 1 2 3 3
:gutter: 3

:::{grid-item-card} Guided tutorial: from a gene list to a pathway network
:link: demo
:link-type: doc
:img-top: ../images/vignettes/guided_pathway_direction.png

Enrichment, the pathway network, the gene networks and the exports, one step at a time.
:::

:::{grid-item-card} Preranked GSEA
:link: gsea
:link-type: doc
:img-top: ../images/vignettes/gsea_nes.png

Rank all measured genes, run GSEA and color the network by enrichment score.
:::

:::{grid-item-card} Starting from other tools' results or a pathway list
:link: other_tools
:link-type: doc
:img-top: ../images/vignettes/pathways_cnet.png

Draw the network for Enrichr, clusterProfiler, g:Profiler or GSEA results, or for pathways you choose.
:::
::::

## Comparing conditions

::::{grid} 1 2 3 3
:gutter: 3

:::{grid-item-card} Comparing conditions in one network
:link: compare
:link-type: doc
:img-top: ../images/vignettes/compare_condition.png

Several gene lists in one network, with a pie-chart node for each pathway.
:::

:::{grid-item-card} Single-cell marker genes
:link: singlecell
:link-type: doc
:img-top: ../images/vignettes/singlecell_library.png

From a Seurat or scanpy marker table to one pie slice per cell cluster.
:::

:::{grid-item-card} Cl13 vs Arm: gene groups of two conditions
:link: cl13_vs_arm
:link-type: doc
:img-top: ../images/tutorial/pooled_GOBP.png

A full comparison on real data: exhausted versus memory CD8 T cells.
:::
::::

## Whole libraries and figures

::::{grid} 1 2 3 3
:gutter: 3

:::{grid-item-card} The pathway map
:link: pathway_map
:link-type: doc
:img-top: ../../docs/images/pathway_map_umap.png

Place a whole gene-set library on a fixed 2-D map and draw every analysis on it.
:::

:::{grid-item-card} Grouping pathways and styling figures
:link: grouping
:link-type: doc
:img-top: ../images/vignettes/grouping_category.png

Other clustering methods, grouping by database or by your own themes, colors and labels.
:::
::::

## Interactive examples

Results you can explore in the browser without installing anything:

* [Cl13 vs Arm on the GO biological process map](https://helenhuangmath.github.io/SimpNetworkViz/examples/cl13-vs-arm.html)
* [Six Cl13 gene groups on the GO biological process map](https://helenhuangmath.github.io/SimpNetworkViz/examples/within-cl13.html)
* [Six Cl13 gene groups on the Hallmark map](https://helenhuangmath.github.io/SimpNetworkViz/examples/hallmark-cl13-subgroups.html)

## Where to go next

* The [standard workflow](../workflow.md) is a short reference for the usual single-list analysis.
* The [subcommand pages](../commands/index.md) list every option.
* The [user guide](../guide/index.md) covers the pathway map, grouping and figures in full.

```{toctree}
:maxdepth: 1
:hidden:

demo
gsea
other_tools
compare
singlecell
cl13_vs_arm
pathway_map
grouping
```
