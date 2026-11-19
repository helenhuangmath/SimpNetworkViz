# How the network is built

The defaults follow the Cytoscape apps EnrichmentMap and AutoAnnotate. This page describes each step and the option that controls it.

## Enrichment test

When you start from a gene list, each gene set is tested by over-representation analysis.

* **Test.** A one-sided hypergeometric test: given the background, is the number of your genes in the gene set larger than expected by chance?
* **Background.** The genes given with `--background`. Without it, every gene in the gene sets. Your genes and the gene sets are both restricted to the background first.
* **Size limits.** Gene sets with fewer than 5 or more than 1000 genes in the background are not tested. Pathways with fewer than two of your genes are not reported.
* **Multiple testing.** Benjamini–Hochberg FDR, over the tested gene sets.
* **Direction.** When fold changes are given, genes with positive and negative values are tested as two separate lists, and each result is marked `up` or `down`. The FDR is computed within each list.

## Pathway network

1. **Filter.** Keep pathways with FDR ≤ 0.05 (`--fdr`), up to the 100 most significant (`--max-pathways`).

2. **Similarity.** Compare every pair of pathways by their full gene sets, with the combined similarity

   ```text
   similarity = ½ · Jaccard + ½ · overlap coefficient

   Jaccard             = |A ∩ B| / |A ∪ B|
   overlap coefficient = |A ∩ B| / min(|A|, |B|)
   ```

   An edge is drawn when the similarity is at least 0.375 (`--cutoff`). `--similarity jaccard` or `overlap` uses one coefficient alone. These are the EnrichmentMap defaults.

3. **Clusters.** Cluster the network with Markov clustering, MCL, with inflation 2.0. `--cluster louvain` uses Louvain community detection, and `--cluster components` makes each connected group a cluster. A pathway left on its own is not a cluster and is drawn without an outline. Clusters are numbered from the most significant.

4. **Labels.** A cluster of three or more pathways is labeled with up to three informative words that occur in at least two of its pathway names, in the order they usually appear. Common words such as "regulation" or "pathway" are skipped. Smaller clusters, and clusters with no shared word, are labeled with their most significant pathway. This is similar to the WordCloud labels of AutoAnnotate.

5. **Layout.** Each cluster is laid out on its own. The clusters are then packed so that they do not overlap and each label has room, like AutoAnnotate's "layout to prevent cluster overlap".

6. **Styling.**

   | Element | Shows |
   |---|---|
   | Node size | Gene-set size |
   | Node color | Direction and −log10 FDR, or cluster |
   | Edge width | Similarity |
   | Outline | Cluster |

## Gene networks

Gene networks are built from the pathway network.

* **Genes.** Your genes that fall in the drawn pathways. If you give only a pathway list and no genes, the full gene sets are used. At most 150 genes are drawn (`--max-genes`).
* **Cluster of a gene.** Each gene is assigned to the pathway cluster it belongs to most. The gene module then carries the same color and label as that cluster in the pathway figure.
* **Edges.** They depend on the type of network:

  | `--gene-network` | Edge between two nodes |
  |---|---|
  | `shared` | The two genes belong to at least one common pathway. The weight is the number of shared pathways. |
  | `string` | STRING reports an interaction between the two proteins, with a combined score of at least 400. |
  | `cnet` | A gene and a pathway it belongs to. |

* **Node color.** The log2FC from your gene list, if given.

## Differences from Cytoscape

SimpNetworkViz reimplements the method; it does not run Cytoscape. The clusters and labels will usually be close to what EnrichmentMap and AutoAnnotate give for the same cutoffs, but they are not guaranteed to be identical, and the layout is different. To continue in Cytoscape itself, see [Fine-tune in Cytoscape](tutorial/figures.md#fine-tune-in-cytoscape).

## References

* Merico D. et al. Enrichment Map: a network-based method for gene-set enrichment visualization and interpretation. *PLoS ONE* 5, e13984 (2010).
* Kucera M. et al. AutoAnnotate: a Cytoscape app for summarizing networks with semantic annotations. *F1000Research* 5, 1717 (2016).
* Enright A.J., Van Dongen S., Ouzounis C.A. An efficient algorithm for large-scale detection of protein families. *Nucleic Acids Res.* 30, 1575–1584 (2002).
