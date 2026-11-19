# Publication figures and Cytoscape

This tutorial covers what to change once the network itself is right: formats, text, colors, and finishing by hand.

## Choose the format

| Format | Use it for |
|---|---|
| PDF, SVG | Journal figures. Vector files look sharp at any size, and the text stays editable in Illustrator or Inkscape. |
| PNG | Slides and quick viewing. Saved at 300 dpi. |
| HTML | Exploring the network yourself, or sharing it with collaborators. |

Choose formats on the command line with `--formats`, for example `--formats pdf,svg,csv`. In Python, the file extension you pass to `plot()` sets the format.

For a sharper PNG:

```bash
snv run degs.csv --library MSigDB_Hallmark_2020 --dpi 600
```

```python
pnet.plot("pathway_network.png", dpi=600)
```

## Choose the colors

| Option | Look |
|---|---|
| `--color-by direction` | Red for up, blue for down, shaded by −log10 FDR. The classic EnrichmentMap look. |
| `--color-by cluster` | One pastel color per theme. |
| `--color-by significance` | One color scale, by −log10 FDR. |
| `--color-by nes` | By normalized enrichment score, for GSEA tables. |

Cluster outlines use the cluster colors, so that the pathway figure and the gene figure match. For grey outlines:

```python
pnet.plot("pathway_network.pdf", hull_style="neutral")
```

## Change the text

All figure text is Arial at regular weight. Where Arial is not installed, Liberation Sans is used; it has the same metrics, so the layout does not change.

The default sizes are chosen for a figure up to about 16 inches wide: cluster labels 26 pt, node labels 13 pt, gene labels 18 pt, legend 15–17 pt and title 24 pt. To change them for one figure:

```python
pnet.plot("pathway_network.pdf", cluster_label_size=20, node_label_size=10)
```

To choose which nodes are labeled:

```python
pnet.plot("pathway_network.pdf", node_labels="top")     # the best pathway of each cluster
pnet.plot("pathway_network.pdf", node_labels="none")    # cluster labels only
```

`node_labels` also takes `"all"`, `"auto"` or a list of node names.

To set the largest figure size, in inches:

```python
pnet.plot("pathway_network.pdf", max_width=10, max_height=12)
```

## Rename the clusters

The automatic labels are a starting point. Replace any of them with your own wording, in Python:

```python
print(pnet.cluster_table()[["cluster", "label"]])
pnet.group_labels[1] = "ECM remodelling"
pnet.group_labels[4] = "Interferon signalling"
pnet.plot("pathway_network.pdf")
```

Build the gene network after renaming, so that it carries the new labels.

## Draw fewer things

A figure for print usually needs fewer nodes than a figure for exploring.

* Fewer pathways: `--max-pathways 40`, or a stricter `--fdr`.
* Fewer genes: `--max-genes 80`.
* Only some clusters in the gene network, in Python: `build_gene_network(pnet, gene_sets, clusters=[1, 3])`.
* No legend: `pnet.plot(..., legend=False)`.

## Edit the PDF or SVG

Open the PDF or SVG in Illustrator, Inkscape or Affinity Designer. Text is stored as text, not as outlines, so you can retype labels and change fonts. Nodes, edges and outlines are separate vector objects, so you can move a label that sits too close to an edge.

## Fine-tune in Cytoscape

Every run writes a `.cyjs` network file and a `_style.xml` style file. To continue in [Cytoscape](https://cytoscape.org/) desktop:

1. Import the network: *File ▸ Import ▸ Network from File*, and choose `pathway_network.cyjs`.
2. Import the style: *File ▸ Import ▸ Styles from File*, and choose `pathway_network_style.xml`.
3. In the *Style* panel, select the **SimpNetworkViz** style.
4. Optionally run AutoAnnotate. Choose the existing `cluster` column as the cluster source, so that the clusters match the figure.

The node positions, colors, sizes and cluster numbers come across, and you can then use any Cytoscape layout or app.

## Share an interactive figure

The HTML viewer is a single file. By default it loads the Cytoscape.js library from the internet when opened. To make a file that works anywhere, including as a journal supplement, download `cytoscape.min.js` once and embed it:

```bash
snv run degs.csv --library MSigDB_Hallmark_2020 --cytoscape-js cytoscape.min.js
```

## Report the method

For the methods section of a paper, report:

* the SimpNetworkViz version, from `snv --version`
* the gene-set library and its release
* the background gene list
* the options you changed from the defaults: FDR cutoff, similarity measure and cutoff, clustering method

A sentence you can adapt:

> Over-representation analysis was run with SimpNetworkViz (version X) using a hypergeometric test with Benjamini–Hochberg correction, against [library, release], with all expressed genes as the background. Pathways with FDR ≤ 0.05 were connected when their gene-set similarity (the mean of the Jaccard and overlap coefficients) was at least 0.375, and clustered with MCL.

The steps are described in [How the network is built](../method.md).
