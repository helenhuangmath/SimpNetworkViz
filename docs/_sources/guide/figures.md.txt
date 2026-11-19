# Customizing the figures

* **Colors.** `--color-by direction` (red/blue, shade = FDR), `cluster` (one pastel color per cluster), `significance` (single scale), `nes` (GSEA), `condition` (comparisons). Gene networks are colored by log2FC when available.
* **Formats, size and text.** Figures are written as PDF and SVG (vector, editable text); add `png` to `--formats` for a raster copy, with `--dpi` setting its resolution (300 by default). Cluster labels 14 pt, node labels 9 pt, legend 11–13 pt, all Arial (Liberation Sans where Arial is missing); the layout leaves room for them. `--title "..."` adds a title. In Python: `plot(cluster_label_size=16, node_label_size=8, max_width=12)`.
* **Node labels.** Shown on pathway networks with ≤ 30 nodes; `plot(node_labels="top")` labels only the best pathway per cluster, `"all"` or `"none"` as named.
* **Cluster outlines.** Each cluster has its own color so pathway and gene figures match; `plot(hull_style="neutral")` gives grey outlines.
* **Rename a cluster.** In Python, `pnet.group_labels[3] = "Interferon signalling"` then `pnet.plot(...)`; or edit the text in the PDF/SVG.
* **Gene network.** `--gene-network shared` (default, offline), `string` (STRING interactions; `--species 10090` for mouse), `cnet` (gene–pathway only) or `none`. `--max-genes` limits the genes shown (default 150).
* **Reproducible layout.** The layout is deterministic for a given `--seed` (default 0).
* **Offline HTML.** `--cytoscape-js path/to/cytoscape.min.js` embeds the JavaScript library in the viewer.

## Fine-tuning in Cytoscape

1. *File ▸ Import ▸ Network from File*: choose `pathway_network.cyjs`. Positions are kept.
2. *File ▸ Import ▸ Styles from File*: choose `pathway_network_style.xml`, then apply the **SimpNetworkViz** style.
3. Optionally run AutoAnnotate on the existing `cluster` column.

For comparisons, the `fdr__<condition>` columns are exported so charts can be added in Cytoscape.
