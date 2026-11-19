# gui

## Overview

`snv gui` starts a point-and-click web app in your browser. It runs the same workflow as [`snv run`](run.md), for people who prefer not to use the command line. The app runs on your own computer; your data are not uploaded anywhere by the app itself.

The app needs streamlit, which is an optional dependency:

```bash
pip install -e ".[gui]"
snv gui
```

A browser page opens. Stop the app with Ctrl+C in the terminal.

## Using the app

The sidebar has three steps and a button.

1. **Input.** Choose where to start from:

   * **Gene list:** paste genes, one per line. An optional second column, separated by a tab or a comma, is the log2FC. Tick **Use demo genes** to try the bundled data.
   * **Enrichment table:** upload a CSV, TSV or XLSX result from Enrichr, clusterProfiler, g:Profiler or GSEA. You can also paste a gene list with log2FC to color the gene nodes.
   * **Pathway list:** paste pathway names, one per line. You can also paste genes of interest to restrict the gene network.

2. **Gene sets.** Pick `demo`, one of the popular libraries, **Upload GMT** or **Other Enrichr library**. Paste background genes here; this is optional but recommended.

3. **Network.** Set the FDR cutoff, the maximum number of pathways, the similarity measure and cutoff, the clustering method, the pathway colors, the gene network type and the maximum number of genes. These are the same settings as the [`run` options](run.md#options-controlling-the-pathway-network).

4. Click **Build networks**.

The results appear in five tabs:

| Tab | Contents |
|---|---|
| Pathway network | The figure, and an interactive view you can expand |
| Gene network | The gene network figure and its interactive view |
| Gene–pathway network | The gene–pathway (cnet) figure |
| Tables | The enrichment table and the pathway node table |
| Download | Every figure and table as one zip, or file by file |

## Notes

* Choosing a library by name downloads it from Enrichr, and the `string` gene network sends the gene names to STRING. Everything else runs locally.
* To use the app on a remote server, forward the port that streamlit prints (8501 by default) over SSH and open the forwarded address in your local browser.
