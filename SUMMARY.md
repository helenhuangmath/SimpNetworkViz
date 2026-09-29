# SimpNetworkViz — project summary

**Repository:** `helenhuangmath/SimpNetworkViz`
**Status:** 25 automated tests pass. Offline demo runs end to end.

## Goal

SimpNetworkViz is an easy-to-use tool for two jobs:

1. Take a **gene list** and run pathway enrichment analysis.
2. Take **pathways**, either from that enrichment or from a list you give it, and draw two kinds of plot:
   * a clustered **pathway network** in the Cytoscape EnrichmentMap/AutoAnnotate style (like Fig. 6a of Nat Cell Biol 2024, s41556-024-01570-0)
   * a matching **gene network**

## What was built

| Part | File | What it does |
|---|---|---|
| Gene sets | `simpnetworkviz/genesets.py` | Loads a GMT file, any Enrichr library (downloaded and cached) or the bundled `demo` set. Reads gene lists, with or without log2FC. |
| Enrichment | `simpnetworkviz/enrichment.py` | Hypergeometric test with Benjamini–Hochberg FDR. Up- and down-regulated genes can be tested separately. Also reads Enrichr, clusterProfiler, g:Profiler and GSEA result tables, and pathway-name lists (with fuzzy matching). |
| Networks | `simpnetworkviz/network.py` | Pathway similarity: ½ Jaccard + ½ overlap coefficient, edges at a cutoff of 0.375. Clusters with MCL (or Louvain). Labels each cluster automatically. Builds gene networks by shared pathways, STRING PPI or gene–pathway (cnet). |
| Layout | `simpnetworkviz/layout.py` | Lays out each cluster separately, then packs clusters in rows with room for labels (like AutoAnnotate). Keeps nodes and labels from overlapping. |
| Colors | `simpnetworkviz/style.py` | Red/blue by direction with shade by FDR, or by cluster, significance, NES or log2FC. Same cluster colors in every figure. |
| Plotting | `simpnetworkviz/plotting.py` | matplotlib figures (PNG/PDF/SVG, editable text): convex-hull cluster outlines, legends, cluster bar plot. |
| Export | `simpnetworkviz/export.py` | Interactive HTML viewer (Cytoscape.js), Cytoscape desktop `.cyjs` plus a style XML, GraphML, CSV tables. |
| Workflow | `simpnetworkviz/pipeline.py` | `run_pipeline()`: gene list, enrichment table or pathway list in; all outputs written. |
| CLI | `simpnetworkviz/cli.py` | `snv run / enrich / demo / libraries / gui` |
| GUI | `simpnetworkviz/app.py` | Streamlit point-and-click app (`snv gui`) |
| Tests / CI | `tests/`, `.github/workflows/tests.yml` | pytest, run on Python 3.9 and 3.12 |

## Change history (requests in this session)

1. **Initial build.** Enrichment, pathway network, gene networks, CLI, GUI, exports, tests, README.
2. **Arial, no bold.** All figure text is Arial at regular weight. Liberation Sans (identical metrics) is the fallback where Arial is missing. The HTML viewer uses Arial too.
3. **Gene network tied to the pathway network.** Added a two-panel (a | b) figure with matching cluster colors.
4. **Separate figures, bigger text.**
   * The pathway and gene networks are saved as separate files again. The two-panel function is still in the API.
   * Cluster outlines use cluster colors by default, so the separate figures match.
   * Text sizes: cluster labels 14 pt, node labels 9 pt, legend 11–13 pt, title 18 pt.
   * Outlines now enclose node labels.

## How to run

```bash
pip install -e .
snv demo                                        # offline example -> demo_output/
snv run genes.txt --library MSigDB_Hallmark_2020
snv run --pathways my_pathways.txt --library Reactome_2022
snv gui                                         # needs: pip install -e ".[gui]"
```

Outputs:

* `pathway_network.*` and `gene_network.*`: png, pdf, svg, html, cyjs, graphml, csv
* `gene_pathway_network.*`
* `cluster_barplot.png`
* `enrichment.csv`

## Not yet verified (the build sandbox had no internet)

* Downloading Enrichr libraries. The code is written but was not tested live.
* STRING PPI queries. Tested only with a mocked response.
* Loading Cytoscape.js from a CDN in the HTML viewer. The viewer was verified with a local copy of the library.
* Importing the generated `_style.xml` into Cytoscape desktop.
* Embedding real Arial in PDFs. The sandbox had no Arial, so the PDFs there contain Liberation Sans.

## Notes

* The bundled `demo` gene sets are **shortened versions of real pathways, for testing only**. Use full libraries (Enrichr or MSigDB GMT) for real analyses.
* Pass a background gene list (`--background`) for correct enrichment statistics.
