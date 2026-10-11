# SimpNetworkViz — summary

**Repository:** `helenhuangmath/SimpNetworkViz`

## Goal

SimpNetworkViz is an easy-to-use tool for these jobs:

1. Take a **gene list** (or a ranked table) and run pathway enrichment analysis (ORA or preranked GSEA).
2. Take **several conditions** (contrasts, time points, single-cell clusters) and compare them in one network.
3. Take **pathways**, either from that enrichment or from a list you give it, and draw two kinds of plot:
   * a clustered **pathway network** in the Cytoscape EnrichmentMap/AutoAnnotate style (like Fig. 6a of Nat Cell Biol 2024, s41556-024-01570-0), with pie-chart nodes for comparisons
   * a matching **gene network**

## What was built

| Part | File | What it does |
|---|---|---|
| Gene sets | `simpnetworkviz/genesets.py` | Loads a GMT file, any Enrichr library (downloaded and cached) or the bundled `demo` set. Reads gene lists, with or without log2FC. |
| Enrichment | `simpnetworkviz/enrichment.py` | Hypergeometric test with Benjamini–Hochberg FDR. Up- and down-regulated genes can be tested separately. Also reads Enrichr, clusterProfiler, g:Profiler and GSEA result tables, and pathway-name lists (with fuzzy matching). |
| GSEA | `simpnetworkviz/gsea.py` | Preranked GSEA: weighted KS statistic, gene-permutation null shared per set size, NES, empirical p, BH FDR (or the classic GSEA q-value), leading-edge genes. ES matches gseapy to 1e-16. |
| Comparison | `simpnetworkviz/compare.py` | `compare_conditions()` runs enrichment per condition and stacks the tables; `merge_conditions()` collapses to one row per pathway with `fdr__<cond>` columns; `enrichment_matrix()` gives a pathway × condition matrix. |
| Pathway map | `simpnetworkviz/pathmap.py` | `build_pathway_map()`: TF-IDF gene-set vectors → SVD → UMAP/t-SNE, Louvain regions with library-weighted word labels; `map_hits()`, `region_summary()`, `plot_pathway_map()` overlay any enrichment or comparison on the fixed map. `snv map`, `--map`. |
| Single cell | `simpnetworkviz/singlecell.py` | Reads Seurat `FindAllMarkers`, scanpy `rank_genes_groups` (long and wide) and AnnData; filters markers and returns `{cluster: {gene: log2FC}}`. |
| Networks | `simpnetworkviz/network.py` | Pathway similarity: ½ Jaccard + ½ overlap coefficient, edges at a cutoff of 0.375 (full sets or hit genes). Grouping: MCL, Louvain, greedy modularity, label propagation, hierarchical, components, by source category, by keyword themes. Labels each cluster automatically. `build_library_network()` draws the whole library with enriched pathways highlighted. Builds gene networks by shared pathways, STRING PPI or gene–pathway (cnet). |
| Layout | `simpnetworkviz/layout.py` | Lays out each cluster separately, then packs clusters in rows with room for labels (like AutoAnnotate). Keeps nodes and labels from overlapping. |
| Colors | `simpnetworkviz/style.py` | Red/blue by direction with shade by FDR, by cluster, significance, NES or log2FC. Per-condition pie slices (condition colors or red/blue). Grey for non-enriched library pathways. Same cluster colors in every figure. |
| Plotting | `simpnetworkviz/plotting.py` | matplotlib figures (PNG/PDF/SVG, editable text): convex-hull cluster outlines, pie-chart nodes, legends (with slice key), cluster bar plot. The figure widens when the legend needs room. |
| Export | `simpnetworkviz/export.py` | Interactive HTML viewer (Cytoscape.js, pie nodes, per-condition details), Cytoscape desktop `.cyjs` plus a style XML, GraphML, CSV tables. |
| Workflow | `simpnetworkviz/pipeline.py` | `run_pipeline()`: gene list, ranked table, enrichment table, pathway list or conditions in; all outputs plus `run_params.json` written. |
| CLI | `simpnetworkviz/cli.py` | `snv run / compare / markers / enrich / demo / libraries / gui` |
| GUI | `simpnetworkviz/app.py` | Streamlit point-and-click app (`snv gui`) with the same input modes |
| Tests / CI | `tests/`, `.github/workflows/tests.yml` | pytest, run on Python 3.9 and 3.12 |

## Validation (`tests/test_validation.py`)

* ORA p-values equal a one-sided Fisher exact test (scipy) and match `gseapy.enrich` with an explicit background.
* BH FDR matches `statsmodels` (or R `p.adjust` reference values).
* GSEA enrichment scores and leading-edge genes match `gseapy.prerank` exactly; NES within 15 % (different permutations).
* Similarity follows the EnrichmentMap definition.
* The README numbers for `snv demo` (17 of 20 pathways, 17 nodes, 16 edges, 5 clusters) are pinned.
* Live: `MSigDB_Hallmark_2020` downloads from Enrichr with 50 sets and the expected E2F genes; the demo genes recover E2F/G2M/OXPHOS hallmarks; STRING returns CDK1–CCNB1 and TP53–MDM2. Skipped automatically when offline.

## How to run

```bash
pip install -e .
snv demo --compare                              # offline example -> demo_output/ and demo_output/compare/
snv run genes.txt --library MSigDB_Hallmark_2020
snv run ranked.csv --method gsea --library MSigDB_Hallmark_2020
snv compare KO=ko.csv WT=wt.csv --library Reactome_2022
snv markers markers.csv --library MSigDB_Hallmark_2020 --whole-library
snv gui                                         # needs: pip install -e ".[gui]"
```

## Notes

* The bundled `demo` gene sets are **shortened versions of real pathways, for testing only**. Use full libraries (Enrichr or MSigDB GMT) for real analyses.
* Pass a background gene list (`--background`) for correct enrichment statistics.
* GSEA needs all measured genes with a statistic, not a filtered DEG list.
