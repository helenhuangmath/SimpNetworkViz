# Global pathway maps of the MSigDB collections

Two sets of maps: **mouse** (`mouse/`, MSigDB 2025.1 native mouse collections, the right choice for the
mouse `Cl13_vs_Arm` example) and **human** (the other folders, MSigDB v7.4). Both were built with
`snv map --by-category --method both`; the mouse ones straight from the package's MSigDB fetcher:

```bash
for c in Hallmark GO CP CellType TFT; do snv map --msigdb $c --species mouse --by-category --method both -o examples/msigdb_maps/mouse/$c; done
```

## Mouse (MSigDB 2025.1)

| Folder | Collection | Sets on the map | Regions |
|---|---|---|---|
| `mouse/Hallmark/` | Hallmark (mh.all) | 50 | 4 |
| `mouse/GO/` | GO biological process (m5.go.bp) | 5489 | 35 |
| | GO cellular component (m5.go.cc) | 679 | 15 |
| | GO molecular function (m5.go.mf) | 1181 | 22 |
| `mouse/CP/` | Reactome (m2.cp.reactome) | 1017 | 22 |
| | WikiPathways (m2.cp.wikipathways) | 169 | 8 |
| | BioCarta (m2.cp.biocarta) | 184 | 8 |
| `mouse/CellType/` | Cell-type signatures (m8.all) | 195 | 6 |
| `mouse/TFT/` | GTRD transcription-factor targets (m3.gtrd) | 140 | 5 |
| `mouse/overall/` | GO BP + GO CC + GO MF + Reactome + WikiPathways in one embedding | 8535 | 36 |

MSigDB has no mouse KEGG, PID, HPO, oncogenic (C6) or immunologic (C7) collections; for those use the
human maps below (symbols are matched case-insensitively, with some loss) or an Enrichr mouse library
(`--library KEGG_2019_Mouse`). The mouse GO BP map is the same as `../Cl13_vs_Arm/genesets/gobp_map.csv`.

## Human (MSigDB v7.4)

One t-SNE and one UMAP map per collection, built from the **human MSigDB v7.4 symbol GMT files** (`h.all`, `c2.cp`, `c5.all`, `c3.tft`, `c6.all`, `c7.all`,
downloaded from <https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp>). Every map embeds all
sets of a collection with 10–500 genes from their gene-set similarity (TF-IDF gene vectors → 50-d SVD →
t-SNE or UMAP), divided into Louvain regions in the 50-d space, so the t-SNE and UMAP maps of one
collection share their regions and labels and differ only in the 2-D layout.

| Folder | Collection(s) | Sets on the map | Regions |
|---|---|---|---|
| `h.all/` | Hallmark | 50 | 5 |
| `c5.all/` | GO biological process | 5161 | 37 |
| | GO cellular component | 643 | 15 |
| | GO molecular function | 1064 | 21 |
| | Human phenotype (HPO) | 3468 | 31 |
| `c3.tft/` | TF targets, GTRD (`*_TARGET_GENES`) | 272 | 10 |
| | TF targets, legacy motifs (`c3.tft`) | 571 | 19 |
| `c6.all/` | Oncogenic signatures (`c6.all`) | 189 | 10 |
| `c7.all/` | ImmuneSigDB (`GSE*`) | 4854 | 23 |
| | Vaccine response (`c7.all`) | 274 | 11 |
| `c2.cp/` | Reactome | 1260 | 28 |
| | WikiPathways | 537 | 13 |
| | BioCarta | 226 | 9 |
| | PID | 196 | 7 |
| | KEGG (legacy) | 186 | 8 |
| `overall/` | GO BP + GO CC + GO MF + KEGG + Reactome + WikiPathways in **one** embedding | 8851 | 36 |

The `overall/` maps embed the six pathway databases together (`overall_map_tsne.csv`, `overall_map_umap.csv`, built
with plain `snv map` on the concatenated GMT; the PDFs color the points by database with
`plot_pathway_map(..., color_by_category=True)`). Use them when an analysis mixes databases, or to see how the
databases relate: GO BP forms the bulk, Reactome and KEGG sit inside the matching GO regions, GO MF and GO CC
occupy their own areas because they describe functions and locations rather than processes.

Each folder contains, per collection, `<collection>_map_tsne.csv`, `<collection>_map_umap.csv` and a PDF
of each, plus `category_maps.csv` (index) and `category_maps.pdf` (all panels side by side).

## Use

Overlay any analysis run with the same collection on a map:

```bash
snv run genes.txt --gmt c5.all.v7.4.symbols.gmt --map examples/msigdb_maps/c5.all/go_biological_process_map_umap.csv
snv compare KO=ko.csv WT=wt.csv --gmt h.all.v7.4.symbols.gmt --map examples/msigdb_maps/h.all/hallmark_map_tsne.csv
```

Keep the map with the GMT release it was built from: a set that is not in the map (other release, other
species, outside 10–500 genes) is dropped from the overlay. To rebuild for another release or for mouse
(e.g. `msigdb.v2025.1.Mm.symbols.gmt`, which holds all collections in one file):

```bash
snv map --gmt msigdb.v2025.1.Mm.symbols.gmt --by-category --method both -o examples/msigdb_maps/msigdb_Mm
```

Note that HPO region labels are generic ("Abnormality abnormal ...") because nearly every HPO term name
starts with those words; edit the `cluster_label` column of the CSV if you want better names.
