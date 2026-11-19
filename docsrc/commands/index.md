# Subcommands

The command-line tool is `snv` (`simpnetworkviz` is a longer alias). `snv -h` lists the subcommands, and `snv run -h` prints every option of one.

| Subcommand | Description |
|---|---|
| [`run`](run.md) | Gene list, ranked table, enrichment table or pathway list in; pathway and gene networks out. |
| [`compare`](compare.md) | Several conditions in one network, with pie-chart nodes. |
| [`markers`](markers.md) | Single-cell marker table in; one pie slice per cell cluster. |
| [`map`](map.md) | Build a fixed 2-D map of a whole gene-set library, to reuse with `--map`. |
| [`enrich`](enrich.md) | Run only the enrichment test and write the table. |
| [`demo`](demo.md) | Run the full workflow on the bundled demo data, offline. |
| [`libraries`](libraries.md) | List the MSigDB collections and common Enrichr libraries. |
| [`gui`](gui.md) | Start the point-and-click app. |

```{toctree}
:maxdepth: 1
:hidden:

run
compare
markers
map
enrich
demo
libraries
gui
```

## Options shared by `run`, `compare` and `markers`

| Option | Default | Meaning |
|---|---|---|
| `--library NAME[,NAME]`, `--gmt FILE` | | gene sets |
| `--background FILE` | all library genes | universe for the test |
| `--method ora\|gsea` | ora | hypergeometric test or preranked GSEA |
| `--no-split` | off | do not test up/down genes separately |
| `--min-size`, `--max-size` | 5 / 1000 (ORA), 15 / 500 (GSEA) | gene-set size limits |
| `--permutations` | 1000 | GSEA gene permutations |
| `--fdr` | 0.05 | FDR cutoff for pathways in the network |
| `--max-pathways` | 100 | most significant pathways kept |
| `--similarity combined\|jaccard\|overlap` | combined | gene-set similarity for edges |
| `--similarity-genes full\|hits` | full | full gene sets or only hit genes |
| `--cutoff` | 0.375 | minimum similarity for an edge |
| `--cluster METHOD` | mcl | grouping method ([Grouping pathways](../guide/grouping.md)) |
| `--inflation`, `--resolution`, `--linkage-threshold` | 2 / 1 / 1 − cutoff | tuning of mcl / louvain, greedy / hierarchical |
| `--categories FILE`, `--themes FILE` | | for `--cluster category` / `keyword` |
| `--whole-library` | off | draw all pathways, highlight enriched |
| `--map MAP.csv` | | also draw the hits on a library map from `snv map` |
| `--color-by` | auto | cluster, direction, significance, nes, condition |
| `--gene-network shared\|string\|cnet\|none` | shared | gene network type |
| `--max-genes` | 150 | genes in the gene network |
| `--msigdb NAME` | | MSigDB collection for `--species`, downloaded once (repeatable; `snv libraries` lists them) |
| `--species` | human | human or mouse: MSigDB collections for `--msigdb` and the STRING taxon |
| `--formats` | pdf,svg,html,cyjs,graphml,csv | outputs to write (add `png` for raster figures) |
| `--title`, `--dpi`, `--seed`, `-o` | | title, PNG resolution if requested (300), layout seed (0), output folder |
| `snv markers` only: `--padj`, `--min-logfc`, `--top-genes`, `--direction`, `--clusters` | 0.05 / 0.25 / 200 / up / all | marker filters |
