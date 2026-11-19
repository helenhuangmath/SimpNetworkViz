# map

## Overview

`snv map` embeds every gene set of a library in two dimensions and saves the result as a CSV. Pass that file with `--map` to `snv run`, `snv compare` or `snv markers` to draw any analysis on the same coordinates.

## Examples of command-line usage

1. Build the map of one library (a PDF with the same name is written too):

   ```bash
   snv map --gmt m5.go.bp.v2025.1.Mm.symbols.gmt -o gobp_map.csv
   ```

2. Build it from an MSigDB collection fetched by the package:

   ```bash
   snv map --msigdb GO_BP --species mouse -o gobp_map.csv
   ```

3. One t-SNE and one UMAP map per MSigDB collection:

   ```bash
   snv map --gmt msigdb.v2025.1.Mm.symbols.gmt --by-category -o msigdb_maps
   ```

4. Use a map in an analysis:

   ```bash
   snv run genes.txt --gmt m5.go.bp.v2025.1.Mm.symbols.gmt --map gobp_map.csv
   ```

UMAP and t-SNE need the `map` extra: `pip install "simpnetworkviz[map]"`.

How the map is built, how to read it and how to use the interactive version are described in [The whole-library view and the pathway map](../guide/library_map.md).

## Command-line options

```{eval-rst}
.. argparse::
   :module: _cli_for_docs
   :func: build_parser
   :prog: snv
   :path: map
   :noepilog:
```

## Output files

* `NAME.csv`: the map, with one row per gene set, its coordinates and its region. Keep it with the library it was built from.
* `NAME.pdf`: the map drawn without hits.
* `NAME.html`: the interactive map.
* With `--by-category`: a folder with `<collection>_map_tsne.csv` and `_umap.csv`, a PDF of each, the index `category_maps.csv` and the overview `category_maps.pdf`.
