# demo

## Overview

`snv demo` runs the full workflow on data bundled with the package: 88 genes with log2 fold changes and 32 shortened gene sets. It needs no input files and no internet connection, so it is the quickest way to check an install and to see every kind of output.

:::{warning}
The bundled `demo` gene sets are shortened versions of real pathways, meant only for testing and demos. Use full libraries for real analyses.
:::

## Examples of command-line usage

```bash
snv demo                     # writes to demo_output/
snv demo --compare           # also a three-condition comparison in demo_output/compare/
```

With `--compare`, the up-regulated genes, the down-regulated genes and all genes are compared as three conditions on the whole demo library.

The [guided tutorial](../tutorial/demo.md) goes through the same analysis step by step.

## Command-line options

```{eval-rst}
.. argparse::
   :module: _cli_for_docs
   :func: build_parser
   :prog: snv
   :path: demo
   :noepilog:
```
