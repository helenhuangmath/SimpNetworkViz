# demo

## Overview

`snv demo` runs the full workflow on data bundled with the package. It needs no input files and no internet connection, so it is the quickest way to check an install and to see every kind of output.

The demo data are 88 genes with log2FC (53 up, 35 down) and 32 gene sets.

:::{warning}
The demo gene sets are shortened versions of real pathways, meant only for testing and demos. Use full libraries for real analyses.
:::

## Examples of command-line usage

```bash
snv demo                    # writes to demo_output/
snv demo -o my_demo --dpi 150
```

`snv demo` is equivalent to running [`snv run`](run.md) on the bundled gene list with `--library demo --cutoff 0.3`.

## Command-line options

* `-o, --outdir DIR`

  The output folder. Default: `demo_output`.

* `--dpi INT`

  PNG resolution. Default: 300.

* `--compare`

  Also compare the up-regulated, down-regulated, and combined demo gene lists in the `compare` subdirectory.

* `--cytoscape-js FILE`

  A local copy of `cytoscape.min.js` to embed in the HTML viewer.

## Output

The same files as [`snv run`](run.md#output-files). The [demo walkthrough](../tutorial/demo.md) explains each one.
