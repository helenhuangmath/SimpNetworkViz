# Subcommands

The command-line tool is `snv` (`simpnetworkviz` is a longer alias). It has five subcommands:

| Subcommand | Description |
|---|---|
| [`run`](run.md) | Main workflow: gene list, enrichment table or pathway list in; pathway and gene networks out. |
| [`enrich`](enrich.md) | Run only the enrichment test and write the table. |
| [`demo`](demo.md) | Run the full workflow on the bundled demo data, offline. |
| [`libraries`](libraries.md) | List popular gene-set libraries. |
| [`gui`](gui.md) | Start the point-and-click web app. |

General help and the version:

```bash
snv -h              # list the subcommands, with examples
snv run -h          # all options of one subcommand
snv --version
```

```{toctree}
:maxdepth: 1
:hidden:

run
enrich
demo
libraries
gui
```
