# Grouping pathways

The default grouping is Markov clustering (MCL) on the gene-overlap network, as in AutoAnnotate, with labels built from the words the member pathways share. Alternatives, chosen with `--cluster`:

| Method | What it does | Tuning |
|---|---|---|
| `mcl` (default) | Markov clustering on the similarity graph | `--inflation` (higher = more, smaller clusters; default 2) |
| `louvain`, `greedy` | Community detection | `--resolution` (higher = more clusters; default 1) |
| `label_propagation` | Fast community detection | |
| `hierarchical` | Average-linkage clustering on gene-set distance; also groups pathways related through intermediates | `--linkage-threshold` (default 1 − cutoff) |
| `components` | Connected components of the graph | |
| `category` | By source database: Hallmark, KEGG, Reactome, GO ..., read from the name, or from `--categories FILE` (`pathway<TAB>category`, e.g. MSigDB sub-collections) | |
| `keyword` | By your own themes from `--themes FILE` | |
| `none` | No grouping | |

A themes file has one theme per line with its keywords:

```
Proliferation	cycle, mitotic, E2F, replication, checkpoint
Immune	interferon, cytokine, antigen
Matrix	matrix, collagen, adhesion
```

The first matching theme wins, so list specific themes before general ones. Pathways matching nothing stay unclustered.

Two related switches: `--similarity jaccard|overlap|combined` changes the similarity measure, and `--similarity-genes hits` computes it from your hit genes only instead of the full gene sets (useful when you care about which of *your* genes the pathways share).
