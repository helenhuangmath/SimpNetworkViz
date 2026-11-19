# libraries

## Overview

`snv libraries` prints the MSigDB collections available for each species and a short list of widely used Enrichr libraries. It takes no options.

```bash
snv libraries
```

```text
MSigDB 2025.1 collections (--msigdb NAME --species human|mouse; downloaded once):

  collection       human                  mouse
  Hallmark         h.all                  mh.all
  GO               c5.go                  m5.go
  GO_BP            c5.go.bp               m5.go.bp
  GO_CC            c5.go.cc               m5.go.cc
  GO_MF            c5.go.mf               m5.go.mf
  Reactome         c2.cp.reactome         m2.cp.reactome
  WikiPathways     c2.cp.wikipathways     m2.cp.wikipathways
  BioCarta         c2.cp.biocarta         m2.cp.biocarta
  KEGG             c2.cp.kegg_medicus     -
  KEGG_legacy      c2.cp.kegg_legacy      -
  PID              c2.cp.pid              -
  CP               c2.cp                  m2.cp
  CGP              c2.cgp                 m2.cgp
  TFT              c3.tft                 m3.gtrd
  miRNA            c3.mir                 m3.mirdb
  HPO              c5.hpo                 -
  TumorPhenotype   -                      m5.mpt
  Oncogenic        c6.all                 -
  ImmuneSigDB      c7.immunesigdb         -
  Vaccine          c7.vax                 -
  CellType         c8.all                 m8.all
  Positional       c1.all                 m1.all
  All              msigdb                 msigdb
  (mouse collections are MSigDB's native mouse sets; KEGG, PID, HPO, C6, C7 exist for human only)

Popular Enrichr libraries (any Enrichr library name works with --library):

  MSigDB_Hallmark_2020           MSigDB Hallmark gene sets (50 broad processes)
  KEGG_2021_Human                KEGG pathways (human)
  Reactome_2022                  Reactome pathways (human)
  WikiPathways_2023_Human        WikiPathways (human)
  GO_Biological_Process_2023     Gene Ontology biological process
  GO_Molecular_Function_2023     Gene Ontology molecular function
  GO_Cellular_Component_2023     Gene Ontology cellular component
  KEGG_2019_Mouse                KEGG pathways (mouse)
  WikiPathways_2019_Mouse        WikiPathways (mouse)
  BioPlanet_2019                 NCATS BioPlanet pathways

Full list: https://maayanlab.cloud/Enrichr/#libraries
Local GMT files: --gmt FILE (repeatable)
Offline test set: --library demo
```

## Using a library

| Source | Option | Notes |
|---|---|---|
| MSigDB collection | `--msigdb NAME --species human\|mouse` | Downloaded once from the Broad release server. Repeat `--msigdb` for several collections. |
| Enrichr library | `--library NAME[,NAME]` | Any [Enrichr library name](https://maayanlab.cloud/Enrichr/#libraries). Mostly human. |
| Local file | `--gmt FILE` | Repeatable. |
| Bundled test set | `--library demo` | Works offline. |

**Species.** Enrichr libraries are mostly human. For mouse data, or to use MSigDB directly, choose an MSigDB collection with `--msigdb NAME --species mouse|human` (downloaded once from the Broad release server, release 2025.1): `Hallmark`, `GO_BP`, `GO_CC`, `GO_MF`, `Reactome`, `WikiPathways`, `BioCarta`, `TFT`, `CellType`, `All` exist for both species; `KEGG`, `PID`, `HPO`, `Oncogenic`, `ImmuneSigDB` and `Vaccine` for human only (`snv libraries` prints the table). Mouse collections are MSigDB's native mouse sets, not orthology-mapped human ones. `--species` also sets the STRING taxon for `--gene-network string`.
