# API reference

Everything on this page is available from the top-level package:

```python
import simpnetworkviz as snv
```

For a guided tour, see [Step-by-step in Python](tutorial/python.md).

## Workflow

```{eval-rst}
.. autofunction:: simpnetworkviz.run_pipeline
.. autofunction:: simpnetworkviz.analyze_genes
```

## Gene sets and gene lists

```{eval-rst}
.. autofunction:: simpnetworkviz.load_gene_sets
.. autofunction:: simpnetworkviz.read_gene_list
.. autofunction:: simpnetworkviz.read_gmt
.. autofunction:: simpnetworkviz.write_gmt
.. autofunction:: simpnetworkviz.fetch_enrichr_library
.. autofunction:: simpnetworkviz.fetch_msigdb
```

## Enrichment

```{eval-rst}
.. autofunction:: simpnetworkviz.enrich
.. autofunction:: simpnetworkviz.enrich_updown
.. autofunction:: simpnetworkviz.bh_fdr
.. autofunction:: simpnetworkviz.read_enrichment
.. autofunction:: simpnetworkviz.read_pathway_list
.. autofunction:: simpnetworkviz.table_from_pathways
.. autofunction:: simpnetworkviz.gsea_prerank
.. autofunction:: simpnetworkviz.compare_conditions
```

## Networks

```{eval-rst}
.. autofunction:: simpnetworkviz.build_pathway_network
.. autofunction:: simpnetworkviz.build_gene_network
.. autofunction:: simpnetworkviz.build_library_network
.. autoclass:: simpnetworkviz.Network
   :members: node_table, edge_table, cluster_table, plot, to_html, to_cytoscape, to_graphml, save
.. autofunction:: simpnetworkviz.annotate_cluster
.. autofunction:: simpnetworkviz.clean_pathway_name
```

## Plotting

```{eval-rst}
.. autofunction:: simpnetworkviz.plot_network
.. autofunction:: simpnetworkviz.plot_pathway_gene_network
.. autofunction:: simpnetworkviz.plot_cluster_barplot
```

## Export

```{eval-rst}
.. autofunction:: simpnetworkviz.save_all
.. autofunction:: simpnetworkviz.write_html
.. autofunction:: simpnetworkviz.write_cyjs
.. autofunction:: simpnetworkviz.write_graphml
```

## Pathway maps

```{eval-rst}
.. autofunction:: simpnetworkviz.build_pathway_map
.. autofunction:: simpnetworkviz.build_category_maps
.. autofunction:: simpnetworkviz.map_hits
.. autofunction:: simpnetworkviz.plot_pathway_map
.. autofunction:: simpnetworkviz.write_map_html
```
