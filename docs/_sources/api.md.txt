# API reference

Everything on this page is available from the top-level package:

```python
import simpnetworkviz as snv
```

For an introduction with examples, see [Using the Python API](guide/python.md).

## Workflow

```{eval-rst}
.. autofunction:: simpnetworkviz.run_pipeline
.. autofunction:: simpnetworkviz.analyze_genes
```

## Gene sets and gene lists

```{eval-rst}
.. autofunction:: simpnetworkviz.load_gene_sets
.. autofunction:: simpnetworkviz.fetch_msigdb
.. autofunction:: simpnetworkviz.fetch_enrichr_library
.. autofunction:: simpnetworkviz.read_gene_list
.. autofunction:: simpnetworkviz.read_gene_groups
.. autofunction:: simpnetworkviz.read_gmt
.. autofunction:: simpnetworkviz.write_gmt
.. autofunction:: simpnetworkviz.normalize_species
.. autofunction:: simpnetworkviz.species_taxon
```

## Enrichment

```{eval-rst}
.. autofunction:: simpnetworkviz.enrich
.. autofunction:: simpnetworkviz.enrich_updown
.. autofunction:: simpnetworkviz.gsea_prerank
.. autofunction:: simpnetworkviz.read_ranking
.. autofunction:: simpnetworkviz.bh_fdr
.. autofunction:: simpnetworkviz.read_enrichment
.. autofunction:: simpnetworkviz.read_pathway_list
.. autofunction:: simpnetworkviz.table_from_pathways
```

## Comparing conditions and single-cell markers

```{eval-rst}
.. autofunction:: simpnetworkviz.compare_conditions
.. autofunction:: simpnetworkviz.merge_conditions
.. autofunction:: simpnetworkviz.enrichment_matrix
.. autofunction:: simpnetworkviz.markers_to_conditions
.. autofunction:: simpnetworkviz.read_marker_table
```

## Networks

```{eval-rst}
.. autofunction:: simpnetworkviz.build_pathway_network
.. autofunction:: simpnetworkviz.build_library_network
.. autofunction:: simpnetworkviz.build_gene_network
.. autoclass:: simpnetworkviz.Network
   :members: node_table, edge_table, cluster_table, plot, to_html, to_cytoscape, to_graphml, save
.. autofunction:: simpnetworkviz.annotate_cluster
.. autofunction:: simpnetworkviz.clean_pathway_name
.. autofunction:: simpnetworkviz.group_by_category
.. autofunction:: simpnetworkviz.group_by_keywords
.. autofunction:: simpnetworkviz.pathway_category
```

## Pathway maps

```{eval-rst}
.. autofunction:: simpnetworkviz.build_pathway_map
.. autofunction:: simpnetworkviz.save_map
.. autofunction:: simpnetworkviz.load_map
.. autofunction:: simpnetworkviz.plot_pathway_map
.. autofunction:: simpnetworkviz.map_hits
.. autofunction:: simpnetworkviz.region_summary
.. autofunction:: simpnetworkviz.write_map_html
.. autofunction:: simpnetworkviz.split_by_category
.. autofunction:: simpnetworkviz.build_category_maps
.. autofunction:: simpnetworkviz.save_category_maps
.. autofunction:: simpnetworkviz.plot_category_maps
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
