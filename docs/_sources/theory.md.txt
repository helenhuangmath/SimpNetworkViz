# Theory behind SimpNetworkViz

**Enrichment.** Over-representation uses the hypergeometric test: for a pathway with *K* genes in a universe of *N*, and a list of *n* genes with *k* in the pathway, the p-value is P(X ≥ k). Every gene set of the library that falls within the size limits counts toward the Benjamini–Hochberg correction. Preranked GSEA follows Subramanian et al. (2005): the weighted Kolmogorov–Smirnov running sum (weight 1) gives the enrichment score; the null distribution comes from gene permutations, shared across sets of the same size; NES = ES / mean |null ES| of the same sign; empirical p-values are corrected with BH as in `fgsea` (`fdr_method="gsea"` in Python gives the original NES-based q-value).

**Similarity.** Two pathways are compared by their full gene sets with the combined coefficient of EnrichmentMap, ½·Jaccard + ½·overlap coefficient, and connected when it reaches 0.375. The overlap coefficient lets a small pathway attach to a large one that contains it; Jaccard tempers that for very unequal sizes.

**Clustering and labels.** Markov clustering (inflation 2) on the weighted similarity graph, as in AutoAnnotate's default. Each cluster is labelled with the informative words shared by at least two of its pathways, ordered as they appear in the names, or with its most significant pathway when it has fewer than three members. Clusters with a single pathway are left unlabelled.

**Layout.** Each cluster is laid out separately (Kamada–Kawai, then nodes pushed apart so that circles and labels never overlap) and rotated so its long axis is horizontal. Clusters are then packed in rows, each with room for its label; unclustered nodes go into a grid at the end. This mirrors AutoAnnotate's "layout network to prevent cluster overlap".

**Gene networks.** Genes are those of your list that hit the selected pathways (the leading edge for GSEA; the full gene sets when you start from a pathway list without genes). Each gene is assigned to the pathway cluster it belongs to most, so gene modules carry the same labels and colors as the pathway network. Edges are shared pathways (weight = number shared) or STRING interactions (score ≥ 0.4).

**Comparisons.** Each condition is tested separately. The merged node keeps the best FDR, the union of hit genes, and per-condition FDR, direction and genes. The direction of the node is `up` or `down` when all conditions agree and `mixed` otherwise.

**Whole library.** All pathways of the library are added before similarity and clustering, so the clusters describe the library; enriched pathways are highlighted.

**Validation.** `tests/test_validation.py` checks every ORA p-value against a one-sided Fisher exact test (scipy) and against `gseapy.enrich`; BH FDR against `statsmodels`/R `p.adjust`; GSEA enrichment scores and leading-edge genes against `gseapy.prerank` (identical); and the similarity formula against the EnrichmentMap definition. Live tests download `MSigDB_Hallmark_2020` from Enrichr and query STRING for known interactions, and are skipped when offline. Run them with `pip install -e ".[test]"` and `pytest -q`.

## References

* Merico D. et al. Enrichment Map: a network-based method for gene-set enrichment visualization and interpretation. *PLoS ONE* 5, e13984 (2010).
* Kucera M. et al. AutoAnnotate: a Cytoscape app for summarizing networks with semantic annotations. *F1000Research* 5, 1717 (2016).
* Enright A.J., Van Dongen S., Ouzounis C.A. An efficient algorithm for large-scale detection of protein families. *Nucleic Acids Res.* 30, 1575–1584 (2002).
* Subramanian A. et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 102, 15545–15550 (2005).
