"""Validation against independent implementations and live web services.

The web tests are skipped automatically when the machine is offline.
"""
import socket

import numpy as np
import pandas as pd
import pytest
from scipy.stats import fisher_exact

import simpnetworkviz as snv
from simpnetworkviz.genesets import demo_path
from simpnetworkviz.network import fetch_string_edges, similarity


def _online(host="maayanlab.cloud", port=443, timeout=5) -> bool:
    try:
        socket.create_connection((host, port), timeout=timeout).close()
        return True
    except OSError:
        return False


needs_internet = pytest.mark.skipif(not _online(), reason="no internet connection")


@pytest.fixture(scope="module")
def gene_sets():
    return snv.load_gene_sets("demo")


@pytest.fixture(scope="module")
def gene_values():
    return snv.read_gene_list(demo_path("demo_genes.tsv"))


# ---------------------------------------------------------------- enrichment statistics

def test_enrich_matches_fisher_exact(gene_values, gene_sets):
    """Every p-value equals a one-sided Fisher exact test on the 2x2 table."""
    genes = [g for g, v in gene_values.items() if v is not None and v > 0]
    universe = sorted(set().union(*map(set, gene_sets.values())))
    res = snv.enrich(genes, gene_sets, min_size=1, min_overlap=1)
    q = set(genes) & set(universe)
    for _, r in res.iterrows():
        members = set(gene_sets[r["pathway"]]) & set(universe)
        a = len(q & members)
        b = len(members) - a
        c = len(q) - a
        d = len(universe) - a - b - c
        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        assert r["pvalue"] == pytest.approx(p, rel=1e-9)
        assert r["overlap"] == a and r["set_size"] == len(members)


def test_enrich_matches_gseapy(gene_values, gene_sets):
    """Same p-values as gseapy.enrich when the background is given explicitly."""
    gseapy = pytest.importorskip("gseapy")
    genes = [g for g, v in gene_values.items() if v is not None and v > 0]
    background = sorted(set().union(*map(set, gene_sets.values())) | {f"BG{i}" for i in range(500)})
    mine = snv.enrich(genes, gene_sets, background=background, min_size=1, min_overlap=1).set_index("pathway")
    ref = gseapy.enrich(gene_list=genes, gene_sets={k: list(v) for k, v in gene_sets.items()}, background=background,
                        outdir=None, verbose=False).results.set_index("Term")
    common = mine.index.intersection(ref.index)
    assert len(common) >= 10
    np.testing.assert_allclose(mine.loc[common, "pvalue"].values, ref.loc[common, "P-value"].astype(float).values, rtol=1e-6)
    overlap_ref = ref.loc[common, "Overlap"].str.split("/").str[0].astype(int).values
    np.testing.assert_array_equal(mine.loc[common, "overlap"].values, overlap_ref)


def test_bh_matches_statsmodels_or_reference():
    p = np.array([0.001, 0.01, 0.02, 0.03, 0.2, 0.5, 0.9, 0.04, 0.05, 0.3])
    try:
        from statsmodels.stats.multitest import multipletests
        ref = multipletests(p, method="fdr_bh")[1]
    except ImportError:
        # R: p.adjust(p, "BH")
        ref = np.array([0.01, 0.05, 0.0666667, 0.0714286, 0.2857143, 0.625, 0.9, 0.0714286, 0.0714286, 0.375])
    np.testing.assert_allclose(snv.bh_fdr(p), ref, rtol=1e-5)


def test_similarity_matches_enrichmentmap_definition():
    """EnrichmentMap: combined = k * overlap + (1 - k) * Jaccard with k = 0.5."""
    a, b = set("ABCDEFGH"), set("EFGHIJ")
    inter, union = 4, 10
    jac, ovl = inter / union, inter / min(len(a), len(b))
    assert similarity(a, b, "jaccard") == pytest.approx(jac)
    assert similarity(a, b, "overlap") == pytest.approx(ovl)
    assert similarity(a, b, "combined") == pytest.approx(0.5 * jac + 0.5 * ovl)
    assert similarity(a, set("XYZ")) == 0.0


def test_demo_results_are_stable(gene_values, gene_sets):
    """The numbers quoted in the README for `snv demo`."""
    table = snv.analyze_genes(gene_values, gene_sets)
    assert (table["fdr"] <= 0.05).sum() == 17 and len(table) == 20
    top = table.iloc[0]
    assert top["pathway"] == "HALLMARK_G2M_CHECKPOINT" and top["overlap"] == 23 and top["set_size"] == 30
    assert top["fdr"] == pytest.approx(1.8e-17, rel=0.1)
    net = snv.build_pathway_network(table, gene_sets, similarity_cutoff=0.3)
    assert net.graph.number_of_nodes() == 17 and net.graph.number_of_edges() == 16 and len(net.group_labels) == 5


# ---------------------------------------------------------------- live services

@needs_internet
def test_enrichr_library_download(tmp_path, monkeypatch):
    monkeypatch.setenv("SIMPNETWORKVIZ_CACHE", str(tmp_path))
    gs = snv.fetch_enrichr_library("MSigDB_Hallmark_2020")
    assert len(gs) == 50
    e2f = next(k for k in gs if "E2F" in k.upper())
    assert {"MCM2", "PCNA", "CDK1"} <= set(gs[e2f]) and 150 <= len(gs[e2f]) <= 250
    assert (tmp_path / "MSigDB_Hallmark_2020.gmt").exists()
    again = snv.load_gene_sets("MSigDB_Hallmark_2020")  # served from cache
    assert again == gs
    with pytest.raises(Exception):
        snv.fetch_enrichr_library("Not_A_Real_Library_XYZ")


@needs_internet
def test_enrichr_hallmark_enrichment_recovers_demo_signal(gene_values, tmp_path, monkeypatch):
    monkeypatch.setenv("SIMPNETWORKVIZ_CACHE", str(tmp_path))
    table = snv.analyze_genes(gene_values, "MSigDB_Hallmark_2020")
    sig = table[table["fdr"] <= 0.05]
    names = sig["pathway"].str.upper().str.replace(" ", "_")
    assert names.str.contains("E2F").any() and names.str.contains("G2-M|G2M").any()
    assert names.str.contains("OXIDATIVE").any()


@needs_internet
def test_string_api_returns_known_interactions():
    edges = fetch_string_edges(["CDK1", "CCNB1", "PLK1", "AURKA", "TP53", "MDM2"], species=9606, required_score=700)
    pairs = {frozenset(e[:2]) for e in edges}
    assert frozenset({"CDK1", "CCNB1"}) in pairs and frozenset({"TP53", "MDM2"}) in pairs
    assert all(0 <= s <= 1 for _, _, s in edges)


@needs_internet
def test_string_gene_network_end_to_end(gene_values, gene_sets):
    table = snv.analyze_genes(gene_values, gene_sets)
    pnet = snv.build_pathway_network(table, gene_sets, similarity_cutoff=0.3)
    gnet = snv.build_gene_network(pnet, gene_sets, mode="string", gene_values=gene_values, max_genes=40)
    assert gnet.graph.number_of_edges() > 10
    assert all("string_score" in d for _, _, d in gnet.graph.edges(data=True))
