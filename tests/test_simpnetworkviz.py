import json
import os

import numpy as np
import pandas as pd
import pytest
from scipy.stats import hypergeom

import simpnetworkviz as snv
from simpnetworkviz.cli import main
from simpnetworkviz.genesets import demo_path, parse_gmt_lines


@pytest.fixture(scope="module")
def gene_sets():
    return snv.load_gene_sets("demo")


@pytest.fixture(scope="module")
def gene_values():
    return snv.read_gene_list(demo_path("demo_genes.tsv"))


@pytest.fixture(scope="module")
def pnet(gene_sets, gene_values):
    table = snv.analyze_genes(gene_values, gene_sets)
    return snv.build_pathway_network(table, gene_sets, similarity_cutoff=0.3)


# ---------------------------------------------------------------- gene sets / input

def test_parse_gmt_handles_enrichr_weights():
    gs = parse_gmt_lines(["SET_A\t\tTP53,1.0\tmdm2\tTP53", "", "BAD_LINE"])
    assert gs == {"SET_A": ["TP53", "MDM2"]}


def test_read_gene_list_formats(tmp_path):
    plain = tmp_path / "a.txt"
    plain.write_text("tp53\nMDM2\n\n# comment\nCDKN1A\n")
    assert list(snv.read_gene_list(plain)) == ["TP53", "MDM2", "CDKN1A"]
    table = tmp_path / "b.csv"
    table.write_text("symbol,baseMean,log2FoldChange\nTP53,10,1.5\nMDM2,5,-2\n")
    assert snv.read_gene_list(table) == {"TP53": 1.5, "MDM2": -2.0}
    headerless = tmp_path / "c.tsv"
    headerless.write_text("TP53\t1.5\nMDM2\t-2\n")
    assert snv.read_gene_list(headerless) == {"TP53": 1.5, "MDM2": -2.0}


# ---------------------------------------------------------------- enrichment

def test_bh_fdr_matches_reference():
    p = np.array([0.01, 0.04, 0.03, 0.2])
    # R: p.adjust(c(.01,.04,.03,.2), "BH") = 0.04 0.05333 0.05333 0.2
    np.testing.assert_allclose(snv.bh_fdr(p), [0.04, 0.0533333, 0.0533333, 0.2], rtol=1e-5)


def test_enrich_hypergeometric():
    gs = {"A": [f"G{i}" for i in range(10)], "B": [f"G{i}" for i in range(10, 30)]}
    background = [f"G{i}" for i in range(100)]
    res = snv.enrich(["G0", "G1", "G2", "G3", "G50"], gs, background=background, min_size=1)
    row = res.set_index("pathway").loc["A"]
    assert row["overlap"] == 4 and row["set_size"] == 10
    assert row["pvalue"] == pytest.approx(hypergeom.sf(3, 100, 10, 5))
    assert row["genes"] == "G0;G1;G2;G3"
    assert "B" not in res["pathway"].values  # no overlap


def test_updown_split(gene_sets, gene_values):
    res = snv.analyze_genes(gene_values, gene_sets)
    assert set(res["direction"]) == {"up", "down"}
    top = res.set_index(["pathway", "direction"])
    assert ("HALLMARK_E2F_TARGETS", "up") in top.index
    assert ("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "down") in top.index


def test_read_enrichment_clusterprofiler_and_enrichr(tmp_path):
    cp = pd.DataFrame({"ID": ["GO:1"], "Description": ["cell cycle"], "GeneRatio": ["3/50"], "BgRatio": ["100/20000"],
                       "pvalue": [1e-5], "p.adjust": [1e-3], "qvalue": [1e-3], "geneID": ["CDK1/CCNB1/PLK1"], "Count": [3]})
    out = snv.read_enrichment(cp)
    assert out.loc[0, "pathway"] == "cell cycle" and out.loc[0, "fdr"] == 1e-3
    assert out.loc[0, "genes"] == "CDK1;CCNB1;PLK1" and out.loc[0, "set_size"] == 100

    er = pd.DataFrame({"Term": ["E2F Targets"], "Overlap": ["5/200"], "P-value": [1e-6], "Adjusted P-value": [1e-4],
                       "Combined Score": [50.0], "Genes": ["MCM2;MCM3;PCNA;CDK1;E2F1"]})
    out = snv.read_enrichment(er)
    assert out.loc[0, "overlap"] == 5 and out.loc[0, "set_size"] == 200 and out.loc[0, "fdr"] == 1e-4
    assert out.loc[0, "direction"] == ""  # combined score must not be read as a direction

    gsea = pd.DataFrame({"NAME": ["X", "Y"], "SIZE": [20, 30], "NES": [2.1, -1.8], "FDR q-val": [0.01, 0.02]})
    out = snv.read_enrichment(gsea).set_index("pathway")
    assert out.loc["X", "direction"] == "up" and out.loc["Y", "direction"] == "down"


def test_table_from_pathways_fuzzy_names(gene_sets):
    with pytest.warns(UserWarning):
        t = snv.table_from_pathways(["G2M checkpoint", "hallmark e2f targets", "nope"], gene_sets)
    assert list(t["pathway"]) == ["HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS"]


# ---------------------------------------------------------------- network

def test_clean_names_and_annotation():
    assert snv.clean_pathway_name("HALLMARK_INTERFERON_GAMMA_RESPONSE") == "Interferon gamma response"
    assert snv.clean_pathway_name("HALLMARK_TNFA_SIGNALING_VIA_NFKB") == "TNFA signaling via NFKB"
    assert snv.clean_pathway_name("cell cycle (GO:0007049)") == "Cell cycle"
    assert snv.clean_pathway_name("Cell Cycle R-HSA-1640170") == "Cell Cycle"
    label = snv.annotate_cluster(["HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                                  "REACTOME_INTERFERON_SIGNALING"])
    assert label.lower().startswith("interferon")


def test_pathway_network_structure(pnet):
    G = pnet.graph
    assert G.number_of_nodes() > 10 and G.number_of_edges() > 5
    assert set(pnet.pos) == set(G.nodes)
    # interferon pathways cluster together, apart from ECM pathways
    g = pnet.groups
    ifn = {g["HALLMARK_INTERFERON_ALPHA_RESPONSE"], g["HALLMARK_INTERFERON_GAMMA_RESPONSE"], g["REACTOME_INTERFERON_SIGNALING"]}
    assert len(ifn) == 1 and 0 not in ifn
    assert g["REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"] not in ifn
    ct = pnet.cluster_table()
    assert list(ct["cluster"]) == list(range(1, len(ct) + 1))


def test_layout_has_no_overlapping_nodes(pnet):
    nodes = list(pnet.graph)
    P = np.array([pnet.pos[n] for n in nodes])
    R = np.array([pnet.radius[n] for n in nodes])
    d = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1))
    np.fill_diagonal(d, np.inf)
    assert (d >= (R[:, None] + R[None]) * 0.999).all()


def test_mcl_separates_two_cliques():
    import networkx as nx
    G = nx.disjoint_union(nx.complete_graph(4), nx.complete_graph(4))
    G.add_edge(0, 4)
    nx.set_edge_attributes(G, 1.0, "weight")
    clusters = snv.network.mcl_clusters(G)
    assert sorted(map(sorted, clusters)) == [[0, 1, 2, 3], [4, 5, 6, 7]]


@pytest.mark.parametrize("mode", ["shared", "cnet"])
def test_gene_networks(pnet, gene_sets, gene_values, mode):
    gnet = snv.build_gene_network(pnet, gene_sets, mode=mode, gene_values=gene_values, max_genes=40)
    genes = [n for n, d in gnet.graph.nodes(data=True) if d["type"] == "gene"]
    assert 0 < len(genes) <= 40
    assert all(gnet.graph.nodes[g]["logfc"] is not None for g in genes)
    if mode == "cnet":
        assert any(d["type"] == "pathway" for _, d in gnet.graph.nodes(data=True))


def test_string_mode_uses_api(monkeypatch, pnet, gene_sets):
    calls = {}

    def fake(genes, species, required_score):
        calls["n"] = len(genes)
        return [(genes[0], genes[1], 0.9)]

    monkeypatch.setattr(snv.network, "fetch_string_edges", fake)
    gnet = snv.build_gene_network(pnet, gene_sets, mode="string", max_genes=20)
    assert calls["n"] == 20 and gnet.graph.number_of_edges() == 1


# ---------------------------------------------------------------- output

def test_save_all_outputs(pnet, tmp_path):
    files = pnet.save(tmp_path, "pn", ("png", "svg", "html", "cyjs", "graphml", "csv"), title="t")
    for f in files.values():
        assert os.path.getsize(f) > 0
    doc = json.load(open(files["cyjs"]))
    assert len(doc["elements"]["nodes"]) == pnet.graph.number_of_nodes()
    assert all("position" in n for n in doc["elements"]["nodes"])
    html = open(files["html"]).read()
    assert "cytoscape" in html and "HALLMARK_E2F_TARGETS" in html


@pytest.mark.parametrize("color_by", ["cluster", "direction", "significance"])
def test_color_modes(pnet, color_by):
    import matplotlib.pyplot as plt
    fig = pnet.plot(color_by=color_by, node_labels="top")
    plt.close(fig)


def test_cli_run_and_pathways(tmp_path, capsys):
    out = tmp_path / "o"
    assert main(["run", demo_path("demo_genes.tsv"), "--library", "demo", "-o", str(out), "--cutoff", "0.3",
                 "--formats", "png,csv"]) == 0
    assert (out / "pathway_network.png").exists() and (out / "gene_network.png").exists()
    assert (out / "enrichment.csv").exists()
    pw = tmp_path / "pw.txt"
    pw.write_text("HALLMARK_E2F_TARGETS\nKEGG_CELL_CYCLE\nHALLMARK_G2M_CHECKPOINT\n")
    assert main(["run", "--pathways", str(pw), "--library", "demo", "-o", str(tmp_path / "p"), "--formats", "png",
                 "--gene-network", "cnet"]) == 0
    assert "Pathway network: 3 pathways" in capsys.readouterr().out


@pytest.mark.parametrize("arrangement", ["auto", "horizontal", "vertical"])
def test_pathway_gene_panel_figure(pnet, gene_sets, gene_values, tmp_path, arrangement):
    import matplotlib.pyplot as plt
    gnet = snv.build_gene_network(pnet, gene_sets, gene_values=gene_values, clusters=[1, 2])
    out = tmp_path / "panel.pdf"
    fig = snv.plot_pathway_gene_network(pnet, gnet, out, arrangement=arrangement)
    texts = {t.get_text() for t in fig.texts}
    assert {"a", "b", "Pathway network", "Gene network"} <= texts
    assert out.stat().st_size > 0
    plt.close(fig)


def test_panels_share_cluster_colors(pnet, gene_sets):
    from simpnetworkviz.style import hull_colors
    gnet = snv.build_gene_network(pnet, gene_sets)
    a = hull_colors(pnet, "direction", "cluster")
    b = hull_colors(gnet, "logfc", "cluster")
    for gid in set(a) & set(b):
        assert a[gid] == b[gid]


def test_pipeline_writes_separate_figures(tmp_path, gene_values):
    res = snv.run_pipeline(gene_values, gene_sets="demo", outdir=str(tmp_path), similarity_cutoff=0.3, formats=("png",))
    assert os.path.exists(res["files"]["pathway_network"]["png"])
    assert os.path.exists(res["files"]["gene_network"]["png"])
    assert "pathway_gene_network" not in res["files"]


def test_hull_encloses_node_labels(pnet):
    from matplotlib.path import Path
    from simpnetworkviz.plotting import _Panel, _label_box
    panel = _Panel(pnet, node_labels="all")
    for gid, poly in panel.hull_polys.items():
        path = Path(poly)
        for n in (n for n, g in pnet.groups.items() if g == gid):
            assert path.contains_points(_label_box(pnet, n, 9.0), radius=1e-6).all()
