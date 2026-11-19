"""Tests for comparison, GSEA, single-cell input, grouping options and library networks."""
import json
import os

import numpy as np
import pandas as pd
import pytest

import simpnetworkviz as snv
from simpnetworkviz.cli import main
from simpnetworkviz.genesets import demo_path
from simpnetworkviz.gsea import gsea_prerank
from simpnetworkviz.network import cluster_graph, pathway_category, read_category_file, read_theme_file
from simpnetworkviz.singlecell import markers_to_conditions, read_marker_table


@pytest.fixture(scope="module")
def gene_sets():
    return snv.load_gene_sets("demo")


@pytest.fixture(scope="module")
def gene_values():
    return snv.read_gene_list(demo_path("demo_genes.tsv"))


@pytest.fixture(scope="module")
def conditions(gene_values):
    up = {g: v for g, v in gene_values.items() if v is not None and v > 0}
    down = {g: v for g, v in gene_values.items() if v is not None and v < 0}
    return {"up": up, "down": down, "all": gene_values}


@pytest.fixture(scope="module")
def cmp_table(conditions, gene_sets):
    return snv.compare_conditions(conditions, gene_sets)


@pytest.fixture(scope="module")
def cmp_net(cmp_table, gene_sets):
    return snv.build_pathway_network(cmp_table, gene_sets, similarity_cutoff=0.3)


@pytest.fixture(scope="module")
def ranking(gene_sets):
    """Synthetic ranking: E2F targets pushed up, OXPHOS pushed down, 2000 noise genes."""
    rng = np.random.default_rng(1)
    genes = sorted(set().union(*map(set, gene_sets.values()))) + [f"NOISE{i}" for i in range(2000)]
    score = rng.normal(size=len(genes))
    idx = {g: i for i, g in enumerate(genes)}
    for g in gene_sets["HALLMARK_E2F_TARGETS"]:
        score[idx[g]] += 2.0
    for g in gene_sets["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]:
        score[idx[g]] -= 2.0
    return pd.Series(score, index=genes)


# ---------------------------------------------------------------- compare conditions

def test_compare_table_has_condition_column(cmp_table):
    assert set(cmp_table["condition"]) == {"up", "down", "all"}
    assert {"pathway", "fdr", "genes"} <= set(cmp_table.columns)
    sub = cmp_table[cmp_table["condition"] == "up"].set_index("pathway")
    assert sub.loc["HALLMARK_E2F_TARGETS", "fdr"] < 1e-10
    assert "HALLMARK_OXIDATIVE_PHOSPHORYLATION" not in sub.index or sub.loc["HALLMARK_OXIDATIVE_PHOSPHORYLATION", "fdr"] > 0.05


def test_merge_conditions(cmp_table):
    merged = snv.merge_conditions(cmp_table)
    assert merged["pathway"].is_unique
    row = merged.set_index("pathway").loc["HALLMARK_E2F_TARGETS"]
    assert row["fdr__up"] < 1e-10 and pd.isna(row["fdr__down"])
    assert row["conditions"] == "up;all" and row["n_conditions"] == 2
    assert row["direction"] == "up"
    assert row["fdr"] == row["fdr__up"] or row["fdr"] == row["fdr__all"]


def test_enrichment_matrix(cmp_table):
    mat = snv.enrichment_matrix(cmp_table)
    assert list(mat.columns) == ["up", "down", "all"]
    assert mat.loc["HALLMARK_OXIDATIVE_PHOSPHORYLATION", "all"] < 0  # signed by direction
    assert mat.loc["HALLMARK_E2F_TARGETS", "up"] > 10
    assert pd.isna(mat.loc["HALLMARK_E2F_TARGETS", "down"])


def test_comparison_network_pies(cmp_net, tmp_path):
    from simpnetworkviz.style import node_colors

    assert cmp_net.conditions == ["up", "down", "all"]
    d = cmp_net.graph.nodes["HALLMARK_E2F_TARGETS"]
    assert d["fdr__up"] is not None and d["fdr__down"] is None
    colors, leg = node_colors(cmp_net, "auto")
    assert leg["mode"] == "condition"
    w = leg["wedges"]["HALLMARK_E2F_TARGETS"]
    assert len(w) == 3 and w[0] and w[1] is None and w[2]
    # direction mode on a comparison network colors slices red/blue
    _, leg2 = node_colors(cmp_net, "direction")
    assert leg2["wedges"]["HALLMARK_E2F_TARGETS"][0] != leg2["wedges"]["HALLMARK_OXIDATIVE_PHOSPHORYLATION"][1]
    files = cmp_net.save(tmp_path, "cmp", ("png", "html", "cyjs", "csv"))
    doc = json.load(open(files["cyjs"]))
    pies = [n for n in doc["elements"]["nodes"] if n["data"].get("pie_n")]
    assert len(pies) == cmp_net.graph.number_of_nodes()
    assert all(n["data"]["pie_n"] == 3 and "pie3" in n["data"] for n in pies)
    html = open(files["html"]).read()
    assert "pie-1-background-color" in html or "pie-${i}-background-color" in html
    assert '"n_conditions": 3' in html


def test_comparison_gene_network(cmp_net, gene_sets):
    gnet = snv.build_gene_network(cmp_net, gene_sets, max_genes=30)
    assert gnet.conditions == cmp_net.conditions
    d = gnet.graph.nodes[next(iter(gnet.graph.nodes))]
    assert "conditions" in d and d["n_conditions"] >= 1


def test_compare_from_enrichment_tables(cmp_table, gene_sets):
    parts = {c: cmp_table[cmp_table["condition"] == c].drop(columns="condition") for c in ["up", "down"]}
    out = snv.compare_conditions(parts, gene_sets)
    assert set(out["condition"]) == {"up", "down"}


# ---------------------------------------------------------------- GSEA

def test_gsea_recovers_planted_signal(ranking, gene_sets):
    res = gsea_prerank(ranking, gene_sets, min_size=10, n_perm=500).set_index("pathway")
    assert res.loc["HALLMARK_E2F_TARGETS", "nes"] > 1.5 and res.loc["HALLMARK_E2F_TARGETS", "direction"] == "up"
    assert res.loc["HALLMARK_OXIDATIVE_PHOSPHORYLATION", "nes"] < -1.5
    assert res.loc["HALLMARK_E2F_TARGETS", "fdr"] < 0.05
    lead = res.loc["HALLMARK_E2F_TARGETS", "genes"].split(";")
    assert 0 < len(lead) <= 34 and set(lead) <= set(gene_sets["HALLMARK_E2F_TARGETS"])
    assert res.loc["HALLMARK_E2F_TARGETS", "overlap"] == len(lead)


def test_gsea_matches_gseapy(ranking, gene_sets):
    gseapy = pytest.importorskip("gseapy")
    mine = gsea_prerank(ranking, gene_sets, min_size=10, max_size=500, n_perm=1000).set_index("pathway")
    ref = gseapy.prerank(rnk=ranking, gene_sets={k: list(v) for k, v in gene_sets.items()}, min_size=10, max_size=500,
                         permutation_num=1000, outdir=None, seed=0, threads=1, verbose=False).res2d.set_index("Term")
    ref = ref.loc[mine.index]
    np.testing.assert_allclose(mine["es"].values, ref["ES"].astype(float).values, atol=1e-9)
    np.testing.assert_array_equal(mine["overlap"].values, (ref["Lead_genes"].str.count(";") + 1).values)
    # NES depends on the permutation null; agree to ~10% for the planted signals
    for pw in ["HALLMARK_E2F_TARGETS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"]:
        assert mine.loc[pw, "nes"] == pytest.approx(float(ref.loc[pw, "NES"]), rel=0.15)
    assert (mine["direction"] == np.where(ref["NES"].astype(float) > 0, "up", "down")).all()


def test_gsea_network_and_pipeline(ranking, gene_sets, tmp_path):
    table = gsea_prerank(ranking, gene_sets, min_size=10, n_perm=200)
    pnet = snv.build_pathway_network(table, gene_sets, fdr_cutoff=0.25, similarity_cutoff=0.3)
    from simpnetworkviz.style import resolve_color_by
    assert resolve_color_by(pnet) == "direction"
    assert pnet.graph.nodes["HALLMARK_E2F_TARGETS"]["nes"] > 0
    res = snv.run_pipeline(ranking.to_dict(), gene_sets=gene_sets, outdir=str(tmp_path), method="gsea", fdr_cutoff=0.25,
                           similarity_cutoff=0.3, formats=("png", "csv"), enrich_kwargs={"min_size": 10, "n_perm": 200})
    assert "nes" in res["enrichment"].columns
    genes = {n for n, d in res["gene_network"].graph.nodes(data=True) if d["type"] == "gene"}
    assert not any(g.startswith("NOISE") for g in genes)  # only leading-edge genes
    assert (tmp_path / "run_params.json").exists()


def test_gsea_needs_scores(gene_sets):
    with pytest.raises(ValueError):
        snv.analyze_genes({"MCM2": None, "MCM3": None}, gene_sets, method="gsea")


# ---------------------------------------------------------------- single cell

def _seurat_markers(gene_values):
    rows = []
    up = [g for g, v in gene_values.items() if v is not None and v > 0]
    down = [g for g, v in gene_values.items() if v is not None and v < 0]
    for cl, genes in [("0", up[:30]), ("1", up[30:] + ["NOT_A_GENE"]), ("T cells", down)]:
        for i, g in enumerate(genes):
            rows.append({"": g, "p_val": 1e-10, "avg_log2FC": 1.0 + i * 0.01, "pct.1": 0.9, "pct.2": 0.1,
                         "p_val_adj": 1e-8 if i < len(genes) - 1 else 0.5, "cluster": cl, "gene": g})
    return pd.DataFrame(rows)


def test_read_seurat_markers(gene_values, tmp_path):
    df = _seurat_markers(gene_values)
    p = tmp_path / "markers.csv"
    df.to_csv(p, index=False)
    t = read_marker_table(p)
    assert list(t.columns) == ["cluster", "gene", "logfc", "padj", "score"]
    assert set(t["cluster"]) == {"0", "1", "T cells"}
    conds = markers_to_conditions(p)
    assert set(conds) == {"0", "1", "T cells"}
    assert len(conds["0"]) == 29  # last marker fails padj
    assert all(v > 0 for v in conds["0"].values())
    top = markers_to_conditions(p, top_n=5)
    assert all(len(v) == 5 for v in top.values())
    sub = markers_to_conditions(p, clusters=["1"])
    assert list(sub) == ["1"]
    with pytest.raises(ValueError):
        markers_to_conditions(p, clusters=["99"])


def test_read_scanpy_markers():
    df = pd.DataFrame({"group": ["a", "a", "b"], "names": ["Mcm2", "Mcm3", "Col1a1"], "scores": [5.0, 4.0, 3.0],
                       "logfoldchanges": [2.0, 1.5, -1.0], "pvals": [1e-5] * 3, "pvals_adj": [1e-4] * 3})
    conds = markers_to_conditions(df, direction="both")
    assert conds == {"a": {"MCM2": 2.0, "MCM3": 1.5}, "b": {"COL1A1": -1.0}}
    wide = pd.DataFrame({"0_n": ["MCM2", "MCM3"], "0_s": [5.0, 4.0], "0_l": [2.0, 1.5], "0_p": [1e-4, 1e-4],
                         "1_n": ["COL1A1", "COL1A2"], "1_s": [3.0, 2.0], "1_l": [1.0, 0.5], "1_p": [1e-4, 1e-4]})
    t = read_marker_table(wide)
    assert set(t["cluster"]) == {"0", "1"} and len(t) == 4


def test_markers_cli_whole_library(gene_values, tmp_path, capsys):
    df = _seurat_markers(gene_values)
    p = tmp_path / "markers.csv"
    df.to_csv(p, index=False)
    out = tmp_path / "o"
    assert main(["markers", str(p), "--library", "demo", "-o", str(out), "--whole-library", "--cutoff", "0.3",
                 "--formats", "png,csv", "--min-logfc", "0"]) == 0
    text = capsys.readouterr().out
    assert "3 clusters" in text and "Conditions (pie slices): 0, 1, T cells" in text
    nodes = pd.read_csv(out / "pathway_network_nodes.csv")
    assert len(nodes) == 32  # every demo gene set is a node
    assert nodes["highlight"].sum() > 0 and (~nodes["highlight"]).sum() > 0
    assert (out / "enrichment_matrix.csv").exists()


# ---------------------------------------------------------------- library network

def test_library_network_highlights(gene_values, gene_sets):
    table = snv.analyze_genes(gene_values, gene_sets)
    net = snv.build_library_network(table, gene_sets, similarity_cutoff=0.3)
    assert net.graph.number_of_nodes() == len(gene_sets)
    hl = {n for n, d in net.graph.nodes(data=True) if d["highlight"]}
    assert "HALLMARK_E2F_TARGETS" in hl and "HALLMARK_DNA_REPAIR" not in hl
    from simpnetworkviz.style import CONTEXT, hull_colors, node_colors
    colors, _ = node_colors(net, "direction")
    assert colors["HALLMARK_DNA_REPAIR"] == CONTEXT and colors["HALLMARK_E2F_TARGETS"] != CONTEXT
    hc = hull_colors(net, "direction")
    grey = {gid for gid, (fill, edge) in hc.items() if edge == "#B7B2A6"}
    lit = {net.groups[n] for n in hl} - {0}
    assert grey.isdisjoint(lit)
    with pytest.warns(UserWarning):
        small = snv.build_library_network(table, gene_sets, max_pathways=20, similarity_cutoff=0.3)
    assert small.graph.number_of_nodes() == 20 and all(d["highlight"] for n, d in small.graph.nodes(data=True) if n in hl)


# ---------------------------------------------------------------- grouping options

@pytest.mark.parametrize("method", ["louvain", "greedy", "label_propagation", "hierarchical", "components"])
def test_cluster_methods(gene_values, gene_sets, method):
    table = snv.analyze_genes(gene_values, gene_sets)
    net = snv.build_pathway_network(table, gene_sets, similarity_cutoff=0.3, cluster_method=method)
    g = net.groups
    assert g["HALLMARK_INTERFERON_ALPHA_RESPONSE"] == g["HALLMARK_INTERFERON_GAMMA_RESPONSE"] != 0
    assert g["HALLMARK_INTERFERON_ALPHA_RESPONSE"] != g["REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"]


def test_pathway_category():
    assert pathway_category("HALLMARK_E2F_TARGETS") == "Hallmark"
    assert pathway_category("KEGG_CELL_CYCLE") == "KEGG"
    assert pathway_category("GOBP_DNA_REPLICATION") == "GO biological process"
    assert pathway_category("cell cycle (GO:0007049)") == "GO"
    assert pathway_category("Cell Cycle R-HSA-1640170") == "Reactome"
    assert pathway_category("Reactome_2022: Cell Cycle R-HSA-1640170") == "Reactome_2022"
    assert pathway_category("my pathway", {"my pathway": "Custom"}) == "Custom"
    assert pathway_category("whatever") == "Other"


def test_group_by_category_and_keyword(gene_values, gene_sets, tmp_path):
    table = snv.analyze_genes(gene_values, gene_sets)
    net = snv.build_pathway_network(table, gene_sets, similarity_cutoff=0.3, cluster_method="category")
    labels = set(net.group_labels.values())
    assert {"Hallmark", "KEGG", "Reactome"} <= labels
    assert all(net.graph.nodes[n]["cluster_label"] == "Hallmark" for n in net.graph if n.startswith("HALLMARK"))

    themes = {"Proliferation": ["cycle", "E2F", "mitotic", "replication", "checkpoint"],
              "Immune": ["interferon"], "Matrix": ["matrix", "collagen", "mesenchymal", "adhesion", "ECM"]}
    net = snv.build_pathway_network(table, gene_sets, similarity_cutoff=0.3, cluster_method="keyword", themes=themes)
    by_label = {}
    for n in net.graph:
        by_label.setdefault(net.graph.nodes[n]["cluster_label"], set()).add(n)
    assert {"HALLMARK_E2F_TARGETS", "KEGG_CELL_CYCLE", "HALLMARK_G2M_CHECKPOINT"} <= by_label["Proliferation"]
    assert "REACTOME_INTERFERON_SIGNALING" in by_label["Immune"]
    assert net.groups["HALLMARK_OXIDATIVE_PHOSPHORYLATION"] == 0  # matches no theme -> unclustered

    f = tmp_path / "themes.tsv"
    f.write_text("# theme\tkeywords\nImmune\tinterferon, cytokine\nMatrix\tmatrix;collagen\n")
    assert read_theme_file(f) == {"Immune": ["interferon", "cytokine"], "Matrix": ["matrix", "collagen"]}
    c = tmp_path / "cats.tsv"
    c.write_text("pathway\tcategory\nHALLMARK_E2F_TARGETS\tCell cycle\n")
    assert read_category_file(c) == {"HALLMARK_E2F_TARGETS": "Cell cycle"}
    with pytest.raises(ValueError):
        snv.build_pathway_network(table, gene_sets, cluster_method="keyword")


def test_similarity_on_hit_genes(gene_values, gene_sets):
    table = snv.analyze_genes(gene_values, gene_sets)
    full = snv.build_pathway_network(table, gene_sets, similarity_cutoff=0.3)
    hits = snv.build_pathway_network(table, gene_sets, similarity_cutoff=0.3, use_full_gene_sets=False)
    assert full.graph.number_of_edges() != hits.graph.number_of_edges() or \
        any(full.graph.edges[e]["weight"] != hits.graph.edges[e]["weight"] for e in full.graph.edges if hits.graph.has_edge(*e))


def test_compare_cli(gene_values, tmp_path):
    up = tmp_path / "up.tsv"
    down = tmp_path / "down.tsv"
    up.write_text("gene\tlog2FC\n" + "\n".join(f"{g}\t{v}" for g, v in gene_values.items() if v and v > 0))
    down.write_text("gene\tlog2FC\n" + "\n".join(f"{g}\t{v}" for g, v in gene_values.items() if v and v < 0))
    out = tmp_path / "o"
    assert main(["compare", f"KO={up}", f"WT={down}", "--library", "demo", "-o", str(out), "--cutoff", "0.3",
                 "--formats", "png,csv", "--cluster", "louvain", "--resolution", "1.2"]) == 0
    nodes = pd.read_csv(out / "pathway_network_nodes.csv")
    assert {"fdr__KO", "fdr__WT"} <= set(nodes.columns)
    params = json.load(open(out / "run_params.json"))
    assert params["cluster_method"] == "louvain" and params["cluster_kwargs"]["louvain_resolution"] == "1.2"


# ---------------------------------------------------------------- wide gene-group tables

def test_read_gene_groups_and_compare_cli(gene_values, tmp_path):
    up = [g for g, v in gene_values.items() if v and v > 0]
    down = [g for g, v in gene_values.items() if v and v < 0]
    rows = max(len(up), len(down))
    wide = pd.DataFrame({"Up": up + [""] * (rows - len(up)), "Down": down + [""] * (rows - len(down))})
    p = tmp_path / "groups.csv"
    wide.to_csv(p, index=False)
    groups = snv.read_gene_groups(p)
    assert list(groups) == ["Up", "Down"] and len(groups["Up"]) == len(up) and len(groups["Down"]) == len(down)
    out = tmp_path / "o"
    assert main(["compare", "--groups", str(p), "--library", "demo", "-o", str(out), "--cutoff", "0.3",
                 "--formats", "pdf,csv"]) == 0
    nodes = pd.read_csv(out / "pathway_network_nodes.csv")
    assert {"fdr__Up", "fdr__Down"} <= set(nodes.columns)
    assert (out / "pathway_network.pdf").exists() and not (out / "pathway_network.png").exists()


def test_library_network_keeps_conditions_without_hits(gene_values, gene_sets):
    up = {g: v for g, v in gene_values.items() if v and v > 0}
    table = snv.compare_conditions({"Up": up, "Nothing": {"NOT_A_GENE_1": None, "NOT_A_GENE_2": None}}, gene_sets)
    assert snv.compare.condition_names(table) == ["Up", "Nothing"]
    net = snv.build_library_network(table, gene_sets, similarity_cutoff=0.3)
    assert net.conditions == ["Up", "Nothing"]
    from simpnetworkviz.style import node_colors
    _, leg = node_colors(net, "condition")
    w = leg["wedges"]["HALLMARK_E2F_TARGETS"]
    assert len(w) == 2 and w[0] and w[1] is None


# ---------------------------------------------------------------- pathway map

@pytest.fixture(scope="module")
def pmap(gene_sets):
    pytest.importorskip("sklearn")
    return snv.build_pathway_map(gene_sets, method="tsne", min_size=5, max_size=500, n_neighbors=5, min_cluster_label=2)


def test_pathway_map_build_and_roundtrip(pmap, gene_sets, tmp_path):
    assert len(pmap) == len(gene_sets) and set(pmap.columns) == set(snv.pathmap.MAP_COLUMNS)
    assert pmap.attrs["method"] == "tsne" and pmap[["x", "y"]].notna().all().all()
    assert pmap["cluster"].max() >= 2 and (pmap.loc[pmap["cluster"] > 0, "cluster_label"] != "").all()
    p = tmp_path / "map.csv"
    snv.save_map(pmap, p)
    back = snv.load_map(p)
    assert back.attrs["method"] == "tsne" and list(back["pathway"]) == list(pmap["pathway"])


def test_pathway_map_overlay(pmap, gene_sets, gene_values, tmp_path):
    import matplotlib.pyplot as plt
    table = snv.analyze_genes(gene_values, gene_sets)
    hits = snv.map_hits(pmap, table)
    assert set(hits["pathway"]) == set(table.loc[table["fdr"] <= 0.05, "pathway"])
    assert {"x", "y", "cluster_label"} <= set(hits.columns)
    regions = snv.region_summary(pmap, hits)
    assert regions["n_hits"].sum() == len(hits[hits["cluster"] > 0])
    fig = snv.plot_pathway_map(pmap, table, tmp_path / "map.pdf")
    plt.close(fig)
    assert (tmp_path / "map.pdf").stat().st_size > 0
    # comparison: pies and per-condition region counts
    up = {g: v for g, v in gene_values.items() if v and v > 0}
    down = {g: v for g, v in gene_values.items() if v and v < 0}
    cmp_table = snv.compare_conditions({"up": up, "down": down}, gene_sets)
    hits = snv.map_hits(pmap, cmp_table)
    regions = snv.region_summary(pmap, hits)
    assert {"n_up", "n_down", "fraction_up"} <= set(regions.columns)
    fig = snv.plot_pathway_map(pmap, cmp_table, tmp_path / "cmp.pdf")
    plt.close(fig)


def test_map_cli(pmap, gene_values, tmp_path, capsys):
    pytest.importorskip("sklearn")
    m = tmp_path / "demo_map.csv"
    assert main(["map", "--library", "demo", "-o", str(m), "--method", "tsne", "--min-size", "5", "--neighbors", "5"]) == 0
    assert m.exists() and (tmp_path / "demo_map.pdf").exists()
    assert "labelled regions" in capsys.readouterr().out
    out = tmp_path / "o"
    assert main(["run", demo_path("demo_genes.tsv"), "--library", "demo", "-o", str(out), "--cutoff", "0.3",
                 "--formats", "pdf,csv", "--map", str(m), "--gene-network", "none"]) == 0
    assert (out / "pathway_map.pdf").exists() and (out / "pathway_map_regions.csv").exists()
    hits = pd.read_csv(out / "pathway_map_hits.csv")
    assert "HALLMARK_E2F_TARGETS" in set(hits["pathway"])


# ---------------------------------------------------------------- maps per category

def test_split_by_category(gene_sets):
    groups = snv.split_by_category(gene_sets)
    assert set(groups) == {"Hallmark", "KEGG", "Reactome", "GO biological process"}
    assert sum(len(g) for g in groups.values()) == len(gene_sets)
    assert list(groups)[0] == "Hallmark"  # largest first
    custom = snv.split_by_category(gene_sets, {p: "A" for p in list(gene_sets)[:10]}, min_sets=10)
    assert set(custom) == {"A"} and len(custom["A"]) == 10  # the other categories fall below min_sets


def test_msigdb_collection_prefixes():
    assert pathway_category("GOCC_NUCLEUS") == "GO cellular component"
    assert pathway_category("KEGG_MEDICUS_REFERENCE_X") == "KEGG MEDICUS"
    assert pathway_category("HP_SEIZURE") == "Human phenotype"
    assert pathway_category("chr1p36") == "Positional"
    assert pathway_category("GSE1234_A_VS_B_UP") == "ImmuneSigDB"
    assert pathway_category("FOXO3_TARGET_GENES") == "TF/miRNA targets"


def test_category_maps_build_save_plot(gene_sets, gene_values, tmp_path):
    with pytest.warns(UserWarning, match="skipped"):  # GO BP has only 4 demo sets
        maps = snv.build_category_maps(gene_sets, methods=("tsne",), min_sets=5, min_size=5, n_neighbors=4,
                                       min_cluster_label=2)
    assert set(maps) == {"Hallmark", "KEGG", "Reactome"}
    for cat, per in maps.items():
        assert set(per) == {"tsne"}
        assert per["tsne"].attrs["category"] == cat and per["tsne"].attrs["method"] == "tsne"
        assert set(per["tsne"]["pathway"]) == {p for p in gene_sets if pathway_category(p) == cat}
    index = snv.save_category_maps(maps, tmp_path / "maps")
    assert len(index) == 3 and (tmp_path / "maps" / "category_maps.csv").exists()
    assert (tmp_path / "maps" / "hallmark_map_tsne.csv").exists() and (tmp_path / "maps" / "hallmark_map_tsne.pdf").stat().st_size > 0
    back = snv.load_map(tmp_path / "maps" / "kegg_map_tsne.csv")
    assert len(back) == 8 and back.attrs["method"] == "tsne"
    table = snv.enrich(list(gene_values), gene_sets)
    fig = snv.plot_category_maps(maps, tmp_path / "grid.pdf", table=table)
    assert len(fig.axes) >= 3 and (tmp_path / "grid.pdf").stat().st_size > 0
    import matplotlib.pyplot as plt
    plt.close(fig)
    # a map of one category can be reused by the pipeline like any other map
    res = snv.run_pipeline(genes=gene_values, gene_sets=gene_sets, outdir=str(tmp_path / "run"), similarity_cutoff=0.3,
                           pathway_map=str(tmp_path / "maps" / "hallmark_map_tsne.csv"), formats=("csv", "pdf"))
    hits = res["pathway_map_hits"]
    assert len(hits) and hits["pathway"].str.startswith("HALLMARK_").all()


def test_map_cli_by_category(tmp_path, capsys):
    out = tmp_path / "cat"
    main(["map", "--library", "demo", "--by-category", "--method", "tsne", "--min-size", "5", "--neighbors", "4",
          "--min-sets", "5", "-o", str(out)])
    text = capsys.readouterr().out
    assert "4 categories" in text and "3 maps written" in text
    assert (out / "category_maps.csv").exists() and (out / "category_maps.pdf").exists()
    index = pd.read_csv(out / "category_maps.csv")
    assert set(index["category"]) == {"Hallmark", "KEGG", "Reactome"} and (index["method"] == "tsne").all()


def test_map_cli_one_category_per_gmt(gene_sets, tmp_path, capsys):
    a = tmp_path / "setA.v2025.1.Mm.symbols.gmt"
    b = tmp_path / "setB.gmt"
    names = list(gene_sets)
    snv.write_gmt({n: gene_sets[n] for n in names[:16]}, a)
    snv.write_gmt({n: gene_sets[n] for n in names[16:]}, b)
    out = tmp_path / "cat"
    main(["map", "--gmt", str(a), "--gmt", str(b), "--by-category", "--method", "tsne", "--min-size", "5",
          "--neighbors", "4", "--min-sets", "5", "-o", str(out)])
    text = capsys.readouterr().out
    assert "2 categories" in text
    index = pd.read_csv(out / "category_maps.csv")
    assert set(index["category"]) == {"setA", "setB"} and set(index["n_sets"]) == {16}


def test_plot_map_color_by_category(pmap, gene_sets, gene_values, tmp_path):
    import matplotlib.pyplot as plt
    fig = snv.plot_pathway_map(pmap, None, tmp_path / "cat.pdf", color_by_category=True)
    assert [t.get_text().split(" (")[0] for t in fig.legends[0].texts] == ["Hallmark", "KEGG", "Reactome", "GO biological process"]
    plt.close(fig)
    cmp_table = snv.compare_conditions({"a": list(gene_values)[:40], "b": list(gene_values)[20:]}, gene_sets)
    fig = snv.plot_pathway_map(pmap, cmp_table, None, color_by_category=True, categories={p: "X" for p in gene_sets})
    assert fig.legends[0].get_title().get_text() == "Collection" and fig.axes[0].get_legend().get_title().get_text() == "Enriched in"
    assert (tmp_path / "cat.pdf").stat().st_size > 0
    plt.close(fig)


# ---------------------------------------------------------------- species / MSigDB

def test_species_helpers():
    from simpnetworkviz.genesets import msigdb_file
    assert snv.normalize_species("Mouse") == "mouse" and snv.normalize_species(9606) == "human" and snv.normalize_species("Mm") == "mouse"
    assert snv.species_taxon("mouse") == 10090 and snv.species_taxon(9606) == 9606
    assert msigdb_file("GO_BP", "mouse") == "m5.go.bp.v2025.1.Mm.symbols"
    assert msigdb_file("hallmark", "human") == "h.all.v2025.1.Hs.symbols"
    with pytest.raises(ValueError, match="no KEGG collection for mouse"):
        msigdb_file("KEGG", "mouse")
    with pytest.raises(ValueError, match="Unknown MSigDB collection"):
        msigdb_file("Nope", "human")
    with pytest.raises(ValueError, match="Unknown species"):
        snv.normalize_species("rat")


def test_msigdb_source_uses_fetcher(monkeypatch, gene_sets):
    calls = []

    def fake(collection, species="human", **kw):
        calls.append((collection, species))
        return dict(gene_sets)

    monkeypatch.setattr(snv.genesets, "fetch_msigdb", fake)
    gs = snv.load_gene_sets("msigdb:Hallmark,msigdb:GO_BP", species="mouse")
    assert calls == [("Hallmark", "mouse"), ("GO_BP", "mouse")]
    assert len(gs) == 2 * len(gene_sets)  # identical names from two sources are prefixed, not merged


def test_cli_species_and_msigdb(monkeypatch, gene_sets, gene_values, tmp_path, capsys):
    monkeypatch.setattr(snv.genesets, "fetch_msigdb", lambda collection, species="human", **kw: dict(gene_sets))
    genes = tmp_path / "genes.txt"
    genes.write_text("\n".join(gene_values))
    main(["run", str(genes), "--msigdb", "Hallmark", "--species", "mouse", "-o", str(tmp_path / "out"), "--formats", "csv",
          "--gene-network", "none", "--cutoff", "0.3"])
    params = json.load(open(tmp_path / "out" / "run_params.json"))
    assert params["species"] == "mouse" and params["gene_sets"] == "msigdb:Hallmark"
    with pytest.raises(SystemExit):
        main(["run", str(genes), "--msigdb", "Hallmark", "--species", "rat"])


# ---------------------------------------------------------------- interactive map

def test_write_map_html(pmap, gene_sets, gene_values, tmp_path):
    import json, re
    cmp_table = snv.compare_conditions({"a": list(gene_values)[:40], "b": list(gene_values)[20:]}, gene_sets)
    out = snv.write_map_html(pmap, tmp_path / "map.html", cmp_table, title="Test map")
    text = open(out, encoding="utf-8").read()
    assert "<title>Test map</title>" in text and "plotly" in text.lower()
    data = json.loads(re.search(r"const D = (\{.*?\});\n", text, re.S).group(1))
    assert data["n"] == len(pmap) and len(data["library"]["name"]) == len(pmap)
    assert {h["condition"] for h in data["hits"]} <= {"a", "b"} and data["hits"]
    assert all(len(r["hull"]) >= 4 and r["label"] and r["color"].startswith("#") for r in data["regions"])
    assert all(r["hull"][0] == r["hull"][-1] for r in data["regions"])  # closed, smoothed rings
    h = data["hits"][0]
    assert len(h["allGenes"]) == len(h["name"]) == len(h["pvalue"]) == len(h["setSize"]) and all(v is None or v > 0 for v in h["pvalue"])
    assert "selGroup" in text and "pathwayList" in text and "Significant in this data" in text
    assert 'id="regRows"' in text and 'id="condRows"' in text and "plotly_click" in text  # sidebar panels, click-to-detail
    assert len(data["regions"]) == pmap.loc[pmap["cluster"] > 0, "cluster"].nunique()
    # bare library, and the pipeline writes it next to the PDF
    snv.write_map_html(pmap, tmp_path / "lib.html")
    assert (tmp_path / "lib.html").stat().st_size > 1000
    snv.save_map(pmap, tmp_path / "m.csv")
    res = snv.run_pipeline(genes=gene_values, gene_sets=gene_sets, outdir=str(tmp_path / "run"), similarity_cutoff=0.3,
                           pathway_map=str(tmp_path / "m.csv"), formats=("csv", "html"), gene_network_mode=None)
    assert os.path.exists(res["files"]["pathway_map"]["html"])
