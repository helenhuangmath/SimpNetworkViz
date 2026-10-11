"""Point-and-click web app. Launch with ``snv gui`` (needs ``pip install simpnetworkviz[gui]``)."""
import io
import os
import tempfile
import zipfile

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

import simpnetworkviz as snv
from simpnetworkviz.genesets import demo_path, gene_values_from_frame, normalize_gene, parse_gmt_lines
from simpnetworkviz.network import read_theme_file

st.set_page_config(page_title="SimpNetworkViz", layout="wide")
st.title("SimpNetworkViz")
st.caption("Gene list → pathway enrichment → clustered pathway network & gene network")

# ---------------------------------------------------------------- inputs
with st.sidebar:
    st.header("1. Input")
    mode = st.radio("Start from", ["Gene list", "Compare conditions", "Single-cell markers", "Enrichment table", "Pathway list"])
    gene_text = pathway_text = enrich_file = markers_file = None
    condition_texts = {}
    method = "ora"
    if mode == "Gene list":
        use_demo = st.checkbox("Use demo genes", value=False)
        default = open(demo_path("demo_genes.tsv")).read() if use_demo else ""
        gene_text = st.text_area("Genes (one per line; optional 2nd column = log2FC, tab or comma separated)",
                                 value=default, height=220)
        method = st.radio("Test", ["ora", "gsea"], horizontal=True,
                          help="ORA: hypergeometric test on your gene list. GSEA: preranked GSEA; paste ALL measured genes with a statistic.")
    elif mode == "Compare conditions":
        n_cond = st.number_input("Number of conditions", 2, 8, 2)
        for i in range(int(n_cond)):
            c1, c2 = st.columns([1, 2])
            name = c1.text_input("Name", f"Condition {i + 1}", key=f"cname{i}")
            txt = c2.text_area("Genes (+ optional log2FC)", height=90, key=f"ctext{i}")
            condition_texts[name] = txt
        method = st.radio("Test", ["ora", "gsea"], horizontal=True)
    elif mode == "Single-cell markers":
        markers_file = st.file_uploader("Marker table (Seurat FindAllMarkers / scanpy rank_genes_groups; CSV/TSV/XLSX)")
        sc_padj = st.number_input("Marker adj. p-value cutoff", 0.0, 1.0, 0.05, format="%.3f")
        sc_lfc = st.number_input("Min |log2FC|", 0.0, 5.0, 0.25)
        sc_top = st.number_input("Max markers per cluster", 10, 2000, 200)
        sc_dir = st.selectbox("Marker sign", ["up", "down", "both"])
        sc_clusters = st.text_input("Clusters to include (comma separated, empty = all)")
    elif mode == "Enrichment table":
        enrich_file = st.file_uploader("Enrichment result (CSV/TSV/XLSX from Enrichr, clusterProfiler, g:Profiler, GSEA…)")
        gene_text = st.text_area("Optional gene list with log2FC (colors gene nodes)", height=100)
    else:
        pathway_text = st.text_area("Pathway names (one per line, must match the gene-set library)", height=200)
        gene_text = st.text_area("Optional genes of interest (restricts gene network)", height=100)

    st.header("2. Gene sets")
    lib_choice = st.selectbox("Library", ["demo"] + list(snv.POPULAR_LIBRARIES) + ["Upload GMT", "Other Enrichr library"])
    gmt_upload = other_lib = None
    if lib_choice == "Upload GMT":
        gmt_upload = st.file_uploader("GMT file", type=["gmt", "txt"])
    elif lib_choice == "Other Enrichr library":
        other_lib = st.text_input("Enrichr library name", "WikiPathways_2023_Human")
    bg_text = st.text_area("Background genes (optional, recommended)", height=80)

    st.header("3. Network")
    fdr = st.slider("FDR cutoff", 0.001, 0.25, 0.05)
    max_pw = st.slider("Max pathways", 10, 300, 100)
    whole_library = st.checkbox("Whole library (highlight enriched pathways)", value=False,
                                help="Draw every pathway of the library; enriched ones are colored, the rest grey.")
    metric = st.selectbox("Similarity", ["combined", "jaccard", "overlap"])
    sim_genes = st.selectbox("Similarity on", ["full gene sets", "hit genes only"])
    cutoff = st.slider("Similarity cutoff", 0.1, 0.9, 0.375 if lib_choice != "demo" else 0.3, 0.025)
    cluster = st.selectbox("Grouping", ["mcl", "louvain", "greedy", "label_propagation", "hierarchical", "components",
                                        "category", "keyword", "none"],
                           help="category = by source (Hallmark/KEGG/...); keyword = by your themes below")
    themes_text = ""
    if cluster == "keyword":
        themes_text = st.text_area("Themes: one per line as  Theme<TAB>keyword1, keyword2",
                                   "Proliferation\tcycle, mitotic, E2F, replication\nImmune\tinterferon, cytokine, antigen")
    color_by = st.selectbox("Pathway colors", ["auto", "cluster", "direction", "significance", "nes", "condition"])
    gmode = st.selectbox("Gene network", ["shared", "cnet", "string"],
                         help="shared = genes linked by shared pathways; cnet = gene–pathway; string = STRING PPI (internet)")
    max_genes = st.slider("Max genes", 20, 400, 150)
    go = st.button("Build networks", type="primary", use_container_width=True)


def parse_genes(text):
    if not text or not text.strip():
        return None
    lines = [l for l in text.strip().splitlines() if l.strip()]
    sep = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else None)
    if sep is None:
        return {normalize_gene(t): None for l in lines for t in l.split()}
    df = pd.read_csv(io.StringIO("\n".join(lines)), sep=sep, header=None, dtype=str)
    try:
        float(df.iloc[0, 1])
    except (ValueError, IndexError):
        df.columns = df.iloc[0]
        df = df.iloc[1:]
    for c in df.columns[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return gene_values_from_frame(df)


@st.cache_data(show_spinner="Loading gene sets…")
def load_sets(choice, other, gmt_bytes):
    if gmt_bytes is not None:
        return parse_gmt_lines(gmt_bytes.decode("utf-8", "replace").splitlines())
    return snv.load_gene_sets(other or choice)


def _tmp_upload(upload):
    suffix = os.path.splitext(upload.name)[1] or ".tsv"
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as fh:
        fh.write(upload.getvalue())
    return fh.name


if go:
    try:
        gs = load_sets(lib_choice if lib_choice not in ("Upload GMT", "Other Enrichr library") else None,
                       other_lib, gmt_upload.getvalue() if gmt_upload else None)
        genes = parse_genes(gene_text)
        background = list(parse_genes(bg_text) or {}) or None
        enrichment = pathways = conditions = None
        if mode == "Compare conditions":
            conditions = {n: parse_genes(t) for n, t in condition_texts.items()}
            conditions = {n: g for n, g in conditions.items() if g}
            if len(conditions) < 2:
                st.error("Enter genes for at least two conditions."); st.stop()
        elif mode == "Single-cell markers":
            if markers_file is None:
                st.error("Upload a marker table."); st.stop()
            conditions = snv.markers_to_conditions(
                _tmp_upload(markers_file), padj_cutoff=sc_padj, min_logfc=sc_lfc, top_n=int(sc_top), direction=sc_dir,
                clusters=[c.strip() for c in sc_clusters.split(",") if c.strip()] or None)
            st.info(", ".join(f"{k}: {len(v)} markers" for k, v in conditions.items()))
        elif mode == "Enrichment table":
            if enrich_file is None:
                st.error("Upload an enrichment table."); st.stop()
            enrichment = snv.read_enrichment(_tmp_upload(enrich_file))
        elif mode == "Pathway list":
            pathways = [l.strip() for l in (pathway_text or "").splitlines() if l.strip()]
            if not pathways:
                st.error("Enter at least one pathway."); st.stop()
        elif not genes:
            st.error("Enter a gene list."); st.stop()
        cluster_kwargs = {}
        if cluster == "keyword":
            with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
                fh.write(themes_text)
            cluster_kwargs["themes"] = read_theme_file(fh.name)
        outdir = tempfile.mkdtemp(prefix="snv_")
        with st.spinner("Running enrichment and building networks…"):
            res = snv.run_pipeline(
                genes=genes, enrichment=enrichment, pathways=pathways, conditions=conditions, gene_sets=gs, outdir=outdir,
                background=background, method=method, fdr_cutoff=fdr, max_pathways=max_pw, similarity_cutoff=cutoff,
                similarity_metric=metric, use_full_gene_sets=sim_genes.startswith("full"), cluster_method=cluster,
                cluster_kwargs=cluster_kwargs, library_network=whole_library, gene_network_mode=gmode, color_by=color_by,
                formats=("png", "pdf", "svg", "html", "cyjs", "graphml", "csv"),
                gene_network_kwargs={"max_genes": max_genes},
            )
        st.session_state["res"] = res
        st.session_state["outdir"] = outdir
    except Exception as exc:  # show errors in the UI instead of a stack trace
        st.error(f"{type(exc).__name__}: {exc}")

res = st.session_state.get("res")
if not res:
    st.info("Choose your input in the sidebar and click **Build networks**. Tick *Use demo genes* to try it.")
    st.stop()

outdir = st.session_state["outdir"]
tabs = st.tabs(["Pathway network", "Gene network", "Gene–pathway network", "Tables", "Download"])
with tabs[0]:
    c1, c2 = st.columns([3, 2])
    c1.image(os.path.join(outdir, "pathway_network.png"))
    c2.dataframe(res["pathway_network"].cluster_table()[["cluster", "label", "n_nodes", "best_fdr"]], hide_index=True)
    if res["pathway_network"].conditions:
        c2.caption("Pie slices (clockwise from top): " + ", ".join(
            f"{i + 1} = {c}" for i, c in enumerate(res["pathway_network"].conditions)))
    with st.expander("Interactive view"):
        components.html(open(os.path.join(outdir, "pathway_network.html")).read(), height=750, scrolling=False)
with tabs[1]:
    if "gene_network" in res:
        st.image(os.path.join(outdir, "gene_network.png"))
        with st.expander("Interactive view"):
            components.html(open(os.path.join(outdir, "gene_network.html")).read(), height=750)
with tabs[2]:
    p = os.path.join(outdir, "gene_pathway_network.png")
    if os.path.exists(p):
        st.image(p)
with tabs[3]:
    if isinstance(res["enrichment"], pd.DataFrame):
        st.subheader("Enrichment")
        st.dataframe(res["enrichment"], hide_index=True)
        mp = os.path.join(outdir, "enrichment_matrix.csv")
        if os.path.exists(mp):
            st.subheader("Pathway × condition (signed -log10 FDR)")
            st.dataframe(pd.read_csv(mp, index_col=0))
    st.subheader("Pathway nodes")
    st.dataframe(res["pathway_network"].node_table(), hide_index=True)
with tabs[4]:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in sorted(os.listdir(outdir)):
            zf.write(os.path.join(outdir, f), f)
    st.download_button("Download all results (.zip)", buf.getvalue(), "simpnetworkviz_results.zip", "application/zip")
    for f in sorted(os.listdir(outdir)):
        st.download_button(f, open(os.path.join(outdir, f), "rb").read(), f, key=f)
plt.close("all")
