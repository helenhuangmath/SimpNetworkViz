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

st.set_page_config(page_title="SimpNetworkViz", layout="wide")
st.title("SimpNetworkViz")
st.caption("Gene list → pathway enrichment → clustered pathway network & gene network")

# ---------------------------------------------------------------- inputs
with st.sidebar:
    st.header("1. Input")
    mode = st.radio("Start from", ["Gene list", "Enrichment table", "Pathway list"])
    gene_text = pathway_text = enrich_file = None
    if mode == "Gene list":
        use_demo = st.checkbox("Use demo genes", value=False)
        default = open(demo_path("demo_genes.tsv")).read() if use_demo else ""
        gene_text = st.text_area("Genes (one per line; optional 2nd column = log2FC, tab or comma separated)",
                                 value=default, height=220)
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
    metric = st.selectbox("Similarity", ["combined", "jaccard", "overlap"])
    cutoff = st.slider("Similarity cutoff", 0.1, 0.9, 0.375 if lib_choice != "demo" else 0.3, 0.025)
    cluster = st.selectbox("Clustering", ["mcl", "louvain", "components"])
    color_by = st.selectbox("Pathway colors", ["auto", "cluster", "direction", "significance", "nes"])
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


if go:
    try:
        gs = load_sets(lib_choice if lib_choice not in ("Upload GMT", "Other Enrichr library") else None,
                       other_lib, gmt_upload.getvalue() if gmt_upload else None)
        genes = parse_genes(gene_text)
        background = list(parse_genes(bg_text) or {}) or None
        enrichment = pathways = None
        if mode == "Enrichment table":
            if enrich_file is None:
                st.error("Upload an enrichment table."); st.stop()
            suffix = os.path.splitext(enrich_file.name)[1] or ".tsv"
            with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as fh:
                fh.write(enrich_file.getvalue())
            enrichment = snv.read_enrichment(fh.name)
        elif mode == "Pathway list":
            pathways = [l.strip() for l in (pathway_text or "").splitlines() if l.strip()]
            if not pathways:
                st.error("Enter at least one pathway."); st.stop()
        elif not genes:
            st.error("Enter a gene list."); st.stop()
        outdir = tempfile.mkdtemp(prefix="snv_")
        with st.spinner("Running enrichment and building networks…"):
            res = snv.run_pipeline(
                genes=genes, enrichment=enrichment, pathways=pathways, gene_sets=gs, outdir=outdir,
                background=background, fdr_cutoff=fdr, max_pathways=max_pw, similarity_cutoff=cutoff,
                similarity_metric=metric, cluster_method=cluster, gene_network_mode=gmode, color_by=color_by,
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
