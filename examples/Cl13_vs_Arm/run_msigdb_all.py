"""Cl13 (TEX) vs Arm (TMEM) gene groups against every MSigDB collection, on the pre-built maps.

Default (--species mouse): the MSigDB *mouse* collections, downloaded once by the package
(Hallmark, GO BP/CC/MF, tumor phenotype, Reactome, WikiPathways, BioCarta, GTRD TF targets, cell
types) with their maps in ../msigdb_maps/mouse, plus the overall mouse GO + Reactome + WikiPathways
map and, as the one cross-species exception, human C7 ImmuneSigDB (no mouse version exists;
mouse symbols are matched to it case-insensitively). --species human uses the v7.4 human files
given with --msigdb-dir and the maps in ../msigdb_maps/{h.all,c2.cp,c5.all,c3.tft,c6.all,c7.all,overall}.

Per collection this writes, under figures/msigdb_<species>/:

  <collection>_network.pdf                    comparison network (pie nodes), 80 most significant sets
  <collection>_map_{tsne,umap}.pdf            pooled Cl13 vs Arm hits on the collection's map
  <collection>_within_Cl13_map_tsne.pdf       the six Cl13 subgroups on the t-SNE map
  <collection>_map_{tsne,umap}_regions.csv    hits per region and condition
  summary.csv                                 one row per collection
  ../../docs/images/msigdb_all_<species>_{tsne,umap}.png  one panel per collection

    python examples/Cl13_vs_Arm/run_msigdb_all.py --species mouse --map-fdr "ImmuneSigDB=1e-8"
    python examples/Cl13_vs_Arm/run_msigdb_all.py --species human --msigdb-dir /path/to/v7.4_gmts
"""
import argparse
import glob
import os
import warnings

import matplotlib.pyplot as plt
import pandas as pd

import simpnetworkviz as snv
from simpnetworkviz.pathmap import category_slug

HERE = os.path.dirname(os.path.abspath(__file__))
MAPS = os.path.join(HERE, "..", "msigdb_maps")
IMG = os.path.join(HERE, "..", "..", "docs", "images")

ap = argparse.ArgumentParser()
ap.add_argument("--species", default="mouse", help="mouse (MSigDB mouse collections, downloaded once) or human (--msigdb-dir)")
ap.add_argument("--msigdb-dir", default=None, help="human only: folder with the v7.4 *.symbols.gmt files")
ap.add_argument("--out", default=None, help="default figures/msigdb_<species>")
ap.add_argument("--fdr", type=float, default=0.05, help="FDR cutoff for the map overlays (networks use 0.05)")
ap.add_argument("--map-fdr", default=None, help="per-collection override, e.g. 'ImmuneSigDB=1e-8,c7.all=1e-4'")
args = ap.parse_args()
warnings.simplefilter("ignore", UserWarning)
SPECIES = snv.normalize_species(args.species)
args.out = args.out or os.path.join(HERE, "figures", f"msigdb_{SPECIES}")
os.makedirs(args.out, exist_ok=True)
map_fdr = {k: float(v) for k, v in (kv.split("=") for kv in (args.map_fdr or "").split(",") if kv)}

cl13 = snv.read_gene_groups(os.path.join(HERE, "data", "Cl13_groups.csv"))
arm = snv.read_gene_groups(os.path.join(HERE, "data", "Arm_groups.csv"))
pooled = {"Cl13 (TEX)": {g: None for grp in cl13.values() for g in grp},
          "Arm (TMEM)": {g: None for grp in arm.values() for g in grp}}

# collections: {name: (gene sets, map file stem)}, split like `snv map --by-category`
collections = {}


def add_collection(gs, folder, file_label):
    for cat, sub in snv.split_by_category(gs).items():
        cat = file_label if cat == "Other" else cat
        stem = os.path.join(MAPS, folder, f"{category_slug(cat)}_map")
        if os.path.exists(stem + "_tsne.csv"):
            collections[cat] = (sub, stem)
        else:
            print(f"  no map for {cat} in {folder}, skipped")


if SPECIES == "mouse":
    for c in ["Hallmark", "GO", "CP", "CellType", "TFT"]:
        add_collection(snv.fetch_msigdb(c, species="mouse"), os.path.join("mouse", c), c)
    overall = {k: v for c in ("GO", "CP") for k, v in snv.fetch_msigdb(c, species="mouse").items()
               if snv.pathway_category(k) in {"GO biological process", "GO cellular component", "GO molecular function", "Reactome", "WikiPathways"}}
    if os.path.exists(os.path.join(MAPS, "mouse", "overall", "overall_map_tsne.csv")):
        collections["GO + Reactome + WP (mouse)"] = (overall, os.path.join(MAPS, "mouse", "overall", "overall_map"))
    if os.path.exists(os.path.join(MAPS, "c7.all", "immunesigdb_map_tsne.csv")):  # human only; no mouse C7 exists
        c7 = snv.fetch_msigdb("ImmuneSigDB", species="human")
        collections["ImmuneSigDB (human C7)"] = (c7, os.path.join(MAPS, "c7.all", "immunesigdb_map"))
else:
    if not args.msigdb_dir:
        raise SystemExit("--species human needs --msigdb-dir with the v7.4 *.symbols.gmt files")
    for f in ["h.all", "c5.all", "c2.cp", "c3.tft", "c6.all", "c7.all"]:
        path = os.path.join(args.msigdb_dir, f"{f}.v7.4.symbols.gmt")
        if os.path.exists(path):
            add_collection(snv.load_gene_sets(path), f, f)
        else:
            print(f"missing {path}, skipped")
    keep = {"GO biological process", "GO cellular component", "GO molecular function", "KEGG", "Reactome", "WikiPathways"}
    merged = {k: v for c in keep if c in collections for k, v in collections[c][0].items()}
    if merged and os.path.exists(os.path.join(MAPS, "overall", "overall_map_tsne.csv")):
        collections["GO + KEGG + Reactome + WP"] = (merged, os.path.join(MAPS, "overall", "overall_map"))
print(f"{len(collections)} collections with maps:", ", ".join(collections))

rows, panels = [], {}
for cat, (gs, stem) in collections.items():
    slug = category_slug(cat)
    table = snv.compare_conditions(pooled, gs)
    within = snv.compare_conditions(cl13, gs)
    sig = table[table["fdr"] <= 0.05]
    row = {"collection": cat, "sets": len(gs), "hits_Cl13": int((sig["condition"] == "Cl13 (TEX)").sum()),
           "hits_Arm": int((sig["condition"] == "Arm (TMEM)").sum())}
    try:
        net = snv.build_pathway_network(table, gs, max_pathways=80, similarity_cutoff=0.375)
        fig = net.plot(None, title=f"Cl13 vs Arm, {cat}")
        fig.savefig(os.path.join(args.out, f"{slug}_network.pdf"), facecolor="white", bbox_inches="tight"); plt.close(fig)
        row.update(network_sets=net.graph.number_of_nodes(), edges=net.graph.number_of_edges(), clusters=len(net.group_labels),
                   top_clusters="; ".join(net.group_labels[g] for g in list(net.group_labels)[:3]))
    except ValueError as exc:
        row.update(network_sets=0, edges=0, clusters=0, top_clusters=str(exc)[:60])
    fdr = map_fdr.get(cat, args.fdr)
    row["map_fdr"] = fdr
    for method in ("tsne", "umap"):
        pmap = snv.load_map(f"{stem}_{method}.csv")
        panels.setdefault(cat, {})[method] = pmap
        hits = snv.map_hits(pmap, table, fdr)
        reg = snv.region_summary(pmap, hits)
        reg.to_csv(os.path.join(args.out, f"{slug}_map_{method}_regions.csv"), index=False)
        multi_db = cat.startswith("GO +")
        fig = snv.plot_pathway_map(pmap, table, os.path.join(args.out, f"{slug}_map_{method}.pdf"), fdr_cutoff=fdr, label_top=0,
                                   title=f"Cl13 vs Arm on {cat} ({method.upper()}), FDR ≤ {fdr:g}", color_by_category=multi_db)
        plt.close(fig)
        if method == "tsne":
            fig = snv.plot_pathway_map(pmap, within, os.path.join(args.out, f"{slug}_within_Cl13_map_tsne.pdf"), fdr_cutoff=fdr,
                                       label_top=0, title=f"Six Cl13 subgroups on {cat} (TSNE), FDR ≤ {fdr:g}")
            plt.close(fig)
            row["on_map"] = len(pmap)
            row["hits_on_map"] = int(hits["pathway"].nunique())
            for cond, key in [("Cl13 (TEX)", "top_region_Cl13"), ("Arm (TMEM)", "top_region_Arm")]:
                col = f"n_{cond}"
                best = reg.sort_values(col, ascending=False).iloc[0] if col in reg and reg[col].max() > 0 else None
                row[key] = f"{best['cluster_label']} ({int(best[col])}/{int(best['n_pathways_in_region'])})" if best is not None else ""
    rows.append(row)
    print(f"  {cat:<28} {row['sets']:>5} sets  hits Cl13/Arm {row['hits_Cl13']}/{row['hits_Arm']}  "
          f"network {row['network_sets']} sets {row['edges']} edges {row['clusters']} clusters  top Cl13 region: {row.get('top_region_Cl13', '')}")

summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(args.out, "summary.csv"), index=False)

# overview: one panel per collection, pooled hits, t-SNE and UMAP as separate figures
all_tables = pd.concat([snv.compare_conditions(pooled, gs) for gs, _ in collections.values()], ignore_index=True)
for method in ("tsne", "umap"):
    cats = list(panels)
    ncol = 4
    nrow = -(-len(cats) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 4.4 * nrow), squeeze=False, facecolor="white")
    for i, cat in enumerate(cats):
        ax = axes[i // ncol][i % ncol]
        fdr = map_fdr.get(cat, args.fdr)
        snv.plot_pathway_map(panels[cat][method], all_tables, None, fdr_cutoff=fdr, label_top=0, title=cat, width=4.6, ax=ax,
                             show_region_labels=False)
        if ax.get_legend() and i != ncol - 1:
            ax.get_legend().remove()
    for j in range(len(cats), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle(f"Cl13 (TEX) vs Arm (TMEM) on the MSigDB {SPECIES} collections, {method.upper()} maps", fontsize=14, y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, f"all_collections_{method}.pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(os.path.join(IMG, f"msigdb_all_{SPECIES}_{method}.png"), dpi=100, facecolor="white", bbox_inches="tight")
    plt.close(fig)
print("\nSummary ->", os.path.join(args.out, "summary.csv"))
print(summary[["collection", "sets", "on_map", "hits_Cl13", "hits_Arm", "network_sets", "edges", "clusters", "top_region_Cl13", "top_region_Arm"]].to_string(index=False))
