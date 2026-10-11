"""Cl13 (exhausted, TEX) vs Arm (memory, TMEM) CD8 T-cell gene groups -> pathway networks.

Inputs (data/): one CSV per condition, one column per gene subgroup (mouse symbols).
Gene sets (genesets/): MSigDB mouse Hallmark (bundled). GO Biological Process is
optional: download m5.go.bp.v2025.1.Mm.symbols.gmt from
https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp and pass --gobp PATH.

Outputs (results/), all figures as PDF:
  1_Cl13_vs_Arm_pooled/<library>/   two conditions, all genes of each condition pooled
  2_matched_groups/G<n>/            TEX_G<n> vs TMEM_G<n>
  3_within_Cl13/, 3_within_Arm/     the six subgroups of one condition as pie slices
  4_pathway_map/                    the GO BP library map with the hits of 1 and 3 overlaid
                                    (t-SNE by default; --map-method umap for the UMAP version)
  summary.csv                       one row per analysis

    python examples/Cl13_vs_Arm/run_comparison.py [--gobp path/to/m5.go.bp.v2025.1.Mm.symbols.gmt] [--map-method umap]
"""
import argparse
import os
import warnings

import pandas as pd

import simpnetworkviz as snv

HERE = os.path.dirname(os.path.abspath(__file__))
FORMATS = ("pdf", "html", "csv")

ap = argparse.ArgumentParser()
ap.add_argument("--gobp", default=None, help="optional GO Biological Process mouse GMT file")
ap.add_argument("--out", default=os.path.join(HERE, "results"))
ap.add_argument("--map-method", choices=["tsne", "umap"], default="tsne",
                help="embedding for the GO BP pathway map (section 4); the bundled map is t-SNE, "
                     "UMAP needs umap-learn and writes gobp_map_umap.csv / map_*_umap.pdf")
args = ap.parse_args()
warnings.simplefilter("ignore", UserWarning)

cl13 = snv.read_gene_groups(os.path.join(HERE, "data", "Cl13_groups.csv"))   # {TEX_G1: {gene: None}, ...}
arm = snv.read_gene_groups(os.path.join(HERE, "data", "Arm_groups.csv"))     # {TMEM_G1: ..., ...}
libraries = {"Hallmark": snv.load_gene_sets(os.path.join(HERE, "genesets", "mh.all.v2025.1.Mm.symbols.gmt"))}
if args.gobp and os.path.exists(args.gobp):
    libraries["GO_BP"] = snv.load_gene_sets(args.gobp, min_size=10, max_size=500)
else:
    msg = f"GO BP file not found ({args.gobp}); " if args.gobp else "No GO BP file given; "
    print(msg + "running Hallmark only.")

rows = []


def run(name, conditions, library, outdir, **kw):
    gs = libraries[library]
    try:
        res = snv.run_pipeline(conditions=conditions, gene_sets=gs, outdir=outdir, formats=FORMATS,
                               gene_network_mode="shared", **kw)
    except ValueError as exc:  # nothing significant
        print(f"  {name} [{library}]: {exc}")
        rows.append({"analysis": name, "library": library, "pathways": 0, "clusters": 0, "note": str(exc)})
        return None
    net = res["pathway_network"]
    n_hl = sum(1 for _, d in net.graph.nodes(data=True) if d.get("highlight", True))
    rows.append({"analysis": name, "library": library, "pathways": n_hl, "clusters": len(net.group_labels),
                 "conditions": ";".join(net.conditions), "note": ""})
    print(f"  {name} [{library}]: {n_hl} enriched pathways, {len(net.group_labels)} clusters")
    return res


# 1. Pooled: every Cl13 gene vs every Arm gene -------------------------------------
pooled = {"Cl13 (TEX)": {g: None for grp in cl13.values() for g in grp},
          "Arm (TMEM)": {g: None for grp in arm.values() for g in grp}}
print("1. Cl13 vs Arm, pooled genes")
for lib in libraries:
    run("Cl13_vs_Arm_pooled", pooled, lib, os.path.join(args.out, "1_Cl13_vs_Arm_pooled", lib),
        library_network=(lib == "Hallmark"), max_pathways=80)

# 2. Matched subgroups: TEX_Gn vs TMEM_Gn ---------------------------------------------
print("2. Matched subgroups")
for n in range(1, 7):
    pair = {f"TEX_G{n}": cl13[f"TEX_G{n}"], f"TMEM_G{n}": arm[f"TMEM_G{n}"]}
    for lib in libraries:
        run(f"G{n}_TEX_vs_TMEM", pair, lib, os.path.join(args.out, "2_matched_groups", f"G{n}", lib),
            library_network=(lib == "Hallmark"), max_pathways=60)

# 3. Within each condition: six subgroups as pie slices -------------------------------
print("3. Subgroups within each condition")
for label, groups in [("Cl13", cl13), ("Arm", arm)]:
    for lib in libraries:
        run(f"within_{label}", groups, lib, os.path.join(args.out, f"3_within_{label}", lib),
            library_network=(lib == "Hallmark"), max_pathways=80)

# 4. The GO BP library map: one fixed coordinate system for every analysis -------------
if "GO_BP" in libraries:
    import matplotlib.pyplot as plt

    suffix = "" if args.map_method == "tsne" else f"_{args.map_method}"
    map_csv = os.path.join(HERE, "genesets", f"gobp_map{suffix}.csv")
    if os.path.exists(map_csv):
        pmap = snv.load_map(map_csv)                     # built once with: snv map --gmt <GO BP gmt> --method tsne|umap -o gobp_map.csv
    else:
        pmap = snv.build_pathway_map(libraries["GO_BP"], method=args.map_method, seed=0)
        snv.save_map(pmap, map_csv)
    mdir = os.path.join(args.out, "4_pathway_map")
    os.makedirs(mdir, exist_ok=True)
    print("4. Library map:", len(pmap), "GO BP terms,", int(pmap["cluster"].max()), "regions")
    plt.close(snv.plot_pathway_map(pmap, None, os.path.join(mdir, f"map_GOBP{suffix}.pdf")))
    for name, conditions in [("pooled_Cl13_vs_Arm", pooled), ("within_Cl13", cl13), ("within_Arm", arm)]:
        table = snv.compare_conditions(conditions, libraries["GO_BP"])
        hits = snv.map_hits(pmap, table)
        snv.region_summary(pmap, hits).to_csv(os.path.join(mdir, f"map_{name}{suffix}_regions.csv"), index=False)
        plt.close(snv.plot_pathway_map(pmap, table, os.path.join(mdir, f"map_{name}{suffix}.pdf"), label_top=8))
        print(f"  {name}: {hits['pathway'].nunique()} enriched terms placed on the map")

pd.DataFrame(rows).to_csv(os.path.join(args.out, "summary.csv"), index=False)
print("\nSummary ->", os.path.join(args.out, "summary.csv"))
