"""Export networks to Cytoscape (.cyjs + style), GraphML, CSV and interactive HTML."""
from __future__ import annotations

import html
import json
import os
from pathlib import Path
from typing import Dict, Optional, Sequence

import networkx as nx
import numpy as np
import pandas as pd

from .style import darken, hull_colors, node_colors, resolve_color_by

PX_PER_UNIT = 100.0
CYTOSCAPE_JS_CDN = [
    "https://cdn.jsdelivr.net/npm/cytoscape@3.30.2/dist/cytoscape.min.js",
    "https://unpkg.com/cytoscape@3.30.2/dist/cytoscape.min.js",
]


def _clean(v):
    if v is None:
        return None
    if isinstance(v, (np.floating, float)):
        return None if np.isnan(v) else float(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    return v


def cytoscape_elements(net, color_by: str = "auto", compound: bool = True) -> Dict[str, list]:
    """Cytoscape.js-style element dicts (positions in pixels, y pointing down)."""
    colors, _ = node_colors(net, color_by)
    mode = resolve_color_by(net, color_by)
    hcols = hull_colors(net, mode)
    nodes, edges = [], []
    if compound:
        for gid, label in net.group_labels.items():
            if not any(g == gid for g in net.groups.values()):
                continue
            fill, edge = hcols.get(gid, ("#F2F2F2", "#BBBBBB"))
            nodes.append({"data": {"id": f"cluster_{gid}", "label": label, "is_cluster": True, "cluster": gid,
                                   "fill": fill, "border": edge, "text_color": darken(edge, 0.45)},
                          "classes": "cluster"})
    for n, d in net.graph.nodes(data=True):
        data = {k: _clean(v) for k, v in d.items() if not isinstance(v, (list, set, tuple, dict))}
        data.update(id=str(n), color=colors[n], border_color=darken(colors[n], 0.35),
                    size_px=round(2 * net.radius[n] * PX_PER_UNIT, 2), cluster=net.groups.get(n, 0),
                    cluster_label=net.group_labels.get(net.groups.get(n, 0), ""))
        gid = net.groups.get(n, 0)
        if compound and gid and net.layout_style == "clustered":
            data["parent"] = f"cluster_{gid}"
        x, y = net.pos[n]
        nodes.append({"data": data, "position": {"x": round(x * PX_PER_UNIT, 2), "y": round(-y * PX_PER_UNIT, 2)},
                      "classes": d.get("type", "pathway")})
    ws = [d.get("weight", 1.0) for _, _, d in net.graph.edges(data=True)]
    wmin, wmax = (min(ws), max(ws)) if ws else (0, 1)
    for i, (u, v, d) in enumerate(net.graph.edges(data=True)):
        w = d.get("weight", 1.0)
        t = (w - wmin) / (wmax - wmin) if wmax > wmin else 0.5
        data = {k: _clean(val) for k, val in d.items()}
        color = "#9A9A9A"
        if net.kind == "gene-pathway":
            p = u if net.graph.nodes[u].get("type") == "pathway" else v
            color = colors[p]
        data.update(id=f"e{i}", source=str(u), target=str(v), width_px=round(0.8 + 4 * t if net.kind == "pathway" else 0.6 + 2 * t, 2),
                    color=color)
        edges.append({"data": data})
    return {"nodes": nodes, "edges": edges}


def write_cyjs(net, path, color_by: str = "auto", style_path: Optional[str] = None) -> str:
    """Write a Cytoscape desktop JSON (File ▸ Import ▸ Network from File) plus a style file.

    Cluster ids are stored in the ``cluster`` column, so in Cytoscape you can run
    AutoAnnotate with "Use existing column: cluster" to draw the cluster outlines.
    """
    el = cytoscape_elements(net, color_by, compound=False)
    doc = {
        "format_version": "1.0",
        "generated_by": "simpnetworkviz",
        "target_cytoscapejs_version": "~2.1",
        "data": {"name": Path(path).stem, "shared_name": Path(path).stem},
        "elements": el,
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=1)
    style_path = style_path or str(Path(path).with_suffix("")) + "_style.xml"
    write_cytoscape_style(style_path)
    return path


def write_cytoscape_style(path, name: str = "SimpNetworkViz") -> str:
    """Cytoscape visual style using passthrough mappings of exported columns."""
    xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<vizmap id="VizMap-simpnetworkviz" documentVersion="3.1">
  <visualStyle name="{name}">
    <network>
      <visualProperty name="NETWORK_BACKGROUND_PAINT" default="#FFFFFF"/>
    </network>
    <node>
      <dependency name="nodeSizeLocked" value="true"/>
      <visualProperty name="NODE_SHAPE" default="ELLIPSE"/>
      <visualProperty name="NODE_FILL_COLOR" default="#BDBDBD">
        <passthroughMapping attributeName="color" attributeType="string"/>
      </visualProperty>
      <visualProperty name="NODE_BORDER_PAINT" default="#555555">
        <passthroughMapping attributeName="border_color" attributeType="string"/>
      </visualProperty>
      <visualProperty name="NODE_BORDER_WIDTH" default="1.0"/>
      <visualProperty name="NODE_SIZE" default="30.0">
        <passthroughMapping attributeName="size_px" attributeType="float"/>
      </visualProperty>
      <visualProperty name="NODE_LABEL" default="">
        <passthroughMapping attributeName="label" attributeType="string"/>
      </visualProperty>
      <visualProperty name="NODE_LABEL_FONT_SIZE" default="9"/>
      <visualProperty name="NODE_LABEL_POSITION" default="N,S,c,0.00,3.00"/>
    </node>
    <edge>
      <dependency name="arrowColorMatchesEdge" value="true"/>
      <visualProperty name="EDGE_WIDTH" default="1.5">
        <passthroughMapping attributeName="width_px" attributeType="float"/>
      </visualProperty>
      <visualProperty name="EDGE_STROKE_UNSELECTED_PAINT" default="#9A9A9A">
        <passthroughMapping attributeName="color" attributeType="string"/>
      </visualProperty>
      <visualProperty name="EDGE_TRANSPARENCY" default="140"/>
    </edge>
  </visualStyle>
</vizmap>
"""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(xml)
    return path


def write_graphml(net, path, color_by: str = "auto") -> str:
    colors, _ = node_colors(net, color_by)
    H = nx.Graph()
    for n, d in net.graph.nodes(data=True):
        attrs = {k: v for k, v in ((k, _clean(v)) for k, v in d.items()) if v is not None and not isinstance(v, (list, set, tuple, dict))}
        x, y = net.pos[n]
        attrs.update(x=float(x * PX_PER_UNIT), y=float(-y * PX_PER_UNIT), color=colors[n],
                     cluster=int(net.groups.get(n, 0)), cluster_label=net.group_labels.get(net.groups.get(n, 0), ""))
        H.add_node(str(n), **attrs)
    for u, v, d in net.graph.edges(data=True):
        H.add_edge(str(u), str(v), **{k: _clean(val) for k, val in d.items() if _clean(val) is not None})
    nx.write_graphml(H, path)
    return path


# ---------------------------------------------------------------------------
# Interactive HTML
# ---------------------------------------------------------------------------

_HTML = r"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
:root { --bg:#ffffff; --panel:#f7f7f8; --text:#1f2328; --muted:#667085; --line:#e4e4e7; --accent:#3C5488; }
* { box-sizing:border-box; }
html,body { margin:0; height:100%; background:var(--bg); color:var(--text);
  font-family: Arial, "Liberation Sans", Helvetica, sans-serif; }
#app { display:flex; height:100vh; }
#main { flex:1; display:flex; flex-direction:column; min-width:0; }
header { display:flex; flex-wrap:wrap; gap:8px; align-items:center; padding:10px 14px; border-bottom:1px solid var(--line); }
header h1 { font-size:16px; margin:0 12px 0 0; font-weight:400; }
header input { padding:6px 10px; border:1px solid var(--line); border-radius:6px; min-width:180px; font-size:13px; }
button { padding:6px 10px; border:1px solid var(--line); background:#fff; border-radius:6px; cursor:pointer; font-size:13px; }
button:hover { border-color:var(--accent); }
#cy { flex:1; min-height:300px; }
aside { width:340px; border-left:1px solid var(--line); background:var(--panel); overflow:auto; padding:14px; font-size:13px; }
aside h2 { font-size:14px; margin:0 0 8px; }
aside .muted { color:var(--muted); }
.kv { display:grid; grid-template-columns: 110px 1fr; gap:4px 8px; margin:8px 0; }
.kv div:nth-child(odd) { color:var(--muted); }
.genes { display:flex; flex-wrap:wrap; gap:4px; margin-top:6px; }
.genes span { background:#fff; border:1px solid var(--line); border-radius:4px; padding:1px 5px; font-family:ui-monospace,monospace; font-size:12px; }
ul.members { padding-left:18px; margin:6px 0; }
ul.members li { cursor:pointer; margin:2px 0; }
ul.members li:hover { color:var(--accent); }
.legend .sw { display:inline-block; width:12px; height:12px; border-radius:50%; vertical-align:middle; margin-right:6px; border:1px solid #0002; }
.grad { height:10px; border-radius:3px; margin:4px 0 2px; }
@media (max-width: 800px) { #app { flex-direction:column; } aside { width:auto; height:40vh; border-left:none; border-top:1px solid var(--line); } }
</style>
__SCRIPT__
</head><body>
<div id="app">
  <div id="main">
    <header>
      <h1>__TITLE__</h1>
      <input id="search" placeholder="Search pathway or gene…" autocomplete="off">
      <button id="fit">Fit</button>
      <button id="labels">Node labels</button>
      <button id="png">Download PNG</button>
    </header>
    <div id="cy"></div>
  </div>
  <aside>
    <div id="info"><h2>Details</h2><p class="muted">Click a node or cluster to see its details. Drag nodes or whole clusters to rearrange; scroll to zoom.</p></div>
    <hr style="border:none;border-top:1px solid var(--line);margin:14px 0">
    <div class="legend" id="legend"></div>
    <hr style="border:none;border-top:1px solid var(--line);margin:14px 0">
    <h2>Clusters</h2>
    <ul class="members" id="clusters"></ul>
  </aside>
</div>
<script>
const ELEMENTS = __ELEMENTS__;
const META = __META__;
function esc(s){ return String(s ?? "").replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c])); }
function fmt(v){ if (v === null || v === undefined || v === "") return "–"; if (typeof v === "number") return (Math.abs(v) < 1e-3 && v !== 0) ? v.toExponential(2) : (+v.toFixed(3)).toString(); return esc(v); }
if (typeof cytoscape === "undefined") {
  document.getElementById("cy").innerHTML = '<p style="padding:20px">Could not load cytoscape.js (no internet?). Re-export with <code>js_path=</code> pointing at a local cytoscape.min.js to make an offline file.</p>';
}
const showLabelsInit = META.node_labels;
const cy = cytoscape({
  container: document.getElementById("cy"),
  elements: ELEMENTS,
  layout: { name: "preset", fit: true, padding: 30 },
  wheelSensitivity: 0.25,
  style: [
    { selector: "node", style: {
        "background-color": "data(color)", "border-color": "data(border_color)", "border-width": 1.2,
        "width": "data(size_px)", "height": "data(size_px)", "label": "", "font-size": 9, "color": "#222",
        "text-valign": "bottom", "text-halign": "center", "text-margin-y": 2, "text-wrap": "wrap", "text-max-width": 140,
        "text-outline-color": "#fff", "text-outline-width": 2,
        "font-family": "Arial, Liberation Sans, Helvetica, sans-serif", "font-weight": 400 } },
    { selector: "node.showlabel", style: { "label": "data(label)" } },
    { selector: "node.gene.showlabel", style: { "font-size": 8 } },
    { selector: "node.cluster", style: {
        "background-color": "data(fill)", "background-opacity": 0.85, "border-color": "data(border)", "border-width": 1.5,
        "shape": "round-rectangle", "corner-radius": 18, "padding": 14, "label": "data(label)", "text-valign": "top", "text-halign": "center",
        "font-size": 14, "font-weight": 400, "color": "data(text_color)", "text-margin-y": -4, "text-wrap": "wrap", "text-max-width": 260,
        "text-outline-width": 0 } },
    { selector: "edge", style: { "width": "data(width_px)", "line-color": "data(color)", "opacity": 0.5, "curve-style": "haystack" } },
    { selector: ".faded", style: { "opacity": 0.12 } },
    { selector: "node.hit", style: { "border-width": 4, "border-color": "#111", "label": "data(label)", "z-index": 99 } },
    { selector: "node:selected", style: { "border-width": 4, "border-color": "#111" } },
  ],
});
let labelsOn = showLabelsInit;
function applyLabels(){ cy.nodes().not(".cluster").toggleClass("showlabel", labelsOn); }
applyLabels();
document.getElementById("labels").onclick = () => { labelsOn = !labelsOn; applyLabels(); };
document.getElementById("fit").onclick = () => cy.animate({ fit: { padding: 30 } }, { duration: 300 });
document.getElementById("png").onclick = () => {
  const a = document.createElement("a"); a.href = cy.png({ full: true, scale: 3, bg: "#ffffff" });
  a.download = (META.file || "network") + ".png"; a.click();
};
function nodeInfo(d){
  const genes = (d.genes || "").split(";").filter(Boolean);
  const rows = [["Cluster", d.cluster_label || (d.cluster ? "Cluster " + d.cluster : "–")]];
  if (d.type === "pathway") {
    rows.push(["FDR", fmt(d.fdr)], ["p-value", fmt(d.pvalue)], ["Direction", fmt(d.direction)], ["NES", fmt(d.nes)],
              ["Gene-set size", fmt(d.set_size)], ["Genes in list", fmt(d.overlap)]);
  } else {
    rows.push(["log2FC", fmt(d.logfc)], ["# pathways", fmt(d.n_pathways)]);
  }
  let h = `<h2>${esc(d.label)}</h2>`;
  if (d.name && d.name !== d.label) h += `<div class="muted">${esc(d.name)}</div>`;
  h += `<div class="kv">${rows.map(r => `<div>${r[0]}</div><div>${r[1]}</div>`).join("")}</div>`;
  if (genes.length) h += `<div class="muted">Genes (${genes.length})</div><div class="genes">${genes.map(g => `<span>${esc(g)}</span>`).join("")}</div>`;
  if (d.pathways) h += `<div class="muted" style="margin-top:8px">Pathways</div><div>${esc(d.pathways)}</div>`;
  return h;
}
function clusterInfo(node){
  const kids = node.children().sort((a,b) => (a.data("fdr") ?? 1) - (b.data("fdr") ?? 1));
  return `<h2>${esc(node.data("label"))}</h2><div class="muted">${kids.length} nodes</div><ul class="members">` +
    kids.map(k => `<li data-id="${esc(k.id())}">${esc(k.data("label"))}${k.data("fdr") != null ? ` <span class="muted">(${fmt(k.data("fdr"))})</span>` : ""}</li>`).join("") + "</ul>";
}
function bindMemberClicks(root){
  root.querySelectorAll("li[data-id]").forEach(li => li.onclick = () => {
    const n = cy.getElementById(li.dataset.id); focus(n.isParent() ? n : n); });
}
const info = document.getElementById("info");
function focus(n){
  cy.elements().removeClass("faded");
  if (n.isParent()) {
    const keep = n.union(n.descendants());
    cy.elements().not(keep).not(keep.connectedEdges()).addClass("faded");
    info.innerHTML = clusterInfo(n); bindMemberClicks(info);
  } else {
    const keep = n.closedNeighborhood().union(n.ancestors());
    cy.elements().not(keep).addClass("faded");
    info.innerHTML = nodeInfo(n.data());
  }
  cy.animate({ fit: { eles: n.isParent() ? n : n.closedNeighborhood(), padding: 60 } }, { duration: 300 });
}
cy.on("tap", "node", e => focus(e.target));
cy.on("tap", e => { if (e.target === cy) cy.elements().removeClass("faded"); });
document.getElementById("search").addEventListener("input", e => {
  const q = e.target.value.trim().toUpperCase();
  cy.nodes().removeClass("hit"); cy.elements().removeClass("faded");
  if (q.length < 2) return;
  const hits = cy.nodes().not(".cluster").filter(n =>
    String(n.data("label")).toUpperCase().includes(q) || String(n.data("name") || "").toUpperCase().includes(q) ||
    (n.data("genes") || "").toUpperCase().split(";").includes(q));
  hits.addClass("hit");
  if (hits.length) { cy.elements().not(hits).not(hits.ancestors()).addClass("faded"); cy.animate({ fit: { eles: hits, padding: 80 } }, { duration: 300 }); }
});
// legend + cluster list
const lg = document.getElementById("legend");
lg.innerHTML = `<h2>Legend</h2>` + META.legend_html;
const cl = document.getElementById("clusters");
cl.innerHTML = cy.nodes(".cluster").sort((a,b) => a.data("cluster") - b.data("cluster")).map(c =>
  `<li data-id="${esc(c.id())}"><span class="sw" style="background:${c.data("fill")};border-color:${c.data("border")};border-radius:3px;display:inline-block;width:12px;height:12px;vertical-align:middle;margin-right:6px"></span>${esc(c.data("label"))} <span class="muted">(${c.children().length})</span></li>`).join("") ||
  '<li class="muted">No clusters</li>';
bindMemberClicks(cl);
</script>
</body></html>
"""


def _legend_html(net, leg) -> str:
    mode = leg["mode"]
    size_title = {"pathway": "Node size = gene-set size", "gene": "Node size = number of pathways",
                  "gene-pathway": "Pathway size = genes shown"}[net.kind]
    parts = [f"<p class='muted'>{size_title}. "]
    if net.kind == "pathway":
        parts.append("Edge width = gene-set overlap.")
    parts.append("</p>")
    if mode == "direction":
        parts.append(
            f"<div>Up-regulated</div><div class='grad' style='background:linear-gradient(90deg,#fcbba1,#a50f15)'></div>"
            f"<div>Down-regulated</div><div class='grad' style='background:linear-gradient(90deg,#c6dbef,#08306b)'></div>"
            f"<div class='muted'>Shade: -log10 FDR 0 → {leg['vmax']:.1f}</div>")
    elif mode == "significance":
        parts.append(f"<div>-log10 FDR</div><div class='grad' style='background:linear-gradient(90deg,#fed976,#bd0026)'></div>"
                     f"<div class='muted'>0 → {leg['vmax']:.1f}</div>")
    elif mode == "cluster":
        parts.append("<div class='muted'>Node color = cluster.</div>")
    else:
        v = leg["vmax"]
        parts.append(f"<div>{html.escape(leg.get('label', mode))}</div><div class='grad' style='background:linear-gradient(90deg,#2166ac,#f7f7f7,#b2182b)'></div>"
                     f"<div class='muted'>{-v:.1f} → 0 → {v:.1f}</div>")
    return "".join(parts)


def write_html(net, path, title: Optional[str] = None, color_by: str = "auto", node_labels: Optional[bool] = None,
               js_path: Optional[str] = None) -> str:
    """Write a standalone interactive network viewer (Cytoscape.js).

    ``js_path``: path to a local ``cytoscape.min.js`` to embed, making the file
    work offline. By default the library is loaded from a CDN.
    """
    el = cytoscape_elements(net, color_by)
    _, leg = node_colors(net, color_by)
    if node_labels is None:
        node_labels = bool(net.label_nodes)
    title = title or {"pathway": "Pathway network", "gene": "Gene network", "gene-pathway": "Gene–pathway network"}[net.kind]
    if js_path:
        with open(js_path, encoding="utf-8") as fh:
            script = "<script>" + fh.read().replace("</script", "<\\/script") + "</script>"
    else:
        script = (f'<script src="{CYTOSCAPE_JS_CDN[0]}"></script>\n'
                  f'<script>window.cytoscape || document.write(\'<script src="{CYTOSCAPE_JS_CDN[1]}"><\\/script>\')</script>')
    meta = {"node_labels": bool(node_labels), "legend_html": _legend_html(net, leg), "file": Path(path).stem}
    safe = lambda obj: json.dumps(obj).replace("</", "<\\/")
    doc = (_HTML.replace("__TITLE__", html.escape(title)).replace("__SCRIPT__", script)
           .replace("__ELEMENTS__", safe(el["nodes"] + el["edges"])).replace("__META__", safe(meta)))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(doc)
    return path


# ---------------------------------------------------------------------------


def save_all(net, outdir, prefix: str = "network", formats: Sequence[str] = ("png", "pdf", "html", "cyjs", "graphml", "csv"),
             title: Optional[str] = None, color_by: str = "auto", js_path: Optional[str] = None, **plot_kwargs) -> Dict[str, str]:
    """Save a network in several formats; returns ``{format: path}``."""
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, prefix)
    out = {}
    img = [f for f in formats if f in ("png", "pdf", "svg", "jpg", "tif", "tiff", "eps")]
    if img:
        fig = net.plot(None, color_by=color_by, title=title, **plot_kwargs)
        for f in img:
            p = f"{base}.{f}"
            fig.savefig(p, dpi=plot_kwargs.get("dpi", 300), facecolor="white")
            out[f] = p
        plt.close(fig)
    if "html" in formats:
        out["html"] = write_html(net, f"{base}.html", title=title, color_by=color_by, js_path=js_path)
    if "cyjs" in formats:
        out["cyjs"] = write_cyjs(net, f"{base}.cyjs", color_by=color_by)
        out["cytoscape_style"] = f"{base}_style.xml"
    if "graphml" in formats:
        out["graphml"] = write_graphml(net, f"{base}.graphml", color_by=color_by)
    if "csv" in formats:
        net.node_table().to_csv(f"{base}_nodes.csv", index=False)
        net.edge_table().to_csv(f"{base}_edges.csv", index=False)
        out["nodes"], out["edges"] = f"{base}_nodes.csv", f"{base}_edges.csv"
        if net.kind == "pathway":
            net.cluster_table().to_csv(f"{base}_clusters.csv", index=False)
            out["clusters"] = f"{base}_clusters.csv"
    return out
