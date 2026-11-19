"""Interactive HTML version of the pathway map (plotly.js).

Hover a point for the pathway, its region, set size and, for enriched
pathways, condition, FDR and overlapping genes. Region outlines (convex hulls)
show the group name on hover. Lasso or box-select ("circle") any pathways and
the side panel lists the regions they belong to, with counts, and the pathway
names themselves. The file is self-contained apart from plotly.js, which is
loaded from a CDN unless ``js_path`` embeds a local copy.
"""
from __future__ import annotations

import html
import json
from collections import defaultdict
from typing import Dict, Optional

import numpy as np
import pandas as pd

from .network import clean_pathway_name

PLOTLY_JS_CDN = [
    "https://cdn.plot.ly/plotly-2.35.2.min.js",
    "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js",
]


def _hull(xy: np.ndarray, core: float = 0.6):
    """Closed outline of the dense core of a region: the convex hull of the ``core``
    fraction of points closest to the region's median (regions are compact in
    similarity space, so a plain hull would span half the map). Small groups get a padded box."""
    if len(xy) >= 8 and core < 1:
        d = np.linalg.norm(xy - np.median(xy, axis=0), axis=1)
        xy = xy[d <= np.quantile(d, core)]
    if len(xy) < 3:
        x0, y0 = xy.min(axis=0)
        x1, y1 = xy.max(axis=0)
        pad = 0.02 * max(x1 - x0, y1 - y0, 1.0) + 0.5
        return _chaikin(np.array([[x0 - pad, y0 - pad], [x1 + pad, y0 - pad], [x1 + pad, y1 + pad], [x0 - pad, y1 + pad]]), passes=3)
    try:
        from scipy.spatial import ConvexHull

        h = ConvexHull(xy)
        pts = xy[h.vertices]
    except Exception:  # collinear points
        pts = xy
    c = pts.mean(axis=0)
    pts = c + (pts - c) * 1.08  # padding, so the rounded curve still encloses the points
    return _chaikin(pts, passes=3)


def _chaikin(pts: np.ndarray, passes: int = 3):
    """Chaikin corner cutting of a closed polygon -> smooth closed curve (list of [x, y], closed)."""
    P = np.asarray(pts, dtype=float)
    for _ in range(passes):
        Q = np.roll(P, -1, axis=0)
        P = np.column_stack([np.ravel(np.column_stack([0.75 * P[:, 0] + 0.25 * Q[:, 0], 0.25 * P[:, 0] + 0.75 * Q[:, 0]])),
                             np.ravel(np.column_stack([0.75 * P[:, 1] + 0.25 * Q[:, 1], 0.25 * P[:, 1] + 0.75 * Q[:, 1]]))])
    ring = np.round(P, 3).tolist()
    ring.append(ring[0])
    return ring


def write_map_html(pmap: pd.DataFrame, path, table: Optional[pd.DataFrame] = None, fdr_cutoff: Optional[float] = 0.05,
                   title: Optional[str] = None, js_path: Optional[str] = None, max_genes: int = 12,
                   categories: Optional[Dict[str, str]] = None) -> str:
    """Write the interactive map; returns ``path``.

    Library points are colored by region (switchable to source collection, for
    a map built from several databases, or to plain grey). ``table`` is an
    enrichment or comparison table whose enriched pathways are drawn as larger
    points on top, one color per condition. Without it the bare library is drawn.
    """
    from .pathmap import map_hits
    from .style import condition_color

    pm = pmap.reset_index(drop=True)
    hits = map_hits(pm, table, fdr_cutoff) if table is not None else pm.iloc[0:0]
    pvals: Dict[tuple, float] = {}
    full_genes: Dict[tuple, str] = {}
    if table is not None and len(hits):
        t = table.copy()
        if "condition" not in t.columns:
            t["condition"] = ""
        for _, r in t.iterrows():
            key = (str(r["pathway"]), str(r["condition"]))
            if "pvalue" in t.columns and not pd.isna(r.get("pvalue")):
                pvals.setdefault(key, float(r["pvalue"]))
            if "genes" in t.columns and isinstance(r.get("genes"), str):
                full_genes.setdefault(key, r["genes"])
    conds = list(hits.attrs.get("conditions", [""])) if table is not None and len(hits) else []
    label = {p: clean_pathway_name(p) for p in pm["pathway"]}
    if len(hits):
        hits = hits.merge(pm[["pathway", "set_size"]], on="pathway", how="left")
    region_label = pm.groupby("cluster")["cluster_label"].first().to_dict()
    region_size = pm[pm["cluster"] > 0].groupby("cluster").size().to_dict()

    # library points
    lib = {
        "x": pm["x"].round(3).tolist(), "y": pm["y"].round(3).tolist(),
        "name": pm["pathway"].tolist(), "label": [label[p] for p in pm["pathway"]],
        "size": pm["set_size"].astype(int).tolist(), "region": pm["cluster"].astype(int).tolist(),
    }
    # regions: hull + centre
    regions = []
    for c, grp in pm[pm["cluster"] > 0].groupby("cluster"):
        if not region_label.get(c):
            continue
        xy = grp[["x", "y"]].to_numpy(dtype=float)
        regions.append({"id": int(c), "label": region_label[c], "n": int(len(grp)), "hull": _hull(xy),
                        "cx": float(np.median(xy[:, 0])), "cy": float(np.median(xy[:, 1]))})
    # hits per condition
    hit_traces = []
    if len(hits):
        sig = -np.log10(hits["fdr"].clip(lower=1e-300))
        vmax = max(float(sig.max()), 1.0)
        for i, cond in enumerate(conds):
            sub = hits[hits["condition"].astype(str) == str(cond)]
            if sub.empty:
                continue
            keys = [(str(pth), str(cond)) for pth in sub["pathway"]]
            genes = pd.Series([full_genes.get(k, "") for k in keys], index=sub.index)
            short = [";".join(g.split(";")[:max_genes]) + (" …" if len(g.split(";")) > max_genes else "") for g in genes]
            hit_traces.append({
                "condition": str(cond) or "enriched", "color": condition_color(i),
                "x": sub["x"].round(3).tolist(), "y": sub["y"].round(3).tolist(),
                "name": sub["pathway"].tolist(), "label": [label.get(p, clean_pathway_name(p)) for p in sub["pathway"]],
                "fdr": [float(f"{v:.3g}") for v in sub["fdr"]], "region": sub["cluster"].astype(int).tolist(),
                "direction": sub["direction"].fillna("").astype(str).tolist() if "direction" in sub else [""] * len(sub),
                "overlap": sub["overlap"].fillna(0).astype(int).tolist() if "overlap" in sub else [0] * len(sub),
                "genes": short, "allGenes": [g.replace(";", ", ") for g in genes],
                "pvalue": [float(f"{pvals[k]:.3g}") if k in pvals else None for k in keys],
                "setSize": sub["set_size"].astype(int).tolist() if "set_size" in sub else [0] * len(sub),
                "msize": (7 + 13 * (-np.log10(sub["fdr"].clip(lower=1e-300)) / vmax)).round(1).tolist(),
            })
    from .style import lighten

    palette = [lighten(condition_color(i), 0.15) for i in range(16)] + [lighten(condition_color(i), 0.45) for i in range(16)]
    region_color = {str(r["id"]): palette[i % len(palette)] for i, r in enumerate(regions)}
    for r in regions:
        r["color"] = region_color[str(r["id"])]
    from .network import pathway_category

    cats = [pathway_category(p, categories) for p in pm["pathway"]]
    cat_order = list(pd.Series(cats).value_counts().index)
    lib["category"] = cats
    cat_color = {c: condition_color(i) for i, c in enumerate(cat_order)}
    data = {"library": lib, "regions": regions, "hits": hit_traces, "regionLabel": {str(k): v for k, v in region_label.items()},
            "regionColor": region_color, "catColor": cat_color, "multiCategory": len(cat_order) > 1,
            "regionSize": {str(k): int(v) for k, v in region_size.items()},
            "method": str(pm.attrs.get("method", "")).upper(), "n": int(len(pm)), "fdr": fdr_cutoff}
    title = title or (f"Pathway map, {len(pm)} gene sets" + (f" ({data['method']})" if data["method"] else ""))

    if js_path:
        with open(js_path, encoding="utf-8") as fh:
            script = "<script>" + fh.read().replace("</script", "<\\/script") + "</script>"
    else:
        script = (f'<script src="{PLOTLY_JS_CDN[0]}"></script>\n'
                  f'<script>window.Plotly || document.write(\'<script src="{PLOTLY_JS_CDN[1]}"><\\/script>\')</script>')
    doc = (_HTML.replace("__TITLE__", html.escape(title)).replace("__SCRIPT__", script)
           .replace("__DATA__", json.dumps(data).replace("</", "<\\/")))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(doc)
    return str(path)


_HTML = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>__TITLE__</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
__SCRIPT__
<style>
  :root{--bg:#fff;--line:#e3e3e3;--muted:#6b6b6b;--accent:#1B6AC9}
  *{box-sizing:border-box}
  body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Arial,sans-serif;color:#222;background:var(--bg);font-size:13px}
  #top{display:flex;align-items:center;gap:12px;height:44px;padding:0 14px;border-bottom:1px solid var(--line);background:#fafafa}
  #top h1{font-size:15px;font-weight:600;margin:0;flex:1;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
  #top .stat{color:var(--muted)}
  #top button{font-size:12px;padding:3px 10px;border:1px solid #c9c9c9;border-radius:4px;background:#fff;cursor:pointer}
  #top button:hover{background:#f0f0f0}
  #wrap{display:flex;height:calc(100vh - 44px)}
  #left{width:270px;border-right:1px solid var(--line);overflow:auto;background:#fcfcfc}
  #map{flex:1;min-width:0}
  #right{width:330px;border-left:1px solid var(--line);overflow:auto;padding:10px 12px}
  .panel{border-bottom:1px solid var(--line)}
  .panel > .hd{display:flex;align-items:center;gap:6px;padding:8px 10px;cursor:pointer;user-select:none;font-weight:600;font-size:12px;text-transform:uppercase;letter-spacing:.04em;color:#444}
  .panel > .hd .tri{display:inline-block;width:0;height:0;border-left:5px solid #888;border-top:4px solid transparent;border-bottom:4px solid transparent;transition:transform .15s}
  .panel.open > .hd .tri{transform:rotate(90deg)}
  .panel > .bd{display:none;padding:2px 6px 8px 10px}
  .panel.open > .bd{display:block}
  .row{display:flex;align-items:center;gap:6px;padding:2px 4px;border-radius:3px;cursor:pointer;line-height:1.3}
  .row:hover{background:#eef3fb}
  .row.active{background:#dfe9f8}
  .row .sw{width:10px;height:10px;border-radius:50%;flex:none;border:1px solid rgba(0,0,0,.25)}
  .row .lbl{flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
  .row .n{color:var(--muted);font-variant-numeric:tabular-nums;flex:none}
  .row input[type=checkbox]{margin:0;flex:none}
  .bd label.opt{display:flex;align-items:center;gap:6px;padding:2px 4px;cursor:pointer}
  .bd input[type=search],.bd select{width:100%;font-size:12px;padding:4px 6px;border:1px solid #ccc;border-radius:4px;margin:3px 0}
  #right h2{font-size:12px;margin:10px 0 4px;color:#555;font-weight:600;text-transform:uppercase;letter-spacing:.04em}
  #right h3{font-size:14px;margin:0 0 4px;font-weight:600}
  #right .hint{color:var(--muted)}
  #right .stat{color:var(--muted);font-size:12px}
  #right table{border-collapse:collapse;width:100%}
  #right td{padding:2px 4px;vertical-align:top;border-bottom:1px solid #f0f0f0}
  #right td.n{text-align:right;white-space:nowrap;color:var(--muted)}
  #right ul{margin:4px 0 0 16px;padding:0}
  #right li{margin:1px 0}
  #right details{margin:3px 0 6px}
  #right summary{cursor:pointer}
  #right .genes{color:#444;font-size:12px;word-break:break-word;margin:3px 0 0 0;line-height:1.4}
  #right button{font-size:11px;padding:1px 8px;margin-left:6px;border:1px solid #c9c9c9;border-radius:4px;background:#fff;cursor:pointer}
  #right a{color:var(--accent);text-decoration:none}
  .chip{display:inline-block;padding:0 6px;border-radius:9px;background:#eee;font-size:11px;color:#444;margin-left:4px}
</style></head>
<body>
<div id="top"><h1>__TITLE__</h1>
  <span class="stat" id="counts"></span>
  <button id="btnReset" title="clear selection, highlights, search and zoom">reset</button>
</div>
<div id="wrap">
  <div id="left">
    <div class="panel open"><div class="hd"><span class="tri"></span>Display</div><div class="bd">
      <label class="opt"><input type="checkbox" id="cbLib"> all library points <span class="stat" style="color:#888">(off = significant only)</span></label>
      <label class="opt"><input type="checkbox" id="cbRegions" checked> region outlines</label>
      <label class="opt"><input type="checkbox" id="cbLabels" checked> region names</label>
      <div style="margin-top:6px;color:#555">Color points by</div>
      <label class="opt"><input type="radio" name="col" value="none" checked> none (grey)</label>
      <label class="opt"><input type="radio" name="col" value="region"> region</label>
      <label class="opt" id="optCat"><input type="radio" name="col" value="category"> collection</label>
    </div></div>
    <div class="panel open"><div class="hd"><span class="tri"></span>Find pathway</div><div class="bd">
      <select id="selGroup" title="narrow the list to a region, collection or condition"><option value="">all groups</option></select>
      <select id="selPathway" title="pick a pathway"><option value="">pathway…</option></select>
      <input type="search" id="find" list="pathwayList" placeholder="or type part of a name"><datalist id="pathwayList"></datalist>
    </div></div>
    <div class="panel open" id="pnlCond"><div class="hd"><span class="tri"></span>Conditions <span class="chip" id="nCond"></span></div><div class="bd" id="condRows"></div></div>
    <div class="panel" id="pnlCat"><div class="hd"><span class="tri"></span>Collections <span class="chip" id="nCat"></span></div><div class="bd" id="catRows"></div></div>
    <div class="panel open"><div class="hd"><span class="tri"></span>Regions <span class="chip" id="nReg"></span></div><div class="bd" id="regRows"></div></div>
  </div>
  <div id="map"></div>
  <div id="right"><div id="sel"></div></div>
</div>
<script>
const D = __DATA__;
const RL = D.regionLabel, RS = D.regionSize, L = D.library;
const $ = id => document.getElementById(id);
function esc(s){return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;');}
function regionName(r){return r>0 && RL[String(r)] ? RL[String(r)] : 'no region';}
function hexa(h, a){const n=parseInt(h.slice(1),16); return `rgba(${n>>16},${(n>>8)&255},${n&255},${a})`;}
const GREY = '#C4C4C4';
function libColors(mode){
  if(mode==='region') return L.region.map(r => D.regionColor[String(r)] || GREY);
  if(mode==='category') return L.category.map(c => D.catColor[c] || GREY);
  return L.region.map(() => GREY);
}
const HINT = '<span class="hint">The map shows the significant pathways; tick "all library points" to see every gene set. Click a dot for its pathway. Hover an outline for the region name. Choose lasso or box in the map toolbar and circle pathways to list their regions, the significant ones in this data with genes and p-values, and all selected pathways here. Click a region or condition on the left to highlight it.</span>';
const traces = [];
D.regions.forEach(r => {
  traces.push({type:'scatter', mode:'lines', x:r.hull.map(p=>p[0]), y:r.hull.map(p=>p[1]), fill:'toself',
    fillcolor:hexa(r.color, 0.10), line:{color:hexa(r.color, 0.7), width:1.2}, hoveron:'fills',
    hoverinfo:'text', text:`<b>${esc(r.label)}</b><br>${r.n} gene sets`, name:r.label, showlegend:false, legendgroup:'regions',
    meta:{kind:'region', id:r.id}});
});
const nRegionTraces = traces.length;
const showLib = D.hits.length === 0;  // with results, the main view shows only the significant pathways
traces.push({type:'scatter', mode:'markers', x:L.x, y:L.y, name:'library', showlegend:false, visible: showLib,
  marker:{size:6, color:libColors('none'), opacity:0.85, line:{width:0}},
  customdata:L.name.map((n,i)=>[n, L.label[i], L.size[i], L.region[i], L.category[i]]),
  hovertemplate:'<b>%{customdata[1]}</b><br>%{customdata[0]}<br>%{customdata[2]} genes · region: '+'%{text}'+(D.multiCategory?'<br>%{customdata[4]}':'')+'<extra></extra>',
  text:L.region.map(regionName), meta:{kind:'library'}});
const LIB = nRegionTraces;
D.hits.forEach((h, k) => {
  traces.push({type:'scatter', mode:'markers', x:h.x, y:h.y, name:h.condition, legendgroup:h.condition, showlegend:false,
    marker:{size:h.msize, color:h.color, line:{color:'#222', width:1.2}, opacity:0.95},
    customdata:h.name.map((n,i)=>[n, h.label[i], h.fdr[i], regionName(h.region[i]), h.direction[i], h.overlap[i], h.genes[i], h.region[i]]),
    hovertemplate:'<b>%{customdata[1]}</b><br>%{customdata[0]}<br>'+esc(h.condition)+(h.direction.some(d=>d)?' %{customdata[4]}':'')+
      ' · FDR %{customdata[2]} · %{customdata[5]} genes<br>region: %{customdata[3]}<br><span style="color:#666">%{customdata[6]}</span><extra></extra>',
    meta:{kind:'hit', k}});
});
function wrap(s, w){const out=[]; let line=''; s.split(' ').forEach(t=>{ if((line+' '+t).trim().length>w){out.push(line.trim()); line=t;} else line+=' '+t; }); if(line.trim()) out.push(line.trim()); return out.join('<br>');}
const ann = D.regions.map(r => ({x:r.cx, y:r.cy, text:wrap(r.label, 18), showarrow:false, font:{size:10, color:'#444'}, bgcolor:'rgba(255,255,255,0.7)', borderpad:2}));
const layout = {margin:{l:10,r:10,t:10,b:10}, hovermode:'closest', dragmode:'pan', showlegend:false, annotations:ann,
  xaxis:{visible:false, scaleanchor:'y'}, yaxis:{visible:false}, plot_bgcolor:'#fff', paper_bgcolor:'#fff'};
const cfg = {responsive:true, displaylogo:false, modeBarButtonsToRemove:['autoScale2d'], scrollZoom:true,
  toImageButtonOptions:{format:'png', filename:'pathway_map', scale:2}};
const el = $('map');
function repelAnnotations(){
  const fl = el._fullLayout; if(!fl || !fl.xaxis || !$('cbLabels').checked) return;
  const xa = fl.xaxis, ya = fl.yaxis, W = fl.width, H = fl.height;
  const boxes = ann.map(a => { const lines = a.text.split('<br>'); const w = Math.max(...lines.map(l=>l.length))*5.6+6, h = lines.length*12+4;
    return {a, x: xa.l2p(a.x)+fl._size.l, y: fl._size.t + ya.l2p(a.y), w, h}; });
  for(let it=0; it<200; it++){ let moved=false;
    for(let i=0;i<boxes.length;i++) for(let j=i+1;j<boxes.length;j++){ const A=boxes[i], B=boxes[j];
      const dx=(A.w+B.w)/2 - Math.abs(A.x-B.x), dy=(A.h+B.h)/2 - Math.abs(A.y-B.y); if(dx<=0 || dy<=0) continue; moved=true;
      if(dy*1.5<=dx){ const s=(A.y>=B.y?1:-1)*dy/2*0.7; A.y+=s; B.y-=s; } else { const s=(A.x>=B.x?1:-1)*dx/2*0.7; A.x+=s; B.x-=s; } }
    if(!moved) break; }
  const out = boxes.map(b => Object.assign({}, b.a, {x: xa.p2l(Math.min(Math.max(b.x-fl._size.l, b.w/2), W-b.w/2)), y: ya.p2l(Math.min(Math.max(b.y-fl._size.t, b.h/2), H-b.h/2))}));
  Plotly.relayout(el, {annotations: out});
}
let relayoutBusy = false;
Plotly.newPlot(el, traces, layout, cfg).then(() => {
  el.on('plotly_selected', ev => showSelection(ev));
  el.on('plotly_deselect', () => resetSel());
  el.on('plotly_click', ev => { const p = ev.points && ev.points[0]; if(!p || !p.data.meta) return;
    if(p.data.meta.kind === 'library') highlight([p.pointIndex], false, true);
    else if(p.data.meta.kind === 'hit'){ const i = L.name.indexOf(p.customdata[0]); if(i >= 0) highlight([i], false, true); }
    else if(p.data.meta.kind === 'region') showRegion(p.data.meta.id); });
  el.on('plotly_relayout', ev => { if(relayoutBusy || !ev || !Object.keys(ev).some(k=>k.startsWith('xaxis')||k.startsWith('yaxis')||k==='autosize')) return;
    relayoutBusy = true; Plotly.relayout(el, {annotations: ann}).then(() => { repelAnnotations(); relayoutBusy = false; }); });
  relayoutBusy = true; repelAnnotations(); relayoutBusy = false;
});
function resetSel(){ $('sel').innerHTML = HINT; setCounts(0); }
function setCounts(nSel){ const nSig = new Set(D.hits.flatMap(h=>h.name)).size;
  $('counts').textContent = `${D.n} gene sets${D.method?' · '+D.method:''}${D.hits.length?' · '+nSig+' significant':''}${nSel?' · '+nSel+' selected':''}`; }
function resetAll(){
  highlight([]); clearActive();
  $('find').value = ''; $('selGroup').value = ''; fillPathways();
  Plotly.restyle(el, {selectedpoints: [null]});
  relayoutBusy = true;
  Plotly.relayout(el, {'xaxis.autorange': true, 'yaxis.autorange': true, selections: [], dragmode: 'pan', annotations: $('cbLabels').checked ? ann : []})
    .then(() => { repelAnnotations(); relayoutBusy = false; });
  resetSel();
}
function sigRows(sig){
  return sig.map(h => `<details><summary><b>${esc(h.label)}</b> <span class="stat">${esc(h.cond)}${h.dir?' '+esc(h.dir):''} · FDR ${h.fdr}${h.p!==null&&h.p!==undefined?' · p '+h.p:''} · ${h.overlap}/${h.size} genes</span></summary>
      <div class="stat">${esc(h.name)} · region: ${esc(h.region)}</div><div class="genes">${esc(h.genes||'(gene list not available)')}</div></details>`).join('');
}
function hitRecord(h, i){ return {name:h.name[i], label:h.label[i], cond:h.condition, fdr:h.fdr[i], p:h.pvalue[i], overlap:h.overlap[i], size:h.setSize[i], dir:h.direction[i], genes:h.allGenes[i], region:regionName(h.region[i])}; }
function showSelection(ev){
  if(!ev || !ev.points || !ev.points.length){resetSel(); return;}
  const perRegion = {}, names = new Map(), sig = [];
  ev.points.forEach(p => {
    const kind = p.data.meta && p.data.meta.kind; if(kind==='region' || kind==='found') return;
    const cd = p.customdata; const reg = kind==='hit' ? cd[7] : cd[3];
    if(!names.has(cd[0])) perRegion[String(reg)] = (perRegion[String(reg)]||0)+1;
    names.set(cd[0], names.get(cd[0]) || {label:cd[1]});
    if(kind==='hit') sig.push(hitRecord(D.hits[p.data.meta.k], p.pointIndex));
  });
  sig.sort((a,b)=>a.fdr-b.fdr);
  const rows = Object.entries(perRegion).sort((a,b)=>b[1]-a[1]).map(([k,n]) =>
    `<tr><td><a href="#" onclick="showRegion(${+k});return false;">${esc(regionName(+k))}</a></td><td class="n">${n}${RS[k]?' / '+RS[k]:''}</td></tr>`).join('');
  const sigNames = new Set(sig.map(h=>h.name));
  const others = [...names.entries()].filter(([n])=>!sigNames.has(n)).slice(0,400).map(([n,v]) => `<li><a href="#" onclick="highlight([${L.name.indexOf(n)}],false,true);return false;">${esc(v.label)}</a></li>`).join('');
  window._selTSV = ['pathway\tcondition\tdirection\tFDR\tp\toverlap\tset_size\tregion\tgenes'].concat(sig.map(h=>[h.name,h.cond,h.dir,h.fdr,h.p,h.overlap,h.size,h.region,h.genes].join('\t'))).join('\n');
  setCounts(names.size);
  $('sel').innerHTML =
    `<h3>${names.size} pathways selected <span class="stat">in ${Object.keys(perRegion).length} regions</span><button onclick="resetAll()">reset</button></h3>
     <h2>Regions (selected / in region)</h2><table>${rows}</table>
     <h2>Significant in this data (${sig.length})${sig.length?'<button onclick="copySel()">copy table</button>':''}</h2>${sig.length ? sigRows(sig) : '<span class="hint">none of the selected pathways is significant in this data</span>'}
     <details><summary><b>Other selected pathways (${names.size - sigNames.size})</b> <span class="stat">not significant in this data</span></summary><ul>${others}</ul>${names.size - sigNames.size>400?'<div class="hint">(first 400 shown)</div>':''}</details>`;
}
function copySel(){ navigator.clipboard && navigator.clipboard.writeText(window._selTSV||'').then(()=>{}, ()=>{}); }
// ---- left sidebar panels
document.querySelectorAll('.panel > .hd').forEach(h => h.onclick = () => h.parentElement.classList.toggle('open'));
let activeRow = null;
function clearActive(){ if(activeRow){ activeRow.classList.remove('active'); activeRow = null; } }
function rowEl(color, label, n, onClick, withCheck){
  const d = document.createElement('div'); d.className = 'row';
  d.innerHTML = (withCheck ? '<input type="checkbox" checked>' : '') + `<span class="sw" style="background:${color}"></span><span class="lbl" title="${esc(label)}">${esc(label)}</span><span class="n">${n}</span>`;
  d.querySelector('.lbl').onclick = () => { if(activeRow === d){ clearActive(); highlight([]); resetSel(); return; } clearActive(); activeRow = d; d.classList.add('active'); onClick(); };
  return d;
}
(function buildPanels(){
  // conditions
  if(D.hits.length){ $('nCond').textContent = D.hits.length;
    D.hits.forEach((h, k) => { const r = rowEl(h.color, h.condition, h.name.length, () => showCondition(k), true);
      r.querySelector('input').onchange = e => Plotly.restyle(el, {visible: e.target.checked}, [LIB + 1 + k]); $('condRows').appendChild(r); });
  } else $('pnlCond').style.display = 'none';
  // collections
  if(D.multiCategory){ const counts = {}; L.category.forEach(c => counts[c] = (counts[c]||0)+1); $('nCat').textContent = Object.keys(counts).length;
    Object.entries(counts).sort((a,b)=>b[1]-a[1]).forEach(([c,n]) => $('catRows').appendChild(rowEl(D.catColor[c]||GREY, c, n, () => showCategory(c), false)));
  } else { $('pnlCat').style.display = 'none'; $('optCat').style.display = 'none'; }
  // regions, largest first
  $('nReg').textContent = D.regions.length;
  [...D.regions].sort((a,b)=>b.n-a.n).forEach(r => $('regRows').appendChild(rowEl(r.color, r.label, r.n, () => showRegion(r.id), false)));
})();
function sigIn(pred){ const out = []; D.hits.forEach(h => h.name.forEach((n,i) => { if(pred(n, h, i)) out.push(hitRecord(h, i)); })); return out.sort((a,b)=>a.fdr-b.fdr); }
function showRegion(id){
  const r = D.regions.find(x => x.id === id); if(!r) return;
  const idx = [...L.name.keys()].filter(i => L.region[i] === id);
  const sig = sigIn((n, h, i) => h.region[i] === id); const sigSet = new Set(sig.map(x=>x.name));
  ringPoints(idx, r.color);
  relayoutBusy = true; zoomTo(idx);
  $('sel').innerHTML = `<h3>${esc(r.label)} <span class="stat">region</span><button onclick="resetAll()">reset</button></h3>
    <div class="stat">${r.n} gene sets${D.hits.length ? ' · ' + new Set(sig.map(s=>s.name)).size + ' significant in this data' : ''}</div>
    ${D.hits.length ? '<h2>Significant in this region ('+sig.length+')</h2>' + (sig.length ? sigRows(sig) : '<span class="hint">none</span>') : ''}
    <details${D.hits.length?'':' open'}><summary><b>${D.hits.length ? 'Other pathways in region' : 'Pathways in region'} (${idx.length - new Set(sig.map(x=>x.name)).size})</b></summary><ul>${idx.filter(i => !sigSet.has(L.name[i])).sort((a,b)=>L.label[a].localeCompare(L.label[b])).map(i => `<li><a href="#" onclick="highlight([${i}],false,true);return false;">${esc(L.label[i])}</a></li>`).join('')}</ul></details>`;
  setCounts(idx.length);
}
function showCondition(k){
  const h = D.hits[k]; const sig = h.name.map((n,i) => hitRecord(h, i)).sort((a,b)=>a.fdr-b.fdr);
  const perRegion = {}; h.region.forEach(r => perRegion[String(r)] = (perRegion[String(r)]||0)+1);
  highlight([]); ringPoints(h.name.map(n => L.name.indexOf(n)).filter(i=>i>=0), h.color);
  const rows = Object.entries(perRegion).sort((a,b)=>b[1]-a[1]).map(([r,n]) => `<tr><td><a href="#" onclick="showRegion(${+r});return false;">${esc(regionName(+r))}</a></td><td class="n">${n}${RS[r]?' / '+RS[r]:''}</td></tr>`).join('');
  window._selTSV = ['pathway\tcondition\tdirection\tFDR\tp\toverlap\tset_size\tregion\tgenes'].concat(sig.map(s=>[s.name,s.cond,s.dir,s.fdr,s.p,s.overlap,s.size,s.region,s.genes].join('\t'))).join('\n');
  $('sel').innerHTML = `<h3>${esc(h.condition)} <span class="stat">condition</span><button onclick="resetAll()">reset</button></h3>
    <div class="stat">${sig.length} significant pathways (FDR ≤ ${D.fdr})</div><h2>By region</h2><table>${rows}</table>
    <h2>Pathways<button onclick="copySel()">copy table</button></h2>${sigRows(sig)}`;
  setCounts(sig.length);
}
function showCategory(c){
  const idx = [...L.name.keys()].filter(i => L.category[i] === c); const sig = sigIn((n, h, i) => L.category[L.name.indexOf(n)] === c);
  ringPoints(idx, D.catColor[c] || GREY);
  const sigSet = new Set(sig.map(x=>x.name));
  $('sel').innerHTML = `<h3>${esc(c)} <span class="stat">collection</span><button onclick="resetAll()">reset</button></h3><div class="stat">${idx.length} gene sets on the map${D.hits.length?' · '+sigSet.size+' significant':''}</div>
    ${D.hits.length ? '<h2>Significant in this data ('+sig.length+')</h2>' + (sig.length ? sigRows(sig) : '<span class="hint">none</span>') : ''}
    <details><summary><b>Other pathways in collection (${idx.length - sigSet.size})</b></summary><ul>${idx.filter(i => !sigSet.has(L.name[i])).sort((a,b)=>L.label[a].localeCompare(L.label[b])).slice(0,2000).map(i => `<li><a href="#" onclick="highlight([${i}],false,true);return false;">${esc(L.label[i])}</a></li>`).join('')}</ul></details>`;
  setCounts(idx.length);
}
// ---- find
const selGroup = $('selGroup'), selPathway = $('selPathway'), findBox = $('find'), dl = $('pathwayList');
(function fillGroups(){
  const regs = D.regions.map(r => ({v:'r:'+r.id, t:`region: ${r.label} (${r.n})`}));
  const cats = D.multiCategory ? Object.keys(D.catColor).map(c => ({v:'c:'+c, t:`collection: ${c}`})) : [];
  const conds = D.hits.map(h => ({v:'h:'+h.condition, t:`significant in ${h.condition} (${h.name.length})`}));
  [...conds, ...cats, ...regs].forEach(o => { const e=document.createElement('option'); e.value=o.v; e.textContent=o.t; selGroup.appendChild(e); });
})();
function pathwaysInGroup(g){
  let idx = [...L.name.keys()];
  if(g.startsWith('r:')) idx = idx.filter(i => String(L.region[i]) === g.slice(2));
  else if(g.startsWith('c:')) idx = idx.filter(i => L.category[i] === g.slice(2));
  else if(g.startsWith('h:')){ const h = D.hits.find(t => t.condition === g.slice(2)); const set = new Set(h ? h.name : []); idx = idx.filter(i => set.has(L.name[i])); }
  return idx.sort((a,b) => L.label[a].localeCompare(L.label[b]));
}
function fillPathways(){
  const idx = pathwaysInGroup(selGroup.value);
  selPathway.innerHTML = '<option value="">pathway… (' + idx.length + ')</option>' + idx.map(i => `<option value="${esc(L.name[i])}">${esc(L.label[i])}</option>`).join('');
  dl.innerHTML = idx.slice(0, 3000).map(i => `<option value="${esc(L.label[i])}"></option>`).join('');
}
fillPathways();
selGroup.onchange = () => { fillPathways(); findBox.value=''; highlight([]); };
selPathway.onchange = () => { findBox.value = ''; const i = L.name.indexOf(selPathway.value); highlight(i >= 0 ? [i] : [], true, true); };
let findTimer;
findBox.oninput = () => { clearTimeout(findTimer); findTimer = setTimeout(() => {
  const q = findBox.value.trim().toLowerCase(); selPathway.value = '';
  if(!q){ highlight([]); return; }
  highlight(pathwaysInGroup(selGroup.value).filter(i => L.name[i].toLowerCase().includes(q) || L.label[i].toLowerCase().includes(q)), false, true);
}, 250); };
// ---- highlight / details
let ringOn = false;
function ringPoints(idx, color, withText){
  if(ringOn){ Plotly.deleteTraces(el, el.data.length-1); ringOn = false; }
  if(!idx.length) return;
  Plotly.addTraces(el, {type:'scatter', mode: withText ? 'markers+text' : 'markers', x:idx.map(i=>L.x[i]), y:idx.map(i=>L.y[i]), text:idx.map(i=>L.label[i]),
    textposition:'top center', textfont:{size:10, color:'#B00'}, marker:{size: idx.length > 50 ? 9 : 14, color:'rgba(0,0,0,0)', line:{color: color || '#E41A1C', width:2}},
    name:'highlight', hoverinfo:'skip', meta:{kind:'found'}});
  ringOn = true;
}
function zoomTo(idx){
  if(!idx.length){ relayoutBusy = false; return; }
  const xs = idx.map(i=>L.x[i]), ys = idx.map(i=>L.y[i]);
  const span = Math.max(Math.max(...xs)-Math.min(...xs), Math.max(...ys)-Math.min(...ys), 1) * 0.6 + 2;
  const cx = (Math.max(...xs)+Math.min(...xs))/2, cy = (Math.max(...ys)+Math.min(...ys))/2;
  Plotly.relayout(el, {'xaxis.range':[cx-span, cx+span], 'yaxis.range':[cy-span, cy+span]}).then(() => { repelAnnotations(); relayoutBusy = false; });
}
function highlight(idx, zoom, details){
  ringPoints(idx, '#E41A1C', idx.length <= 50);
  if(!idx.length) return;
  if(idx.length === 1){ showPathwayDetails(idx[0]); if(zoom){ relayoutBusy = true; zoomTo(idx); } }
  else if(details) showList(idx);
}
function showPathwayDetails(i){
  const n = L.name[i]; const found = sigIn((hn) => hn === n);
  const stats = found.length ? found.map(f => `<div><b>${esc(f.cond)}</b>${f.dir?' '+esc(f.dir):''}: FDR ${f.fdr}${f.p!==null&&f.p!==undefined?' · p '+f.p:''} · ${f.overlap}/${f.size} genes<div class="genes">${esc(f.genes||'')}</div></div>`).join('')
    : '<span class="hint">not significant in this data</span>';
  $('sel').innerHTML = `<h3>${esc(L.label[i])}<button onclick="resetAll()">reset</button></h3><div class="stat">${esc(n)}</div>
    <table><tr><td>Gene-set size</td><td class="n">${L.size[i]}</td></tr><tr><td>Region</td><td class="n"><a href="#" onclick="showRegion(${L.region[i]});return false;">${esc(regionName(L.region[i]))}</a></td></tr>${D.multiCategory?`<tr><td>Collection</td><td class="n">${esc(L.category[i])}</td></tr>`:''}</table>
    <h2>In this data</h2>${stats}`;
  setCounts(1);
}
function showList(idx){
  const sigSet = new Set(D.hits.flatMap(h=>h.name));
  const li = i => `<li><a href="#" onclick="highlight([${i}],true,true);return false;">${esc(L.label[i])}</a></li>`;
  const sigIdx = idx.filter(i => sigSet.has(L.name[i])), rest = idx.filter(i => !sigSet.has(L.name[i]));
  $('sel').innerHTML = `<h3>${idx.length} matching pathways<button onclick="resetAll()">reset</button></h3>` +
    (D.hits.length ? `<h2>Significant in this data (${sigIdx.length})</h2>` + (sigIdx.length ? '<ul>'+sigIdx.map(li).join('')+'</ul>' : '<span class="hint">none</span>') +
      `<details><summary><b>Other matches (${rest.length})</b></summary><ul>${rest.map(li).join('')}</ul></details>` : '<ul>'+idx.map(li).join('')+'</ul>');
  setCounts(idx.length);
}
// ---- display controls
$('btnReset').onclick = resetAll;
$('cbLib').checked = showLib;
$('cbRegions').onchange = e => Plotly.restyle(el, {visible: e.target.checked}, [...Array(nRegionTraces).keys()]);
$('cbLabels').onchange = e => { relayoutBusy = true; Plotly.relayout(el, {annotations: e.target.checked ? ann : []}).then(() => { if(e.target.checked) repelAnnotations(); relayoutBusy = false; }); };
$('cbLib').onchange = e => Plotly.restyle(el, {visible: e.target.checked}, [LIB]);
document.querySelectorAll('input[name=col]').forEach(r => r.onchange = e => Plotly.restyle(el, {'marker.color': [libColors(e.target.value)]}, [LIB]));
resetSel();
</script></body></html>
"""
