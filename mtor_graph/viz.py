"""Visualisations: static layered figure (matplotlib) and interactive HTML (vis-network).

Colour follows the accompanying visual-design description:
module groups take a fixed eight-colour categorical palette (never
cycled); edge polarity uses the reserved status colours (green = activation,
red = inhibition, grey = neutral) reinforced by arrowhead shape; evidence level
is encoded by line style (solid / dashed / dotted). Every node is directly
labelled, so colour never carries identity alone.
"""

from __future__ import annotations

import json
from html import escape
import math
import os
import textwrap
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import networkx as nx

from .loader import to_simple

# Fixed categorical order for light and dark themes; see docs/visual-design.md
PALETTE_LIGHT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
PALETTE_DARK = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"]

GROUP_LABELS: dict[str, str] = {
    "nutrient_sensing": "Lysosomal nutrient sensing",
    "upstream_inputs": "Growth-factor, energy and stress inputs",
    "core_complexes": "Core complexes (mTORC1 / mTORC2)",
    "mtorc1_outputs": "mTORC1 outputs",
    "mtorc2": "mTORC2 regulation and outputs",
    "circadian": "Circadian coupling",
    "tissue_disease": "Tissue physiology and disease",
    "pharmacology": "Pharmacology",
}

EFFECT_COLOURS_LIGHT = {"activates": "#0ca30c", "required_for": "#0ca30c", "inhibits": "#d03b3b", "neutral": "#898781"}
EFFECT_COLOURS_DARK = {"activates": "#0ca30c", "required_for": "#0ca30c", "inhibits": "#e66767", "neutral": "#898781"}

SURFACE = {"light": "#fcfcfb", "dark": "#1a1a19"}
INK = {"light": "#0b0b0b", "dark": "#ffffff"}
INK2 = {"light": "#52514e", "dark": "#c3c2b7"}

MARKERS = {
    "protein": "o",
    "complex": "s",
    "metabolite": "D",
    "lipid": "D",
    "drug": "h",
    "process": "o",
    "stimulus": "^",
    "phenotype": "*",
    "tissue": "s",
    "disease": "*",
    "organelle": "s",
    "rna": "v",
    "transcription_program": "o",
    "gene_variant": "v",
}


def group_of(module: str, modules_meta: dict[str, dict[str, Any]] | None) -> str:
    if modules_meta and module in modules_meta:
        return modules_meta[module].get("group", "core_complexes")
    fallback = {
        "lysosomal_nutrient_sensing": "nutrient_sensing",
        "growth_factor_input": "upstream_inputs",
        "energy_stress_input": "upstream_inputs",
        "core_complexes": "core_complexes",
        "mtorc1_output_translation_ribosome": "mtorc1_outputs",
        "mtorc1_output_autophagy_metabolism_feedback": "mtorc1_outputs",
        "mtorc2_regulation_and_outputs": "mtorc2",
        "circadian": "circadian",
        "tissue_physiology_disease": "tissue_disease",
        "pharmacology": "pharmacology",
    }
    return fallback.get(module, "core_complexes")


def group_colour(group: str, dark: bool = False) -> str:
    pal = PALETTE_DARK if dark else PALETTE_LIGHT
    return pal[list(GROUP_LABELS).index(group)] if group in GROUP_LABELS else "#898781"


# --------------------------------------------------------------------------- layout

ROW_CAPACITY = 28  # nodes per row before a layer is split into staggered sub-rows


def layered_positions(
    G: nx.MultiDiGraph, x_gap: float = 1.0, y_gap: float = 2.2, sweeps: int = 8, row_capacity: int = ROW_CAPACITY
) -> tuple[dict[str, tuple[float, float]], dict[int, float]]:
    """Sugiyama-style layered layout.

    Returns ``(pos, layer_y)``: node positions and the y coordinate of every
    layer present in the graph. Empty layers are skipped (the vertical rank is
    the order among layers that actually occur), and layers wider than
    ``row_capacity`` are split into staggered sub-rows so that labels stay
    legible. x order is refined by barycentre sweeps to limit crossings.
    """
    modules = G.graph.get("modules", {})
    order = list(GROUP_LABELS)
    layers: dict[int, list[str]] = defaultdict(list)
    for n, d in G.nodes(data=True):
        layers[int(d.get("layer", 4))].append(n)
    adj = {n: set(G.predecessors(n)) | set(G.successors(n)) for n in G}
    for l in layers:
        layers[l].sort(key=lambda n: (order.index(group_of(G.nodes[n].get("module", ""), modules)), G.nodes[n].get("label", n)))
    x: dict[str, float] = {}

    def assign(l: int) -> None:
        m = len(layers[l])
        for i, n in enumerate(layers[l]):
            x[n] = i - (m - 1) / 2

    for l in layers:
        assign(l)
    keys = sorted(layers)
    for s in range(sweeps):
        seq = keys if s % 2 == 0 else list(reversed(keys))
        for l in seq:
            def bary(n: str) -> float:
                nb = [x[m] for m in adj[n] if int(G.nodes[m].get("layer", 4)) != l]
                return sum(nb) / len(nb) if nb else x[n]
            layers[l].sort(key=lambda n: (bary(n), x[n]))
            assign(l)

    pos: dict[str, tuple[float, float]] = {}
    layer_y: dict[int, float] = {}
    y = 0.0
    for l in keys:
        members = layers[l]
        rows = max(1, min(3, math.ceil(len(members) / row_capacity)))
        layer_y[l] = y
        row_gap = y_gap * 0.38
        for i, n in enumerate(members):
            r = i % rows
            # spread sub-rows horizontally so that stacked nodes do not share x
            pos[n] = (x[n] * x_gap * (1.0 if rows == 1 else 1.0), y - r * row_gap)
        y -= y_gap + (rows - 1) * row_gap
    return pos, layer_y


# --------------------------------------------------------------------------- matplotlib

def draw_layered(
    G: nx.MultiDiGraph,
    path: str | os.PathLike,
    *,
    title: str | None = None,
    subtitle: str | None = None,
    dark: bool = False,
    collapse: bool = True,
    label_wrap: int = 16,
    font_size: float = 6.5,
    dpi: int = 200,
    show_sites: bool = False,
    figsize: tuple[float, float] | None = None,
) -> Path:
    """Draw the graph top-down by layer and save to ``path`` (png / svg / pdf)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    import matplotlib.patheffects as pe

    mode = "dark" if dark else "light"
    modules = G.graph.get("modules", {})
    pos, layer_y = layered_positions(G)
    H = to_simple(G) if collapse else G
    layer_sizes = defaultdict(int)
    for _, d in G.nodes(data=True):
        layer_sizes[int(d.get("layer", 4))] += 1
    max_n = max(layer_sizes.values()) if layer_sizes else 1
    y_span = (max(layer_y.values()) - min(layer_y.values())) if layer_y else 1.0
    if figsize is None:
        width = max(12.0, 0.7 * min(max_n, 3 * ROW_CAPACITY) + 2)
        figsize = (width, max(8.0, 1.05 * y_span + 4.5))
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    fig.patch.set_facecolor(SURFACE[mode])
    ax.set_facecolor(SURFACE[mode])
    ax.axis("off")

    deg = {n: G.in_degree(n) + G.out_degree(n) for n in G}
    eff_col = EFFECT_COLOURS_DARK if dark else EFFECT_COLOURS_LIGHT
    ls = {"established": "-", "emerging": (0, (4, 2)), "contested": (0, (1, 2))}

    # edges first (recessive)
    for u, v, d in H.edges(data=True):
        effect = d.get("effect", "modulates")
        colour = eff_col.get(effect, eff_col["neutral"])
        if effect in ("activates", "required_for"):
            style = "-|>,head_length=0.35,head_width=0.18"
        elif effect == "inhibits":
            style = "-[,widthB=0.45,lengthB=0.12"
        elif effect == "component_of":
            style = "-"
        else:
            style = "-|>,head_length=0.25,head_width=0.12"
        rad = 0.12 if pos[u][0] <= pos[v][0] else -0.12
        ax.annotate(
            "",
            xy=pos[v],
            xytext=pos[u],
            arrowprops=dict(
                arrowstyle=style,
                color=colour,
                lw=0.7,
                alpha=0.55 if effect != "component_of" else 0.35,
                connectionstyle=f"arc3,rad={rad}",
                shrinkA=6,
                shrinkB=6,
                linestyle=ls.get(d.get("evidence"), "-"),
            ),
            zorder=1,
        )
        if show_sites and not collapse and d.get("sites"):
            mx, my = (pos[u][0] + pos[v][0]) / 2, (pos[u][1] + pos[v][1]) / 2
            ax.text(mx, my, ", ".join(d["sites"]), fontsize=font_size * 0.75, color=INK2[mode], ha="center", va="center", zorder=2)

    # nodes grouped by (group, marker) for one scatter call each
    buckets: dict[tuple[str, str], list[str]] = defaultdict(list)
    for n, d in G.nodes(data=True):
        buckets[(group_of(d.get("module", ""), modules), MARKERS.get(d.get("type", "protein"), "o"))].append(n)
    for (g, marker), nodes in buckets.items():
        xs = [pos[n][0] for n in nodes]
        ys = [pos[n][1] for n in nodes]
        sizes = [min(320, 36 + 10 * deg[n]) for n in nodes]
        ax.scatter(xs, ys, s=sizes, marker=marker, c=group_colour(g, dark), edgecolors=SURFACE[mode], linewidths=0.8, zorder=3)

    # direct labels, staggered above / below to limit collisions
    halo = [pe.withStroke(linewidth=2.2, foreground=SURFACE[mode])]
    by_row: dict[float, list[str]] = defaultdict(list)
    for n in G.nodes:
        by_row[round(pos[n][1], 3)].append(n)
    for yv, members in by_row.items():
        members.sort(key=lambda n: pos[n][0])
        for i, n in enumerate(members):
            label = "\n".join(textwrap.wrap(G.nodes[n].get("label", n), label_wrap)) or n
            above = i % 2 == 0
            ax.text(
                pos[n][0],
                pos[n][1] + (0.2 if above else -0.2),
                label,
                fontsize=font_size,
                color=INK[mode],
                ha="center",
                va="bottom" if above else "top",
                zorder=4,
                path_effects=halo,
                linespacing=0.95,
            )

    # layer captions on the left
    from .schema import LAYERS

    from matplotlib.transforms import blended_transform_factory

    cap_tf = blended_transform_factory(fig.transFigure, ax.transData)  # x in figure fraction, y in data
    for l, yv in layer_y.items():
        ax.text(0.008, yv, "\n".join(textwrap.wrap(f"L{l}  {LAYERS.get(l, '')}", 22)), transform=cap_tf, fontsize=font_size + 1,
                color=INK2[mode], ha="left", va="center", style="italic", clip_on=False)

    # legend
    handles = [Patch(facecolor=group_colour(g, dark), edgecolor="none", label=GROUP_LABELS[g]) for g in GROUP_LABELS if any(group_of(d.get("module", ""), modules) == g for _, d in G.nodes(data=True))]
    handles += [
        Line2D([0], [0], color=eff_col["activates"], lw=1.2, label="activates / required for (arrow)"),
        Line2D([0], [0], color=eff_col["inhibits"], lw=1.2, label="inhibits (bar)"),
        Line2D([0], [0], color=eff_col["neutral"], lw=1.2, label="binding / component / modulates"),
        Line2D([0], [0], color=INK2[mode], lw=1.0, ls="-", label="established"),
        Line2D([0], [0], color=INK2[mode], lw=1.0, ls=(0, (4, 2)), label="emerging"),
        Line2D([0], [0], color=INK2[mode], lw=1.0, ls=(0, (1, 2)), label="contested"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=INK2[mode], markersize=6, label="protein / process"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=INK2[mode], markersize=6, label="complex / tissue / organelle"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=INK2[mode], markersize=6, label="metabolite / lipid"),
        Line2D([0], [0], marker="h", color="none", markerfacecolor=INK2[mode], markersize=7, label="drug"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=INK2[mode], markersize=6, label="stimulus"),
        Line2D([0], [0], marker="v", color="none", markerfacecolor=INK2[mode], markersize=6, label="gene variant / RNA"),
        Line2D([0], [0], marker="*", color="none", markerfacecolor=INK2[mode], markersize=8, label="disease / phenotype"),
    ]
    ncol = 4 if figsize[0] >= 16 else 3
    leg = fig.legend(handles=handles, loc="lower center", fontsize=font_size + 0.5, frameon=False, ncol=ncol, labelcolor=INK[mode], bbox_to_anchor=(0.5, 0.0), borderaxespad=0.3)
    leg.set_zorder(5)

    t = title or G.graph.get("name", "mTOR signalling interactome")
    n_pairs = len({(u, v) for u, v in G.edges()})
    st = subtitle or f"{G.number_of_nodes()} nodes · {n_pairs} source–target pairs ({G.number_of_edges()} claim records) · layers top (inputs) → bottom (physiology)"
    ax.set_title(t + "\n" + st, fontsize=font_size + 6, color=INK[mode], loc="left", pad=12)
    ax.margins(x=0.02, y=0.04)
    legend_rows = math.ceil(len(handles) / ncol)
    left = min(0.16, 2.4 / figsize[0])  # room for the layer captions
    fig.subplots_adjust(left=left, right=0.99, top=0.95, bottom=min(0.28, 0.03 + 0.022 * legend_rows * (12 / figsize[1])))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    metadata = None
    if path.suffix.lower() == ".svg":
        metadata = {"Date": None}
    elif path.suffix.lower() == ".pdf":
        metadata = {"CreationDate": None, "ModDate": None}
    with matplotlib.rc_context({"svg.hashsalt": "mtor_graph"}):
        fig.savefig(path, dpi=dpi, facecolor=fig.get_facecolor(), metadata=metadata)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- interactive HTML

_HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<script src="https://unpkg.com/vis-network@9.1.9/standalone/umd/vis-network.min.js"></script>
<style>
:root {
  color-scheme: light;
  --surface-1: #fcfcfb; --page: #f9f9f7; --text-primary: #0b0b0b; --text-secondary: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --border: rgba(11,11,11,0.10);
  --act: #0ca30c; --inh: #d03b3b; --neu: #898781;
  --s1:#2a78d6; --s2:#eb6834; --s3:#1baf7a; --s4:#eda100; --s5:#e87ba4; --s6:#008300; --s7:#4a3aa7; --s8:#e34948;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --surface-1: #1a1a19; --page: #0d0d0d; --text-primary: #ffffff; --text-secondary: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
    --act: #0ca30c; --inh: #e66767; --neu: #898781;
    --s1:#3987e5; --s2:#d95926; --s3:#199e70; --s4:#c98500; --s5:#d55181; --s6:#008300; --s7:#9085e9; --s8:#e66767;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --surface-1: #1a1a19; --page: #0d0d0d; --text-primary: #ffffff; --text-secondary: #c3c2b7; --muted: #898781;
  --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10);
  --act: #0ca30c; --inh: #e66767; --neu: #898781;
  --s1:#3987e5; --s2:#d95926; --s3:#199e70; --s4:#c98500; --s5:#d55181; --s6:#008300; --s7:#9085e9; --s8:#e66767;
}
* { box-sizing: border-box; }
html, body { margin: 0; height: 100%; background: var(--page); color: var(--text-primary); font: 14px/1.4 system-ui, -apple-system, "Segoe UI", sans-serif; }
header { display: flex; flex-wrap: wrap; gap: 8px 16px; align-items: center; padding: 10px 16px; border-bottom: 1px solid var(--grid); background: var(--surface-1); }
header h1 { font-size: 16px; margin: 0 12px 0 0; font-weight: 600; }
header .meta { color: var(--text-secondary); font-size: 12px; }
.controls { display: flex; flex-wrap: wrap; gap: 6px 14px; align-items: center; padding: 8px 16px; background: var(--surface-1); border-bottom: 1px solid var(--grid); font-size: 12px; }
.controls fieldset { border: 0; padding: 0; margin: 0; display: flex; flex-wrap: wrap; gap: 4px 10px; align-items: center; }
.controls legend { float: left; padding: 0; margin-right: 6px; color: var(--text-secondary); font-weight: 600; }
.chip { display: inline-flex; align-items: center; gap: 5px; padding: 2px 8px; border: 1px solid var(--grid); border-radius: 999px; cursor: pointer; user-select: none; background: var(--surface-1); }
.chip input { margin: 0; }
.sw { width: 10px; height: 10px; border-radius: 50%; display: inline-block; }
.line { width: 18px; height: 0; border-top: 2px solid; display: inline-block; vertical-align: middle; }
input[type=search], select { font: inherit; padding: 4px 8px; border: 1px solid var(--axis); border-radius: 6px; background: var(--surface-1); color: var(--text-primary); }
button { font: inherit; font-size: 12px; padding: 4px 10px; border: 1px solid var(--axis); border-radius: 6px; background: var(--surface-1); color: var(--text-primary); cursor: pointer; }
button:hover { border-color: var(--text-secondary); }
main { display: grid; grid-template-columns: 1fr 360px; height: calc(100% - 118px); }
#network { background: var(--surface-1); height: 100%; min-height: 400px; }
aside { border-left: 1px solid var(--grid); overflow: auto; padding: 12px 14px; background: var(--surface-1); font-size: 13px; }
aside h2 { font-size: 15px; margin: 0 0 4px; }
aside h3 { font-size: 12px; margin: 14px 0 4px; color: var(--text-secondary); text-transform: uppercase; letter-spacing: .04em; }
aside .tag { display: inline-block; padding: 1px 6px; border-radius: 4px; border: 1px solid var(--grid); font-size: 11px; margin-right: 4px; color: var(--text-secondary); }
aside .edge { padding: 6px 0; border-top: 1px solid var(--grid); }
aside .edge .who { font-weight: 600; }
aside .edge .mech { color: var(--text-secondary); font-size: 12px; }
aside .edge .refs { font-size: 11px; color: var(--muted); }
aside .edge .refs abbr { text-decoration: none; border-bottom: 1px dotted var(--muted); cursor: help; }
aside a.node { color: var(--s1); cursor: pointer; text-decoration: underline dotted; }
.legend { display: flex; flex-wrap: wrap; gap: 4px 14px; font-size: 12px; color: var(--text-secondary); }
#table { display: none; padding: 12px 16px; overflow: auto; height: 100%; background: var(--surface-1); font-size: 12px; }
#table table { border-collapse: collapse; width: 100%; }
#table th, #table td { text-align: left; padding: 4px 6px; border-bottom: 1px solid var(--grid); vertical-align: top; }
#table th { position: sticky; top: 0; background: var(--surface-1); }
.count { color: var(--text-secondary); }
@media (max-width: 900px) { main { grid-template-columns: 1fr; grid-template-rows: 60vh auto; } aside { border-left: 0; border-top: 1px solid var(--grid); } }
</style>
</head>
<body>
<header>
  <h1>__TITLE__</h1>
  <span class="meta" id="counts"></span>
  <span class="meta">Click a node for details · double-click to focus on its neighbourhood · Esc to clear</span>
</header>
<div class="controls">
  <fieldset id="modules"><legend>Modules</legend></fieldset>
  <fieldset id="evidence"><legend>Evidence</legend>
    <label class="chip"><input type="checkbox" value="established" checked><span class="line" style="border-top-style:solid;border-color:var(--text-secondary)"></span>established</label>
    <label class="chip"><input type="checkbox" value="emerging" checked><span class="line" style="border-top-style:dashed;border-color:var(--text-secondary)"></span>emerging</label>
    <label class="chip"><input type="checkbox" value="contested" checked><span class="line" style="border-top-style:dotted;border-color:var(--text-secondary)"></span>contested</label>
  </fieldset>
  <fieldset id="effects"><legend>Effect</legend>
    <label class="chip"><input type="checkbox" value="activates" checked><span class="line" style="border-color:var(--act)"></span>activates ▶</label>
    <label class="chip"><input type="checkbox" value="inhibits" checked><span class="line" style="border-color:var(--inh)"></span>inhibits ┤</label>
    <label class="chip"><input type="checkbox" value="neutral" checked><span class="line" style="border-color:var(--neu)"></span>binding / component / modulates</label>
  </fieldset>
  <fieldset><legend>Layout</legend>
    <select id="layout"><option value="hier">Hierarchical (tree by layer)</option><option value="force">Force-directed</option></select>
    <select id="depth" title="Focus depth"><option value="1">depth 1</option><option value="2" selected>depth 2</option><option value="3">depth 3</option></select>
  </fieldset>
  <fieldset><legend>Find</legend>
    <input type="search" id="search" placeholder="node, gene, alias…" list="nodelist"><datalist id="nodelist"></datalist>
    <button id="focusBtn">Focus</button><button id="resetBtn">Show all</button><button id="tableBtn">Table view</button><button id="themeBtn">Theme</button>
  </fieldset>
</div>
<main>
  <div id="network"></div>
  <div id="table"></div>
  <aside id="panel">
    <h2>mTOR interaction graph</h2>
    <p class="count">Select a node to see its description, regulators, targets, phosphosites, evidence grade and the publications supporting each interaction.</p>
    <h3>Legend</h3>
    <div class="legend" id="legend"></div>
    <h3>Node shapes</h3>
    <div class="legend">● protein / process &nbsp; ■ complex / tissue &nbsp; ◆ metabolite / lipid &nbsp; ⬢ drug &nbsp; ▲ stimulus &nbsp; ▼ gene variant &nbsp; ★ disease / phenotype</div>
  </aside>
</main>
<script id="kb" type="application/json">__DATA__</script>
<script>
(function(){
  const KB = JSON.parse(document.getElementById('kb').textContent);
  const PAL = { light: ["#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4","#008300","#4a3aa7","#e34948"],
                dark:  ["#3987e5","#d95926","#199e70","#c98500","#d55181","#008300","#9085e9","#e66767"] };
  const EFF = { light: {activates:'#0ca30c', required_for:'#0ca30c', inhibits:'#d03b3b', neutral:'#898781'},
                dark:  {activates:'#0ca30c', required_for:'#0ca30c', inhibits:'#e66767', neutral:'#898781'} };
  const SHAPES = { protein:'dot', complex:'square', metabolite:'diamond', lipid:'diamond', drug:'hexagon', process:'dot', stimulus:'triangle',
                   phenotype:'star', tissue:'square', disease:'star', organelle:'square', rna:'triangleDown', transcription_program:'dot', gene_variant:'triangleDown' };
  const groups = KB.groups; // [{key,label}]
  const groupIndex = Object.fromEntries(groups.map((g,i)=>[g.key,i]));
  const modMeta = Object.fromEntries(KB.modules.map(m=>[m.key,m]));
  const nodesById = Object.fromEntries(KB.nodes.map(n=>[n.id,n]));
  const deg = {}; KB.edges.forEach(e=>{deg[e.from]=(deg[e.from]||0)+1; deg[e.to]=(deg[e.to]||0)+1;});
  const esc = s => String(s==null?'':s).replace(/[&<>"']/g, c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  function mode(){ const t=document.documentElement.getAttribute('data-theme'); if(t) return t; return matchMedia('(prefers-color-scheme: dark)').matches?'dark':'light'; }
  const state = { modules:new Set(KB.modules.map(m=>m.key)), evidence:new Set(['established','emerging','contested']), effects:new Set(['activates','inhibits','neutral']), focus:null, depth:2, layout:'hier', selected:null };

  // controls
  const modFs = document.getElementById('modules');
  KB.modules.forEach(m=>{
    const lab=document.createElement('label'); lab.className='chip';
    lab.innerHTML = `<input type="checkbox" value="${esc(m.key)}" checked><span class="sw" style="background:var(--s${groupIndex[m.group]+1})"></span>${esc(m.title)}`;
    modFs.appendChild(lab);
  });
  const legend=document.getElementById('legend');
  groups.forEach((g,i)=>{ const s=document.createElement('span'); s.innerHTML=`<span class="sw" style="background:var(--s${i+1})"></span> ${esc(g.label)}`; legend.appendChild(s); });
  const dl=document.getElementById('nodelist'); KB.nodes.forEach(n=>{const o=document.createElement('option'); o.value=n.label; dl.appendChild(o);});

  function effClass(e){ return (e.effect==='activates'||e.effect==='required_for')?'activates':(e.effect==='inhibits'?'inhibits':'neutral'); }
  function visibleEdges(){
    return KB.edges.filter(e => state.modules.has(e.module) && state.evidence.has(e.evidence) && state.effects.has(effClass(e)));
  }
  function focusSet(edges){
    if(!state.focus) return null;
    const keep=new Set([state.focus]); let frontier=[state.focus];
    for(let d=0; d<state.depth; d++){ const next=[]; edges.forEach(e=>{ if(frontier.includes(e.from)&&!keep.has(e.to)){keep.add(e.to);next.push(e.to);} if(frontier.includes(e.to)&&!keep.has(e.from)){keep.add(e.from);next.push(e.from);} }); frontier=next; }
    return keep;
  }
  function styleNode(n, m){
    const col = PAL[m][groupIndex[(modMeta[n.module]||{}).group]||0];
    const shape = SHAPES[n.type]||'dot';
    const node = { id:n.id, label:n.label, level:n.layer, shape, size: 8 + Math.min(deg[n.id]||0, 30)*0.5,
      color:{background:col, border:col, highlight:{background:col, border:getComputedStyle(document.documentElement).getPropertyValue('--text-primary').trim()}},
      font:{color:getComputedStyle(document.documentElement).getPropertyValue('--text-primary').trim(), size:12, strokeWidth:0},
      title: n.description ? n.label + ' — ' + n.description : n.label };
    if(n.type==='process'||n.type==='transcription_program'){ node.shapeProperties={borderDashes:[4,2]}; node.borderWidth=2; node.color.background=getComputedStyle(document.documentElement).getPropertyValue('--surface-1').trim(); }
    return node;
  }
  function styleEdge(e, m){
    const cls=effClass(e); const col=EFF[m][cls];
    const ed = { id:e.id, from:e.from, to:e.to, color:{color:col, highlight:col, opacity:0.75}, width:1.2, title:(e.sites&&e.sites.length?e.sites.join(', ')+' · ':'')+e.interaction+' · '+e.mechanism };
    if(cls==='activates') ed.arrows={to:{enabled:true,type:'arrow',scaleFactor:0.6}};
    else if(cls==='inhibits') ed.arrows={to:{enabled:true,type:'bar',scaleFactor:0.8}};
    else if(e.effect==='component_of') ed.arrows={to:{enabled:false}};
    else ed.arrows={to:{enabled:true,type:'circle',scaleFactor:0.5}};
    if(e.evidence==='emerging') ed.dashes=[6,4]; else if(e.evidence==='contested') ed.dashes=[2,4];
    return ed;
  }
  const container=document.getElementById('network');
  let network=null;
  function options(){
    const m=mode();
    const base={ interaction:{hover:true, tooltipDelay:120, navigationButtons:false, keyboard:false, multiselect:false},
      nodes:{borderWidth:1, font:{face:'system-ui, -apple-system, Segoe UI, sans-serif'}},
      edges:{selectionWidth:2, hoverWidth:1.5} };
    if(state.layout==='hier'){
      base.layout={hierarchical:{enabled:true, direction:'UD', sortMethod:'hubsize', levelSeparation:150, nodeSpacing:48, treeSpacing:60, blockShifting:true, edgeMinimization:true, parentCentralization:true}};
      base.physics={enabled:false};
      base.edges.smooth={type:'cubicBezier', forceDirection:'vertical', roundness:0.35};
    } else {
      base.layout={hierarchical:false, improvedLayout:true};
      base.physics={solver:'forceAtlas2Based', forceAtlas2Based:{gravitationalConstant:-45, centralGravity:0.008, springLength:110, springConstant:0.06, avoidOverlap:0.7}, stabilization:{iterations:400, fit:true}};
      base.edges.smooth={type:'continuous', roundness:0.3};
    }
    return base;
  }
  function render(){
    const m=mode();
    let edges=visibleEdges();
    const keep=focusSet(edges);
    if(keep){ edges=edges.filter(e=>keep.has(e.from)&&keep.has(e.to)); }
    const used=new Set(); edges.forEach(e=>{used.add(e.from);used.add(e.to);});
    if(keep) keep.forEach(k=>used.add(k));
    const nodes=KB.nodes.filter(n=>used.has(n.id)).map(n=>styleNode(n,m));
    const vEdges=edges.map(e=>styleEdge(e,m));
    document.getElementById('counts').textContent = `${nodes.length} nodes · ${vEdges.length} interactions shown (of ${KB.nodes.length} / ${KB.edges.length})`;
    const data={nodes:new vis.DataSet(nodes), edges:new vis.DataSet(vEdges)};
    if(network){ network.destroy(); }
    network=new vis.Network(container, data, options());
    network.on('click', p=>{ if(p.nodes.length){ showNode(p.nodes[0]); } else if(p.edges.length){ showEdge(p.edges[0]); } });
    network.on('doubleClick', p=>{ if(p.nodes.length){ state.focus=p.nodes[0]; render(); showNode(p.nodes[0]); } });
    if(state.selected && used.has(state.selected)){ try{ network.selectNodes([state.selected]); }catch(e){} }
    if(!state.focus && !state.selected && used.has('MTORC1') && nodes.length > 120){
      // large graph: open centred on mTORC1 at a readable scale instead of the whole-graph fit
      network.once('afterDrawing', ()=>{ try{ network.focus('MTORC1', {scale:0.75}); }catch(e){} });
    }
    renderTable(edges);
  }
  function refList(e){
    const parts=[];
    (e.refs||[]).forEach(r=>{ const c=KB.references[r]; parts.push(`<abbr title="${esc(c||'')}">[${r}]</abbr>`); });
    return parts.join(' ');
  }
  function edgeHtml(e, dirLabel){
    const other = dirLabel==='in' ? e.from : e.to;
    const glyph = effClass(e)==='activates'?'▶':(effClass(e)==='inhibits'?'┤':'○');
    const who = dirLabel==='in' ? `<a class="node" data-id="${esc(e.from)}">${esc(nodesById[e.from].label)}</a> ${glyph} this` : `this ${glyph} <a class="node" data-id="${esc(e.to)}">${esc(nodesById[e.to].label)}</a>`;
    return `<div class="edge"><div class="who">${who}</div><div class="mech">${esc(e.interaction)}${e.sites&&e.sites.length?' · '+esc(e.sites.join(', ')):''} · <em>${esc(e.evidence)}</em>${e.context?' · '+esc(e.context):''}<br>${esc(e.mechanism)}</div><div class="refs">${refList(e)}</div></div>`;
  }
  function showNode(id){
    state.selected=id; const n=nodesById[id]; if(!n) return;
    const ins=KB.edges.filter(e=>e.to===id), outs=KB.edges.filter(e=>e.from===id);
    const p=document.getElementById('panel');
    p.innerHTML = `<h2>${esc(n.label)}</h2><div><span class="tag">${esc(n.type)}</span><span class="tag">${esc((modMeta[n.module]||{}).title||n.module)}</span><span class="tag">layer ${n.layer}</span>${n.gene?`<span class="tag">gene ${esc(n.gene)}</span>`:''}</div>
      <p>${esc(n.description)}</p>${n.aliases&&n.aliases.length?`<p class="count">Aliases: ${esc(n.aliases.join(', '))}</p>`:''}
      <h3>Regulated by (${ins.length})</h3>${ins.map(e=>edgeHtml(e,'in')).join('')||'<p class="count">—</p>'}
      <h3>Acts on (${outs.length})</h3>${outs.map(e=>edgeHtml(e,'out')).join('')||'<p class="count">—</p>'}`;
    p.querySelectorAll('a.node').forEach(a=>a.addEventListener('click', ()=>{ const t=a.dataset.id; showNode(t); try{network.selectNodes([t]); network.focus(t,{scale:1.2, animation:true});}catch(err){} }));
  }
  function showEdge(id){ const e=KB.edges.find(x=>x.id===id); if(!e) return; document.getElementById('panel').innerHTML=`<h2>${esc(nodesById[e.from].label)} → ${esc(nodesById[e.to].label)}</h2>${edgeHtml(e,'out')}`; }
  function renderTable(edges){
    const t=document.getElementById('table');
    const rows=edges.map(e=>`<tr><td>${esc(nodesById[e.from].label)}</td><td>${esc(e.interaction)}</td><td>${esc(e.effect)}</td><td>${esc(nodesById[e.to].label)}</td><td>${esc((e.sites||[]).join(', '))}</td><td>${esc(e.evidence)}</td><td>${esc((e.refs||[]).map(r=>'['+r+']').join(' '))}</td><td>${esc(e.mechanism)}</td></tr>`).join('');
    t.innerHTML=`<table><thead><tr><th>Source</th><th>Interaction</th><th>Effect</th><th>Target</th><th>Sites</th><th>Evidence</th><th>References</th><th>Mechanism</th></tr></thead><tbody>${rows}</tbody></table>`;
  }
  document.querySelectorAll('#modules input').forEach(i=>i.addEventListener('change', ()=>{ i.checked?state.modules.add(i.value):state.modules.delete(i.value); render(); }));
  document.querySelectorAll('#evidence input').forEach(i=>i.addEventListener('change', ()=>{ i.checked?state.evidence.add(i.value):state.evidence.delete(i.value); render(); }));
  document.querySelectorAll('#effects input').forEach(i=>i.addEventListener('change', ()=>{ i.checked?state.effects.add(i.value):state.effects.delete(i.value); render(); }));
  document.getElementById('layout').addEventListener('change', ev=>{ state.layout=ev.target.value; render(); });
  document.getElementById('depth').addEventListener('change', ev=>{ state.depth=parseInt(ev.target.value,10); if(state.focus) render(); });
  function findNode(q){ q=q.trim().toLowerCase(); if(!q) return null; return KB.nodes.find(n=>n.label.toLowerCase()===q||n.id.toLowerCase()===q) || KB.nodes.find(n=>n.label.toLowerCase().includes(q)||n.id.toLowerCase().includes(q)||(n.gene||'').toLowerCase()===q||(n.aliases||[]).some(a=>a.toLowerCase().includes(q))); }
  document.getElementById('search').addEventListener('change', ev=>{ const n=findNode(ev.target.value); if(n){ showNode(n.id); try{ network.selectNodes([n.id]); network.focus(n.id,{scale:1.2, animation:true}); }catch(e){} } });
  document.getElementById('focusBtn').addEventListener('click', ()=>{ const n=findNode(document.getElementById('search').value)||(state.selected?nodesById[state.selected]:null); if(n){ state.focus=n.id; render(); showNode(n.id);} });
  document.getElementById('resetBtn').addEventListener('click', ()=>{ state.focus=null; render(); });
  document.getElementById('tableBtn').addEventListener('click', ()=>{ const t=document.getElementById('table'); const on=t.style.display==='block'; t.style.display=on?'none':'block'; container.style.display=on?'block':'none'; });
  document.getElementById('themeBtn').addEventListener('click', ()=>{ const cur=mode(); document.documentElement.setAttribute('data-theme', cur==='dark'?'light':'dark'); render(); });
  document.addEventListener('keydown', ev=>{ if(ev.key==='Escape'){ state.focus=null; render(); } });
  matchMedia('(prefers-color-scheme: dark)').addEventListener('change', ()=>{ if(!document.documentElement.getAttribute('data-theme')) render(); });
  render();
})();
</script>
</body>
</html>
"""


def write_html(G: nx.MultiDiGraph, path: str | os.PathLike, title: str = "mTOR signalling interactome") -> Path:
    """Write an HTML explorer; loading vis-network from unpkg requires internet."""
    modules_meta = G.graph.get("modules", {})
    refs = G.graph.get("references", {})
    nodes = []
    for n, d in G.nodes(data=True):
        nodes.append(
            {
                "id": n,
                "label": d.get("label", n),
                "type": d.get("type", "protein"),
                "module": d.get("module", ""),
                "layer": int(d.get("layer", 4)),
                "gene": d.get("gene"),
                "description": d.get("description", ""),
                "aliases": d.get("aliases", []),
            }
        )
    edges = []
    for i, (u, v, d) in enumerate(G.edges(data=True)):
        edges.append(
            {
                "id": f"e{i}",
                "from": u,
                "to": v,
                "interaction": d.get("interaction", "other"),
                "effect": d.get("effect", "modulates"),
                "sites": d.get("sites", []),
                "mechanism": d.get("mechanism", ""),
                "evidence": d.get("evidence", "established"),
                "refs": d.get("refs", []),
                "context": d.get("context"),
                "module": d.get("module", ""),
                "direction": d.get("direction", ""),
            }
        )
    mods = sorted(modules_meta.values(), key=lambda m: m.get("order", 0)) if modules_meta else [
        {"key": k, "title": k, "group": group_of(k, None), "order": i} for i, k in enumerate(sorted({d.get("module", "") for _, d in G.nodes(data=True)}))
    ]
    payload = {
        "nodes": nodes,
        "edges": edges,
        "modules": [{"key": m["key"], "title": m.get("title", m["key"]), "group": m.get("group", "core_complexes")} for m in mods],
        "groups": [{"key": k, "label": v} for k, v in GROUP_LABELS.items()],
        "references": {str(k): v.get("citation", "") for k, v in refs.items()},
    }
    html = (
        _HTML_TEMPLATE.replace("__TITLE__", escape(title))
        .replace("__DATA__", json.dumps(payload, ensure_ascii=False).replace("</", "<\\/"))
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
    return path
