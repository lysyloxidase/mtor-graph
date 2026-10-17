"""Export the pathway graph to standard network formats."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
from typing import Any

import networkx as nx

from .loader import to_simple
from .viz import group_of, GROUP_LABELS, PALETTE_LIGHT, EFFECT_COLOURS_LIGHT

NODE_SHAPE_DOT = {
    "protein": "ellipse",
    "complex": "box",
    "metabolite": "diamond",
    "lipid": "diamond",
    "drug": "hexagon",
    "process": "box",
    "stimulus": "invtrapezium",
    "phenotype": "octagon",
    "tissue": "folder",
    "disease": "octagon",
    "organelle": "component",
    "rna": "parallelogram",
    "transcription_program": "tab",
    "gene_variant": "invtriangle",
}

ARROWHEAD = {
    "activates": "normal",
    "required_for": "normal",
    "inhibits": "tee",
    "modulates": "odot",
    "component_of": "none",
    "converts_to": "vee",
}


def _scalar(v: Any) -> Any:
    if v is None:
        return ""
    if isinstance(v, (list, tuple, set)):
        return "; ".join(str(x) for x in v)
    if isinstance(v, dict):
        return json.dumps(v, ensure_ascii=False)
    return v


def _scalarised(G: nx.MultiDiGraph) -> nx.MultiDiGraph:
    H = nx.MultiDiGraph()
    H.graph["name"] = G.graph.get("name", "")
    for n, d in G.nodes(data=True):
        H.add_node(n, **{k: _scalar(v) for k, v in d.items()})
    for u, v, k, d in G.edges(keys=True, data=True):
        H.add_edge(u, v, key=k, **{kk: _scalar(vv) for kk, vv in d.items()})
    return H


def write_graphml(G: nx.MultiDiGraph, path: str | os.PathLike) -> Path:
    path = Path(path)
    nx.write_graphml(_scalarised(G), path, encoding="utf-8")
    return path


def write_gexf(G: nx.MultiDiGraph, path: str | os.PathLike) -> Path:
    path = Path(path)
    writer = nx.readwrite.gexf.GEXFWriter(_scalarised(G), encoding="utf-8")
    # NetworkX otherwise inserts the build date, making identical data differ.
    meta = writer.xml.find("meta")
    if meta is not None:
        meta.attrib.pop("lastmodifieddate", None)
    writer.write(path)
    return path


def write_cytoscape_json(G: nx.MultiDiGraph, path: str | os.PathLike) -> Path:
    """Cytoscape.js JSON (``elements`` with ``nodes`` and ``edges``)."""
    path = Path(path)
    H = nx.MultiDiGraph()
    for n, d in G.nodes(data=True):
        H.add_node(n, **d)
    for u, v, k, d in G.edges(keys=True, data=True):
        H.add_edge(u, v, key=k, **d)
    data = nx.cytoscape_data(H)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=1)
    return path


def write_sif(G: nx.MultiDiGraph, path: str | os.PathLike) -> Path:
    """Simple Interaction Format: ``source<TAB>interaction<TAB>target``."""
    path = Path(path)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for u, v, d in G.edges(data=True):
            fh.write(f"{u}\t{d.get('interaction', 'other')}:{d.get('effect', '')}\t{v}\n")
    return path


def write_csv(G: nx.MultiDiGraph, nodes_path: str | os.PathLike, edges_path: str | os.PathLike) -> tuple[Path, Path]:
    nodes_path, edges_path = Path(nodes_path), Path(edges_path)
    node_fields = ["id", "label", "type", "gene", "module", "layer", "aliases", "description"]
    with open(nodes_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=node_fields, extrasaction="ignore")
        w.writeheader()
        for n, d in G.nodes(data=True):
            w.writerow({k: _scalar(d.get(k)) for k in node_fields})
    edge_fields = [
        "source", "target", "interaction", "effect", "sites", "evidence", "refs",
        "mechanism", "context", "direction", "module", "notes",
    ]
    with open(edges_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=edge_fields, extrasaction="ignore")
        w.writeheader()
        for u, v, d in G.edges(data=True):
            row = {k: _scalar(d.get(k)) for k in edge_fields}
            row["source"], row["target"] = u, v
            w.writerow(row)
    return nodes_path, edges_path


def write_json(G: nx.MultiDiGraph, path: str | os.PathLike) -> Path:
    path = Path(path)
    data = {
        "name": G.graph.get("name"),
        "nodes": [d for _, d in G.nodes(data=True)],
        "edges": [dict(d, source=u, target=v) for u, v, d in G.edges(data=True)],
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, ensure_ascii=False, indent=1)
    return path


def _dot_escape(s: str) -> str:
    return str(s).replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")


def write_dot(
    G: nx.MultiDiGraph,
    path: str | os.PathLike,
    *,
    rank_by_layer: bool = True,
    edge_labels: bool = True,
    collapse: bool = False,
    rankdir: str = "TB",
) -> Path:
    """Graphviz DOT file (render with ``dot -Tsvg file.dot -o file.svg``).

    Nodes are coloured by module group and shaped by node type; edges carry an
    arrowhead by effect (``normal`` activation, ``tee`` inhibition, ``odot``
    modulation) and a line style by evidence (solid / dashed / dotted). With
    ``rank_by_layer`` nodes of the same layer share a rank, which produces the
    top-down "tree" view of the pathway.
    """
    path = Path(path)
    modules = G.graph.get("modules", {})
    lines = [
        "digraph mTOR {",
        f'  graph [rankdir={rankdir}, splines=true, overlap=false, nodesep=0.25, ranksep=0.9, fontname="Helvetica", label="{_dot_escape(G.graph.get("name", "mTOR signalling interactome"))}", labelloc=t, fontsize=18];',
        '  node [style="filled", fontname="Helvetica", fontsize=9, penwidth=0.6, color="#898781"];',
        '  edge [fontname="Helvetica", fontsize=7, penwidth=0.9, arrowsize=0.7, color="#898781"];',
    ]
    group_colour = {g: PALETTE_LIGHT[i] for i, g in enumerate(GROUP_LABELS)}
    light_fill = {
        "nutrient_sensing": "#cde2fb", "upstream_inputs": "#fde1d6", "core_complexes": "#d2f3e5",
        "mtorc1_outputs": "#fff0c2", "mtorc2": "#fbe0eb", "circadian": "#d6ead6",
        "tissue_disease": "#e0dcf5", "pharmacology": "#fbd9d8",
    }
    for n, d in G.nodes(data=True):
        g = group_of(d.get("module", ""), modules)
        shape = NODE_SHAPE_DOT.get(d.get("type", "protein"), "ellipse")
        style = "filled,rounded" if shape == "box" else "filled"
        tooltip = _dot_escape((d.get("description") or d.get("label", n))[:300])
        lines.append(
            f'  "{n}" [label="{_dot_escape(d.get("label", n))}", shape={shape}, style="{style}", '
            f'fillcolor="{light_fill.get(g, "#eeeeee")}", color="{group_colour.get(g, "#898781")}", tooltip="{tooltip}"];'
        )
    if rank_by_layer:
        layers: dict[int, list[str]] = {}
        for n, d in G.nodes(data=True):
            layers.setdefault(int(d.get("layer", 4)), []).append(n)
        for l in sorted(layers):
            members = " ".join(f'"{n}"' for n in layers[l])
            lines.append(f"  {{ rank=same; {members} }}")
    style_by_evidence = {"established": "solid", "emerging": "dashed", "contested": "dotted"}
    edge_iter = ((u, v, d) for u, v, d in (to_simple(G).edges(data=True) if collapse else G.edges(data=True)))
    for u, v, d in edge_iter:
        effect = d.get("effect", "modulates")
        colour = EFFECT_COLOURS_LIGHT.get(effect, EFFECT_COLOURS_LIGHT["neutral"])
        attrs = [f'arrowhead={ARROWHEAD.get(effect, "normal")}', f'color="{colour}"', f'style={style_by_evidence.get(d.get("evidence"), "solid")}']
        if effect == "component_of":
            attrs.append("arrowhead=none")
        lab = ""
        if edge_labels:
            sites = d.get("sites") or []
            if isinstance(sites, str):
                sites = [sites]
            lab = ", ".join(sites)
        if lab:
            attrs.append(f'label="{_dot_escape(lab)}"')
        tip = _dot_escape((d.get("mechanism") or "")[:300]) if not collapse else _dot_escape("; ".join(d.get("mechanisms", []))[:300])
        attrs.append(f'tooltip="{tip}"')
        lines.append(f'  "{u}" -> "{v}" [{", ".join(attrs)}];')
    lines.append("}")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def export_all(G: nx.MultiDiGraph, outdir: str | os.PathLike, stem: str = "mtor_interactome") -> dict[str, Path]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out: dict[str, Path] = {}
    out["graphml"] = write_graphml(G, outdir / f"{stem}.graphml")
    out["gexf"] = write_gexf(G, outdir / f"{stem}.gexf")
    out["cyjs"] = write_cytoscape_json(G, outdir / f"{stem}.cyjs")
    out["sif"] = write_sif(G, outdir / f"{stem}.sif")
    n, e = write_csv(G, outdir / f"{stem}_nodes.csv", outdir / f"{stem}_edges.csv")
    out["nodes_csv"], out["edges_csv"] = n, e
    out["json"] = write_json(G, outdir / f"{stem}.json")
    out["dot"] = write_dot(G, outdir / f"{stem}.dot")
    return out
