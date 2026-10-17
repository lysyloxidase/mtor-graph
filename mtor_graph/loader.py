"""Load the knowledge base from ``data/`` and build a networkx graph."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Iterable

import networkx as nx

from .schema import Edge, Node, EFFECT_SIGN

_PACKAGE_DIR = Path(__file__).resolve().parent
_REPO_DATA_DIR = _PACKAGE_DIR.parent / "data"
DEFAULT_DATA_DIR = _REPO_DATA_DIR if (_REPO_DATA_DIR / "nodes.json").is_file() else _PACKAGE_DIR / "data"


def data_dir(override: str | os.PathLike | None = None) -> Path:
    """Resolve explicit argument > environment > source or installed data."""
    if override:
        return Path(override)
    env = os.environ.get("MTOR_GRAPH_DATA")
    if env:
        return Path(env)
    return DEFAULT_DATA_DIR


def _read_json(path: Path) -> Any:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def load_nodes(directory: str | os.PathLike | None = None) -> list[Node]:
    return [Node.from_dict(d) for d in _read_json(data_dir(directory) / "nodes.json")]


def load_edges(directory: str | os.PathLike | None = None) -> list[Edge]:
    return [Edge.from_dict(d) for d in _read_json(data_dir(directory) / "edges.json")]


def load_references(directory: str | os.PathLike | None = None) -> dict[int, dict[str, Any]]:
    """Bibliography keyed by the citation number used in ``refs``."""
    return {int(r["id"]): r for r in _read_json(data_dir(directory) / "references.json")}


def load_modules(directory: str | os.PathLike | None = None) -> dict[str, dict[str, Any]]:
    """Module metadata (title, summary, colour group) keyed by module key."""
    return {m["key"]: m for m in _read_json(data_dir(directory) / "modules.json")}


def build_graph(
    directory: str | os.PathLike | None = None,
    *,
    nodes: Iterable[Node] | None = None,
    edges: Iterable[Edge] | None = None,
    evidence: Iterable[str] | None = None,
    modules: Iterable[str] | None = None,
) -> nx.MultiDiGraph:
    """Build a :class:`networkx.MultiDiGraph` of the pathway.

    Node attributes mirror :class:`~mtor_graph.schema.Node`; edge attributes mirror
    :class:`~mtor_graph.schema.Edge` plus ``sign`` (+1 / -1 / 0). Parallel edges
    between the same pair (e.g. phosphorylation at different sites, or binding and
    phosphorylation) are kept as separate keyed edges.

    ``evidence`` / ``modules`` optionally restrict the edges loaded (nodes left
    without edges are dropped).
    """
    nodes = list(nodes) if nodes is not None else load_nodes(directory)
    edges = list(edges) if edges is not None else load_edges(directory)
    ev = set(evidence) if evidence else None
    mods = set(modules) if modules else None

    G = nx.MultiDiGraph()
    G.graph["name"] = "mTOR signalling interactome"
    try:
        G.graph["modules"] = load_modules(directory)
    except FileNotFoundError:
        G.graph["modules"] = {}
    try:
        G.graph["references"] = load_references(directory)
    except FileNotFoundError:
        G.graph["references"] = {}

    for n in nodes:
        if n.id in G:
            raise ValueError(f"duplicate node {n.id!r}")
        G.add_node(n.id, **n.to_dict())
    for e in edges:
        if ev and e.evidence not in ev:
            continue
        if mods and e.module not in mods:
            continue
        if e.source not in G or e.target not in G:
            raise ValueError(f"edge {e.source!r}->{e.target!r} has an unknown endpoint")
        if G.has_edge(e.source, e.target, e.key):
            raise ValueError(f"duplicate edge {e.source!r}->{e.target!r} ({e.key})")
        attrs = e.to_dict()
        attrs["sign"] = EFFECT_SIGN.get(e.effect, 0)
        G.add_edge(e.source, e.target, key=e.key, **attrs)

    if ev or mods:
        isolated = [n for n in G.nodes if G.degree(n) == 0]
        G.remove_nodes_from(isolated)
    return G


def to_simple(G: nx.MultiDiGraph) -> nx.DiGraph:
    """Collapse parallel edges into a :class:`networkx.DiGraph`.

    The collapsed edge keeps ``interactions`` (list), ``effects`` (list), a
    summarising ``effect``, the union of citations and the weakest recorded
    ``evidence``. ``sign`` is nonzero only if all signed claims agree; opposing
    signs produce 0 and ``sign_conflict=True``, regardless of their counts.
    Neutral structural edges do not cancel a signed claim. These attributes
    describe the records, not a causal prediction across biological contexts.
    """
    rank = {"established": 0, "emerging": 1, "contested": 2}
    H = nx.DiGraph()
    H.graph.update(G.graph)
    for n, d in G.nodes(data=True):
        H.add_node(n, **d)
    for u, v, d in G.edges(data=True):
        if H.has_edge(u, v):
            h = H[u][v]
            h["interactions"].append(d["interaction"])
            h["effects"].append(d["effect"])
            h["refs"] = sorted(set(h["refs"]) | set(d.get("refs", [])))
            if rank.get(d.get("evidence"), 9) > rank.get(h["evidence"], 9):
                h["evidence"] = d["evidence"]
            h["mechanisms"].append(d.get("mechanism", ""))
        else:
            H.add_edge(
                u,
                v,
                interactions=[d["interaction"]],
                effects=[d["effect"]],
                effect=d["effect"],
                refs=list(d.get("refs", [])),
                evidence=d.get("evidence", "established"),
                mechanisms=[d.get("mechanism", "")],
                module=d.get("module", ""),
                direction=d.get("direction", ""),
            )
    for _, _, d in H.edges(data=True):
        signs = {EFFECT_SIGN.get(effect, 0) for effect in d["effects"]} - {0}
        d["sign_conflict"] = len(signs) > 1
        d["sign"] = next(iter(signs)) if len(signs) == 1 else 0
        if d["sign_conflict"]:
            d["effect"] = "modulates"
        elif d["sign"] == 1:
            d["effect"] = "activates" if "activates" in d["effects"] else "required_for"
        elif d["sign"] == -1:
            d["effect"] = "inhibits"
        else:
            d["effect"] = d["effects"][0] if len(set(d["effects"])) == 1 else "modulates"
    return H
