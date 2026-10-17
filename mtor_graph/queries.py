"""Query helpers over the pathway graph."""

from __future__ import annotations

from collections import Counter, deque
from itertools import islice
from typing import Any, Iterable

import networkx as nx

from .loader import to_simple


def _bfs(G: nx.MultiDiGraph, start: str, depth: int | None, reverse: bool) -> dict[str, int]:
    if depth is not None and depth < 0:
        raise ValueError("depth must be nonnegative")
    if start not in G:
        raise KeyError(f"unknown node {start!r}")
    nxt = G.predecessors if reverse else G.successors
    dist = {start: 0}
    q = deque([start])
    while q:
        u = q.popleft()
        if depth is not None and dist[u] >= depth:
            continue
        for v in nxt(u):
            if v not in dist:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist


def upstream(G: nx.MultiDiGraph, node: str, depth: int | None = None) -> dict[str, int]:
    """Nodes that can reach ``node`` (within ``depth`` steps) → distance."""
    d = _bfs(G, node, depth, reverse=True)
    d.pop(node, None)
    return d


def downstream(G: nx.MultiDiGraph, node: str, depth: int | None = None) -> dict[str, int]:
    """Nodes reachable from ``node`` (within ``depth`` steps) → distance."""
    d = _bfs(G, node, depth, reverse=False)
    d.pop(node, None)
    return d


def neighborhood(G: nx.MultiDiGraph, node: str, depth: int = 1) -> nx.MultiDiGraph:
    """Ego sub-graph containing everything within ``depth`` steps in either direction."""
    if depth < 0:
        raise ValueError("depth must be nonnegative")
    if node not in G:
        raise KeyError(f"unknown node {node!r}")
    keep = nx.single_source_shortest_path_length(G.to_undirected(as_view=True), node, cutoff=depth)
    return G.subgraph(keep).copy()


def shortest_paths(G: nx.MultiDiGraph, source: str, target: str, k: int = 5) -> list[list[str]]:
    """Up to ``k`` shortest simple directed paths, ranked by hop count.

    All interaction types are traversable. Signs, evidence and compatibility
    of contexts do not affect ranking, so a path is not a causal prediction.
    """
    if k < 0:
        raise ValueError("k must be nonnegative")
    H = nx.DiGraph(G)
    try:
        return list(islice(nx.shortest_simple_paths(H, source, target), k))
    except nx.NetworkXNoPath:
        return []


def path_sign(G: nx.MultiDiGraph, path: list[str]) -> int:
    """Product of recorded signs; 0 means neutral, unknown or conflicting.

    This is only a descriptive calculation. Context compatibility, causal
    sufficiency, feedback strength and kinetics are not modelled.
    """
    if not path:
        raise ValueError("path must contain at least one node")
    for node in path:
        if node not in G:
            raise KeyError(f"unknown node {node!r}")
    H = to_simple(G)
    s = 1
    for u, v in zip(path, path[1:]):
        s *= H[u][v].get("sign", 0)
    return s


def by_module(G: nx.MultiDiGraph, *modules: str, include_node_module: bool = True) -> nx.MultiDiGraph:
    """Sub-graph of edges whose ``module`` is in ``modules`` (plus the nodes they touch)."""
    mods = set(modules)
    H = nx.MultiDiGraph()
    H.graph.update(G.graph)
    for u, v, k, d in G.edges(keys=True, data=True):
        if d.get("module") in mods or (
            include_node_module and (G.nodes[u].get("module") in mods and G.nodes[v].get("module") in mods)
        ):
            for n in (u, v):
                if n not in H:
                    H.add_node(n, **G.nodes[n])
            H.add_edge(u, v, key=k, **d)
    return H


def filter_graph(
    G: nx.MultiDiGraph,
    *,
    evidence: Iterable[str] | None = None,
    effects: Iterable[str] | None = None,
    interactions: Iterable[str] | None = None,
    directions: Iterable[str] | None = None,
    node_types: Iterable[str] | None = None,
    min_refs: int = 0,
    drop_isolated: bool = True,
) -> nx.MultiDiGraph:
    """Return a copy keeping only edges matching all given criteria."""
    ev, ef, it, di = (set(x) if x else None for x in (evidence, effects, interactions, directions))
    nt = set(node_types) if node_types else None
    H = nx.MultiDiGraph()
    H.graph.update(G.graph)
    for n, d in G.nodes(data=True):
        if nt is None or d.get("type") in nt:
            H.add_node(n, **d)
    for u, v, k, d in G.edges(keys=True, data=True):
        if u not in H or v not in H:
            continue
        if ev and d.get("evidence") not in ev:
            continue
        if ef and d.get("effect") not in ef:
            continue
        if it and d.get("interaction") not in it:
            continue
        if di and d.get("direction") not in di:
            continue
        if len(d.get("refs", [])) < min_refs:
            continue
        H.add_edge(u, v, key=k, **d)
    if drop_isolated:
        H.remove_nodes_from([n for n in H.nodes if H.degree(n) == 0])
    return H


def search(G: nx.MultiDiGraph, term: str) -> list[str]:
    """Node ids whose id, label, gene, aliases or description contain ``term`` (case-insensitive)."""
    t = term.lower()
    hits: list[tuple[int, str]] = []
    for n, d in G.nodes(data=True):
        hay = [n, d.get("label", ""), d.get("gene") or "", *d.get("aliases", [])]
        if any(t == h.lower() for h in hay):
            hits.append((0, n))
        elif any(t in h.lower() for h in hay):
            hits.append((1, n))
        elif t in (d.get("description") or "").lower():
            hits.append((2, n))
    return [n for _, n in sorted(hits)]


def edges_between(G: nx.MultiDiGraph, a: str, b: str, both_directions: bool = True) -> list[dict[str, Any]]:
    out = []
    pairs = [(a, b), (b, a)] if both_directions and a != b else [(a, b)]
    for u, v in pairs:
        if G.has_edge(u, v):
            for k, d in G[u][v].items():
                out.append({"source": u, "target": v, "key": k, **d})
    return out


def hubs(G: nx.MultiDiGraph, n: int = 15) -> list[tuple[str, int, int]]:
    """Top-``n`` nodes by total degree → (id, in_degree, out_degree)."""
    deg = sorted(G.nodes, key=lambda x: G.in_degree(x) + G.out_degree(x), reverse=True)
    return [(x, G.in_degree(x), G.out_degree(x)) for x in deg[:n]]


def summary(G: nx.MultiDiGraph) -> dict[str, Any]:
    """Descriptive statistics of the knowledge base."""
    edges = [d for _, _, d in G.edges(data=True)]
    nodes = [d for _, d in G.nodes(data=True)]
    cited = Counter(r for d in edges for r in d.get("refs", []))
    return {
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "node_types": dict(Counter(d.get("type") for d in nodes).most_common()),
        "node_modules": dict(Counter(d.get("module") for d in nodes).most_common()),
        "edge_modules": dict(Counter(d.get("module") for d in edges).most_common()),
        "interactions": dict(Counter(d.get("interaction") for d in edges).most_common()),
        "effects": dict(Counter(d.get("effect") for d in edges).most_common()),
        "evidence": dict(Counter(d.get("evidence") for d in edges).most_common()),
        "directions": dict(Counter(d.get("direction") for d in edges).most_common()),
        "cited_edges": sum(1 for d in edges if d.get("refs")),
        "distinct_references_cited": len(cited),
        "layers": dict(sorted(Counter(d.get("layer") for d in nodes).items())),
        "hubs": hubs(G, 10),
    }
