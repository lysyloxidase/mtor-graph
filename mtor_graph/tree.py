"""Hierarchical ("tree") views of the pathway rooted at any node.

The pathway graph is cyclic, so every node is *expanded* at most once by
default. Later occurrences are leaves marked ``↺``. A repeated node can reflect
convergent branches or parallel claims as well as a feedback cycle.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Any, Iterable

import networkx as nx

GLYPH = {"activates": "──▶", "required_for": "──▶", "inhibits": "──┤", "modulates": "──◇", "component_of": "──■", "converts_to": "──»"}
ASCII_GLYPH = {"activates": "-->", "required_for": "-->", "inhibits": "--|", "modulates": "--o", "component_of": "--#", "converts_to": "-->>"}


@dataclass
class TreeNode:
    id: str
    label: str
    depth: int
    effect: str | None = None
    interaction: str | None = None
    sites: list[str] = field(default_factory=list)
    evidence: str | None = None
    refs: list[int] = field(default_factory=list)
    context: str | None = None
    repeated: bool = False
    children: list["TreeNode"] = field(default_factory=list)


def build_tree(
    G: nx.MultiDiGraph,
    root: str = "MTORC1",
    direction: str = "downstream",
    max_depth: int = 4,
    evidence: Iterable[str] | None = None,
    effects: Iterable[str] | None = None,
    expand_once: bool = True,
    exclude_types: Iterable[str] | None = None,
) -> TreeNode:
    """Breadth-first tree from ``root``.

    ``direction`` is ``"downstream"`` (follow out-edges), ``"upstream"`` (follow
    in-edges; children are the regulators of their parent) or ``"both"``.
    Parallel claims remain separate children with their own sites and citations.
    """
    if direction not in {"downstream", "upstream", "both"}:
        raise ValueError("direction must be downstream, upstream or both")
    if max_depth < 0:
        raise ValueError("max_depth must be nonnegative")
    if root not in G:
        raise KeyError(f"unknown node {root!r}")
    ev = set(evidence) if evidence else None
    ef = set(effects) if effects else None
    ex = set(exclude_types) if exclude_types else set()

    def kids(u: str) -> list[tuple[str, dict[str, Any]]]:
        out: list[tuple[str, dict[str, Any]]] = []
        if direction in ("downstream", "both"):
            for v in G.successors(u):
                out.extend((v, d) for d in G[u][v].values())
        if direction in ("upstream", "both"):
            for v in G.predecessors(u):
                if direction == "both" and v == u:
                    continue
                out.extend((v, d) for d in G[v][u].values())
        res = []
        for v, s in out:
            if G.nodes[v].get("type") in ex:
                continue
            if ev and s.get("evidence") not in ev:
                continue
            if ef and s.get("effect") not in ef:
                continue
            res.append((v, s))
        # stable order: signed edges first, then alphabetical
        res.sort(key=lambda x: (0 if x[1]["effect"] in ("activates", "inhibits") else 1, x[0]))
        return res

    tree = TreeNode(id=root, label=G.nodes[root].get("label", root), depth=0)
    expanded = {root}
    q: deque[TreeNode] = deque([tree])
    while q:
        node = q.popleft()
        if node.depth >= max_depth:
            continue
        for v, s in kids(node.id):
            child = TreeNode(
                id=v,
                label=G.nodes[v].get("label", v),
                depth=node.depth + 1,
                effect=s.get("effect"),
                interaction=s.get("interaction"),
                sites=list(s.get("sites", [])),
                evidence=s.get("evidence"),
                refs=list(s.get("refs", [])),
                context=s.get("context"),
            )
            if v in expanded and expand_once:
                child.repeated = True
                node.children.append(child)
                continue
            expanded.add(v)
            node.children.append(child)
            q.append(child)
    return tree


def _fmt(node: TreeNode, show_sites: bool, show_evidence: bool, show_refs: bool, ascii_only: bool) -> str:
    glyphs = ASCII_GLYPH if ascii_only else GLYPH
    parts = []
    if node.effect:
        parts.append(glyphs.get(node.effect, "--?") + " ")
    parts.append(node.label)
    extra = []
    if show_sites and node.sites:
        extra.append(", ".join(node.sites))
    if node.interaction:
        extra.append(node.interaction)
    if node.context:
        extra.append(node.context)
    if show_evidence and node.evidence:
        extra.append(node.evidence)
    if show_refs and node.refs:
        extra.append("[" + ",".join(str(r) for r in node.refs) + "]")
    if extra:
        parts.append("  (" + "; ".join(extra) + ")")
    if node.repeated:
        parts.append("  (repeated)" if ascii_only else "  ↺")
    return "".join(parts)


def render_ascii(
    tree: TreeNode,
    show_sites: bool = True,
    show_evidence: bool = False,
    show_refs: bool = False,
    ascii_only: bool = False,
) -> str:
    """Render the tree with box-drawing characters (or pure ASCII)."""
    lines = [_fmt(tree, show_sites, show_evidence, show_refs, ascii_only)]
    tee, last, pipe, blank = ("|-- ", "`-- ", "|   ", "    ") if ascii_only else ("├── ", "└── ", "│   ", "    ")

    def rec(node: TreeNode, prefix: str) -> None:
        for i, ch in enumerate(node.children):
            is_last = i == len(node.children) - 1
            lines.append(prefix + (last if is_last else tee) + _fmt(ch, show_sites, show_evidence, show_refs, ascii_only))
            rec(ch, prefix + (blank if is_last else pipe))

    rec(tree, "")
    return "\n".join(lines)


def render_markdown(tree: TreeNode, show_sites: bool = True) -> str:
    """Render the tree as a nested Markdown bullet list."""
    lines: list[str] = []

    def rec(node: TreeNode, depth: int) -> None:
        indent = "  " * depth
        if depth == 0:
            lines.append(f"- **{node.label}**")
        else:
            glyph = {"activates": "→", "required_for": "→", "inhibits": "⊣"}.get(node.effect or "", "–")
            s = f" ({', '.join(node.sites)})" if show_sites and node.sites else ""
            if node.context:
                s += f" [{node.context}]"
            loop = " ↺" if node.repeated else ""
            lines.append(f"{indent}- {glyph} {node.label}{s}{loop}")
        for ch in node.children:
            rec(ch, depth + 1)

    rec(tree, 0)
    return "\n".join(lines)


def tree_to_dict(tree: TreeNode) -> dict[str, Any]:
    return {
        "id": tree.id,
        "label": tree.label,
        "depth": tree.depth,
        "effect": tree.effect,
        "interaction": tree.interaction,
        "sites": tree.sites,
        "evidence": tree.evidence,
        "refs": tree.refs,
        "context": tree.context,
        "repeated": tree.repeated,
        "children": [tree_to_dict(c) for c in tree.children],
    }


def count(tree: TreeNode) -> int:
    return 1 + sum(count(c) for c in tree.children)
