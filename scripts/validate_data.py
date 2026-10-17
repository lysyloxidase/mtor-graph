"""Validate data/nodes.json, data/edges.json and data/references.json; exit 1 on problems."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mtor_graph import load_edges, load_nodes, load_references, validate  # noqa: E402


def main() -> int:
    nodes, edges, refs = load_nodes(), load_edges(), load_references()
    problems = validate(nodes, edges, refs)
    for p in problems:
        print(p)
    print(f"{len(nodes)} nodes, {len(edges)} edges, {len(refs)} references, {len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
