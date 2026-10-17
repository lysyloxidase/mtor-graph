"""Regenerate every derived artefact in ``output/`` from ``data/``.

Run from the repository root::

    python scripts/build_all.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mtor_graph import build_graph, build_tree, render_ascii, render_markdown, by_module, neighborhood, summary, load_modules  # noqa: E402
from mtor_graph import export, viz  # noqa: E402
from mtor_graph.tree import count  # noqa: E402


def main() -> None:
    out = ROOT / "output"
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "trees").mkdir(parents=True, exist_ok=True)
    (out / "exports").mkdir(parents=True, exist_ok=True)
    G = build_graph()
    print(f"graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # exports
    for k, p in export.export_all(G, out / "exports").items():
        print(f"  export {k:10s} -> {p.relative_to(ROOT)}")

    # interactive explorer
    p = viz.write_html(G, out / "mtor_interactome.html")
    print(f"  html -> {p.relative_to(ROOT)}")

    # static figures
    p = viz.draw_layered(G, out / "figures" / "mtor_interactome_full.png", dpi=150)
    print(f"  figure -> {p.relative_to(ROOT)}")
    viz.draw_layered(G, out / "figures" / "mtor_interactome_full.svg")
    viz.draw_layered(G, out / "figures" / "mtor_interactome_full_dark.png", dark=True, dpi=150)
    for key, meta in load_modules().items():
        sub = by_module(G, key)
        if sub.number_of_edges() == 0:
            continue
        p = viz.draw_layered(sub, out / "figures" / f"module_{key}.png", title=meta["title"], dpi=150, show_sites=False)
        print(f"  figure -> {p.relative_to(ROOT)}")
    for core in ("MTORC1", "MTORC2"):
        sub = neighborhood(G, core, depth=1)
        viz.draw_layered(sub, out / "figures" / f"ego_{core}.png", title=f"Direct interactors of {core}", dpi=150, collapse=False, show_sites=True)

    # trees
    specs = [
        ("MTORC1", "downstream", 4), ("MTORC1", "upstream", 4), ("MTORC2", "downstream", 4), ("MTORC2", "upstream", 3),
        ("MTOR", "both", 2), ("LEUCINE", "downstream", 5), ("INSULIN", "downstream", 6), ("RAPAMYCIN", "downstream", 4), ("AMPK", "downstream", 4),
    ]
    index_lines = ["# Tree views of the mTOR pathway", ""]
    for root, direction, depth in specs:
        if root not in G:
            continue
        t = build_tree(G, root, direction=direction, max_depth=depth, exclude_types=None)
        stem = f"tree_{root}_{direction}_d{depth}"
        (out / "trees" / f"{stem}.txt").write_text(render_ascii(t, show_refs=True), encoding="utf-8")
        (out / "trees" / f"{stem}.md").write_text(f"# {root} — {direction}, depth {depth}\n\n" + render_markdown(t), encoding="utf-8")
        index_lines.append(f"- [{root} {direction} (depth {depth})]({stem}.md) — {count(t)} nodes")
        print(f"  tree {stem}: {count(t)} nodes")
    (out / "trees" / "README.md").write_text("\n".join(index_lines) + "\n", encoding="utf-8")

    # stats
    s = summary(G)
    s["source_target_pairs"] = len({(u, v) for u, v in G.edges()})
    s["distinct_interactions"] = len({(u, v, d.get("interaction"), d.get("effect")) for u, v, d in G.edges(data=True)})
    (out / "stats.json").write_text(json.dumps(s, indent=1, ensure_ascii=False), encoding="utf-8")
    lines = [f"# Knowledge-base statistics", "", f"- nodes: {s['nodes']}",
             f"- claim records (edges in the multigraph): {s['edges']}",
             f"- distinct source–target–interaction–effect combinations: {s['distinct_interactions']}",
             f"- distinct source–target pairs: {s['source_target_pairs']}",
             f"- cited claim records: {s['cited_edges']}",
             f"- distinct references cited: {s['distinct_references_cited']}", ""]
    for key in ("node_types", "edge_modules", "interactions", "effects", "evidence"):
        lines.append(f"## {key}")
        lines.append("")
        lines.append("| value | count |")
        lines.append("|---|---|")
        for k, v in s[key].items():
            lines.append(f"| {k} | {v} |")
        lines.append("")
    lines.append("## Hubs (total degree)")
    lines.append("")
    lines.append("| node | in | out |")
    lines.append("|---|---|---|")
    for h in s["hubs"]:
        lines.append(f"| {h[0]} | {h[1]} | {h[2]} |")
    (out / "stats.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # refresh the summary block in README.md between the STATS markers
    readme = ROOT / "README.md"
    if readme.exists():
        txt = readme.read_text(encoding="utf-8")
        mods = load_modules()
        block = [
            "<!--STATS-->",
            f"**{s['nodes']} nodes · {s['source_target_pairs']} directed source–target pairs · {s['edges']} claim records** "
            f"(a pair can carry several records when they differ in interaction type, sign, sites, context or evidence). "
            f"Every record cites the bibliography ({s['distinct_references_cited']} distinct publications); "
            f"recorded evidence grades: {s['evidence'].get('established', 0)} established, {s['evidence'].get('emerging', 0)} emerging, "
            f"{s['evidence'].get('contested', 0)} contested. Node types: "
            + ", ".join(f"{v} {k.replace('_', ' ')}" for k, v in s["node_types"].items())
            + ". Claim records per module: "
            + ", ".join(f"{mods.get(k, {}).get('title', k)} {v}" for k, v in s["edge_modules"].items())
            + ".",
            "<!--/STATS-->",
        ]
        import re

        if "<!--/STATS-->" in txt:
            txt = re.sub(r"<!--STATS-->.*?<!--/STATS-->", "\n".join(block), txt, flags=re.S)
        else:
            txt = txt.replace("<!--STATS-->", "\n".join(block), 1)
        readme.write_text(txt, encoding="utf-8")
        print("  README stats block refreshed")
    print("done")


if __name__ == "__main__":
    main()
