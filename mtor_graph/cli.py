"""Command-line interface: ``mtor-graph <command>`` (or ``python -m mtor_graph.cli``)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import loader, queries, tree as tree_mod
from .schema import EVIDENCE_LEVELS, MODULE_KEYS, validate


def _csv(allowed):
    """A comma-separated list option, so a value never swallows the subcommand name."""
    allowed = sorted(allowed)

    def parse(raw: str) -> list[str]:
        values = [v.strip() for v in raw.split(",") if v.strip()]
        bad = [v for v in values if v not in allowed]
        if bad or not values:
            raise argparse.ArgumentTypeError(
                f"{', '.join(bad) or 'no value'} — choose from {', '.join(allowed)} (comma-separated)")
        return values

    return parse


def _graph(args):
    return loader.build_graph(args.data, evidence=getattr(args, "evidence", None), modules=getattr(args, "modules", None))


def cmd_stats(args):
    G = _graph(args)
    s = queries.summary(G)
    if args.json:
        print(json.dumps(s, indent=1, ensure_ascii=False))
        return
    pairs = len({(u, v) for u, v in G.edges()})
    print(f"{s['nodes']} nodes, {pairs} source-target pairs, {s['edges']} claim records")
    for key in ("node_types", "edge_modules", "interactions", "effects", "evidence", "directions"):
        print(f"\n{key}:")
        for k, v in s[key].items():
            print(f"  {k:45s} {v}")
    print(f"\ncited records: {s['cited_edges']}; distinct references cited: {s['distinct_references_cited']}")
    print("\nhubs (id, in, out):")
    for h in s["hubs"]:
        print(f"  {h[0]:30s} {h[1]:3d} {h[2]:3d}")


def cmd_tree(args):
    G = _graph(args)
    t = tree_mod.build_tree(
        G, root=args.root, direction=args.direction, max_depth=args.depth,
        effects=args.effects, exclude_types=args.exclude_types,
    )
    if args.json:
        print(json.dumps(tree_mod.tree_to_dict(t), indent=1, ensure_ascii=False))
    elif args.markdown:
        print(tree_mod.render_markdown(t))
    else:
        print(tree_mod.render_ascii(t, show_sites=not args.no_sites, show_evidence=args.show_evidence, show_refs=args.refs, ascii_only=args.ascii))
    print(f"\n[{tree_mod.count(t)} tree nodes; root {args.root}; {args.direction}; depth {args.depth}]", file=sys.stderr)


def _print_edge(G, u, v, d, refs):
    sites = f" @ {', '.join(d['sites'])}" if d.get("sites") else ""
    ctx = f" [{d['context']}]" if d.get("context") else ""
    print(f"  {u} --{d['interaction']}/{d['effect']}{sites}--> {v}  ({d['evidence']}){ctx}")
    print(f"      {d.get('mechanism', '')}")
    cites = [f"[{r}]" for r in d.get("refs", [])]
    if cites:
        print(f"      refs: {' '.join(cites)}")
    if refs:
        for r in d.get("refs", []):
            c = G.graph.get("references", {}).get(r, {}).get("citation", "?")
            print(f"        [{r}] {c}")


def cmd_node(args):
    G = _graph(args)
    ids = [args.id] if args.id in G else queries.search(G, args.id)
    if not ids:
        sys.exit(f"no node matches {args.id!r}")
    for n in ids[: args.limit]:
        d = G.nodes[n]
        print(f"{n}: {d.get('label')}  [{d.get('type')}; module {d.get('module')}; layer {d.get('layer')}; gene {d.get('gene')}]")
        print(f"  {d.get('description')}")
        if d.get("aliases"):
            print(f"  aliases: {', '.join(d['aliases'])}")
        print(f"  regulated by ({G.in_degree(n)}):")
        for u, _, dd in G.in_edges(n, data=True):
            _print_edge(G, u, n, dd, args.refs)
        print(f"  acts on ({G.out_degree(n)}):")
        for _, v, dd in G.out_edges(n, data=True):
            _print_edge(G, n, v, dd, args.refs)
        print()


def cmd_search(args):
    G = _graph(args)
    for n in queries.search(G, args.term):
        d = G.nodes[n]
        print(f"{n:32s} {d.get('label'):40s} {d.get('type'):12s} {d.get('module')}")


def cmd_path(args):
    G = _graph(args)
    paths = queries.shortest_paths(G, args.source, args.target, k=args.k)
    if not paths:
        print("no directed path")
        return
    H = loader.to_simple(G)
    for p in paths:
        parts = [p[0]]
        for u, v in zip(p, p[1:]):
            e = H[u][v]
            glyph = {1: "─▶", -1: "─┤"}.get(e.get("sign", 0), "─○")
            parts.append(f" {glyph} {v}")
        net = queries.path_sign(G, p)
        label = '+' if net > 0 else '-' if net < 0 else 'undetermined'
        print("".join(parts) + f"    (recorded sign product: {label})")
    print("Paths show graph connectivity; contexts and causal effects are not modelled.")


def cmd_module(args):
    G = _graph(args)
    H = queries.by_module(G, args.key)
    meta = G.graph.get("modules", {}).get(args.key, {})
    print(f"{args.key}: {meta.get('title', '')}")
    if meta.get("summary"):
        print(f"  {meta['summary']}")
    print(f"{H.number_of_nodes()} nodes, {H.number_of_edges()} interactions\n")
    for u, v, d in H.edges(data=True):
        _print_edge(G, u, v, d, args.refs)


def cmd_export(args):
    from . import export

    G = _graph(args)
    out = export.export_all(G, args.outdir, stem=args.stem)
    for k, p in out.items():
        print(f"{k:10s} {p}")


def cmd_draw(args):
    from . import viz

    G = _graph(args)
    if args.module:
        G = queries.by_module(G, *args.module)
    if args.ego:
        G = queries.neighborhood(G, args.ego, args.depth)
    p = viz.draw_layered(G, args.out, dark=args.dark, collapse=not args.no_collapse, title=args.title, show_sites=args.sites)
    print(p)


def cmd_html(args):
    from . import viz

    G = _graph(args)
    if args.module:
        G = queries.by_module(G, *args.module)
    print(viz.write_html(G, args.out, title=args.title))


def cmd_validate(args):
    nodes = loader.load_nodes(args.data)
    edges = loader.load_edges(args.data)
    refs = loader.load_references(args.data)
    problems = validate(nodes, edges, refs)
    if problems:
        print("\n".join(problems))
        sys.exit(f"{len(problems)} problem(s)")
    print(f"OK: {len(nodes)} nodes, {len(edges)} edges, {len(refs)} references")


def cmd_refs(args):
    refs = loader.load_references(args.data)
    for n in args.numbers:
        r = refs.get(n)
        print(f"[{n}] {r['citation'] if r else '— not in list —'}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="mtor-graph", description="Query, visualise and export the mTOR pathway interaction graph.")
    p.add_argument("--data", help="data directory (default: source/installed data or $MTOR_GRAPH_DATA)")
    p.add_argument("--evidence", type=_csv(EVIDENCE_LEVELS),
                   help="comma-separated evidence labels to keep, e.g. --evidence established,emerging")
    p.add_argument("--modules", type=_csv(MODULE_KEYS),
                   help="comma-separated module keys to keep, e.g. --modules circadian,pharmacology")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("stats", help="descriptive statistics"); s.add_argument("--json", action="store_true"); s.set_defaults(fn=cmd_stats)

    s = sub.add_parser("tree", help="hierarchical tree from a root node")
    s.add_argument("--root", default="MTORC1"); s.add_argument("--direction", choices=["downstream", "upstream", "both"], default="downstream")
    s.add_argument("--depth", type=int, default=3); s.add_argument("--effects", nargs="+"); s.add_argument("--exclude-types", nargs="+")
    s.add_argument("--no-sites", action="store_true"); s.add_argument("--show-evidence", action="store_true"); s.add_argument("--refs", action="store_true")
    s.add_argument("--ascii", action="store_true", help="pure-ASCII glyphs"); s.add_argument("--markdown", action="store_true"); s.add_argument("--json", action="store_true")
    s.set_defaults(fn=cmd_tree)

    s = sub.add_parser("node", help="show a node with all its interactions"); s.add_argument("id"); s.add_argument("--refs", action="store_true", help="print full citations"); s.add_argument("--limit", type=int, default=5); s.set_defaults(fn=cmd_node)
    s = sub.add_parser("search", help="find nodes by id / label / gene / alias / description"); s.add_argument("term"); s.set_defaults(fn=cmd_search)
    s = sub.add_parser("path", help="shortest directed paths between two nodes"); s.add_argument("source"); s.add_argument("target"); s.add_argument("-k", type=int, default=5); s.set_defaults(fn=cmd_path)
    s = sub.add_parser("module", help="list the interactions of one module"); s.add_argument("key"); s.add_argument("--refs", action="store_true"); s.set_defaults(fn=cmd_module)
    s = sub.add_parser("export", help="write GraphML / GEXF / Cytoscape JSON / DOT / SIF / CSV / JSON"); s.add_argument("outdir"); s.add_argument("--stem", default="mtor_interactome"); s.set_defaults(fn=cmd_export)
    s = sub.add_parser("draw", help="static layered figure (png / svg / pdf)"); s.add_argument("out"); s.add_argument("--module", nargs="+"); s.add_argument("--ego"); s.add_argument("--depth", type=int, default=1)
    s.add_argument("--dark", action="store_true"); s.add_argument("--no-collapse", action="store_true"); s.add_argument("--sites", action="store_true"); s.add_argument("--title"); s.set_defaults(fn=cmd_draw)
    s = sub.add_parser("html", help="interactive HTML explorer"); s.add_argument("out"); s.add_argument("--module", nargs="+"); s.add_argument("--title", default="mTOR signalling interactome"); s.set_defaults(fn=cmd_html)
    s = sub.add_parser("validate", help="validate the data files"); s.set_defaults(fn=cmd_validate)
    s = sub.add_parser("refs", help="print bibliography entries by number"); s.add_argument("numbers", nargs="+", type=int); s.set_defaults(fn=cmd_refs)
    return p


def main(argv: list[str] | None = None) -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8")  # type: ignore[attr-defined]
    except Exception:
        pass
    args = build_parser().parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
