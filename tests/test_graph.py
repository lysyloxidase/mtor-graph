import json

import pytest

from mtor_graph import (
    build_graph, to_simple, upstream, downstream, neighborhood, shortest_paths, by_module,
    filter_graph, search, summary, hubs, build_tree, render_ascii, render_markdown, tree_to_dict,
)
from mtor_graph import export, viz
from mtor_graph.tree import count


@pytest.fixture(scope="module")
def G():
    return build_graph()


def test_build_graph_basic(G):
    assert G.number_of_nodes() > 150
    assert G.number_of_edges() > 400
    assert G.graph["modules"]
    assert G.graph["references"]


def test_to_simple_collapses(G):
    H = to_simple(G)
    assert H.number_of_edges() <= G.number_of_edges()
    assert all("interactions" in d for _, _, d in H.edges(data=True))


def test_upstream_downstream(G):
    up = upstream(G, "MTORC1", depth=2)
    down = downstream(G, "MTORC1", depth=2)
    assert "RHEB" in up and "RPS6KB1" in down
    assert all(v <= 2 for v in up.values())


def test_neighborhood(G):
    H = neighborhood(G, "RHEB", depth=1)
    assert "RHEB" in H and "MTORC1" in H and "TSC2" in H


def test_shortest_paths(G):
    paths = shortest_paths(G, "INSULIN", "RPS6", k=3)
    assert paths and paths[0][0] == "INSULIN" and paths[0][-1] == "RPS6"
    assert all(len(path) == len(set(path)) for path in paths)


def test_filters(G):
    H = filter_graph(G, evidence=["established"])
    assert all(d["evidence"] == "established" for _, _, d in H.edges(data=True))
    H2 = filter_graph(G, effects=["inhibits"])
    assert H2.number_of_edges() > 50
    M = by_module(G, "circadian")
    assert "BMAL1" in M


def test_search(G):
    assert "RPTOR" in search(G, "raptor")
    assert "AKT1S1" in search(G, "PRAS40")


def test_summary_and_hubs(G):
    s = summary(G)
    assert s["nodes"] == G.number_of_nodes()
    assert hubs(G, 3)[0][0] in ("MTORC1", "MTORC2", "AKT", "RPS6KB1")


def test_tree(G):
    t = build_tree(G, "MTORC1", direction="downstream", max_depth=3)
    assert t.id == "MTORC1" and t.children
    assert count(t) > 20
    txt = render_ascii(t)
    assert "MTORC1" in txt.splitlines()[0] or "mTORC1" in txt.splitlines()[0]
    md = render_markdown(t)
    assert md.startswith("- **")
    d = tree_to_dict(t)
    json.dumps(d)
    up = build_tree(G, "MTORC1", direction="upstream", max_depth=2)
    assert any(c.id == "RHEB" for c in up.children)


def test_tree_cycle_safe(G):
    t = build_tree(G, "RPS6KB1", direction="downstream", max_depth=8)
    assert count(t) < 5000


def test_exports(G, tmp_path):
    out = export.export_all(G, tmp_path)
    for k, p in out.items():
        assert p.exists() and p.stat().st_size > 0, k
    dot = (tmp_path / "mtor_interactome.dot").read_text(encoding="utf-8")
    assert "digraph" in dot and "arrowhead=tee" in dot


def test_draw_and_html(G, tmp_path):
    sub = by_module(G, "circadian")
    p = viz.draw_layered(sub, tmp_path / "circ.png", dpi=60)
    assert p.exists() and p.stat().st_size > 1000
    h = viz.write_html(sub, tmp_path / "circ.html")
    txt = h.read_text(encoding="utf-8")
    assert "vis-network" in txt and "BMAL1" in txt
