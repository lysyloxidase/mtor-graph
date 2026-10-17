"""Regression checks for claim preservation and graph-query semantics."""

from copy import deepcopy
from itertools import permutations

import networkx as nx
import pytest

from mtor_graph import Edge, Node, build_graph, build_tree, neighborhood, render_ascii, to_simple, validate
from mtor_graph import export, viz
from mtor_graph.loader import data_dir
from mtor_graph.queries import edges_between, path_sign, shortest_paths
from scripts.assemble_kb import canon_id, normalise


def node(nid):
    return Node(nid, nid, "protein", "core_complexes", 4)


def edge(**changes):
    fields = dict(source="A", target="B", interaction="phosphorylation", effect="activates",
                  refs=[1], sites=["Ser1"], module="core_complexes", mechanism="Recorded mechanism")
    fields.update(changes)
    return Edge(**fields)


def graph(*edges):
    ids = sorted({e.source for e in edges} | {e.target for e in edges})
    return build_graph(nodes=[node(n) for n in ids], edges=edges)


@pytest.mark.parametrize("changes", [
    {"effect": "inhibits"}, {"sites": ["Ser2"]}, {"context": "liver"},
    {"evidence": "contested"}, {"mechanism": "A different mechanism"},
    {"direction": "feedback"}, {"module": "growth_factor_input"},
])
def test_assembly_and_loader_preserve_distinct_claims(changes):
    first, second = edge(), edge(**changes)
    nodes, edges, audit = normalise(
        [node("A").to_dict(), node("B").to_dict()],
        [first.to_dict(), second.to_dict()], {1: {}},
    )
    assert len(edges) == 2
    assert audit["merged_duplicate_edges"] == 0
    built = build_graph(nodes=map(Node.from_dict, nodes), edges=map(Edge.from_dict, edges))
    assert built.number_of_edges() == 2


def test_assembly_merges_only_duplicate_claim_citations():
    first, second = edge().to_dict(), edge(refs=[2]).to_dict()
    _, edges, audit = normalise([], [first, second], {1: {}, 2: {}})
    assert len(edges) == 1
    assert audit["merged_duplicate_edges"] == 1
    assert edges[0]["refs"] == [1, 2]


def test_assembly_is_repeatable_without_modifying_inputs():
    nodes = [node("A").to_dict(), node("B").to_dict()]
    edges = [edge().to_dict(), edge(effect="inhibits", context="liver").to_dict()]
    original = deepcopy((nodes, edges))
    assert normalise(nodes, edges, {1: {}}) == normalise(nodes, edges, {1: {}})
    assert (nodes, edges) == original


def test_distinct_gene_is_not_a_global_alias():
    assert canon_id("XPLN") == "ARHGEF3"
    assert canon_id("ARHGEF26") == "ARHGEF26"


def test_loader_rejects_silent_data_loss():
    with pytest.raises(ValueError, match="duplicate edge"):
        graph(edge(), edge())
    with pytest.raises(ValueError, match="unknown endpoint"):
        build_graph(nodes=[node("A")], edges=[edge()])
    with pytest.raises(ValueError, match="duplicate node"):
        build_graph(nodes=[node("A"), node("A")], edges=[])


def test_conflicting_signs_never_vote_or_depend_on_order():
    claims = [edge(sites=["Ser1"]), edge(sites=["Ser2"]), edge(effect="inhibits")]
    for ordered in permutations(claims):
        g = graph(*ordered)
        collapsed = to_simple(g)["A"]["B"]
        assert collapsed["sign"] == 0
        assert collapsed["sign_conflict"] is True
        assert collapsed["effect"] == "modulates"
        assert path_sign(g, ["A", "B"]) == 0


def test_collapsed_evidence_does_not_upgrade_contested_claim():
    collapsed = to_simple(graph(edge(), edge(evidence="contested", refs=[2])))["A"]["B"]
    assert collapsed["evidence"] == "contested"
    assert collapsed["refs"] == [1, 2]


def test_tree_preserves_parallel_claims_and_filters_before_selection():
    g = graph(edge(), edge(effect="inhibits", evidence="contested", refs=[2],
                           sites=["Ser2"], context="liver"))
    tree = build_tree(g, "A", max_depth=1)
    assert len(tree.children) == 2
    assert {(c.effect, tuple(c.sites), tuple(c.refs)) for c in tree.children} == {
        ("activates", ("Ser1",), (1,)), ("inhibits", ("Ser2",), (2,))}
    filtered = build_tree(g, "A", max_depth=1, evidence=["contested"])
    assert len(filtered.children) == 1
    assert filtered.children[0].effect == "inhibits"
    text = render_ascii(filtered, show_refs=True)
    assert "liver" in text and "[2]" in text


def test_tree_repeated_marker_does_not_assert_cycle():
    g = graph(edge(), edge(source="A", target="C"), edge(source="B", target="D"),
              edge(source="C", target="D"))
    tree = build_tree(g, "A", max_depth=2)
    assert nx.is_directed_acyclic_graph(g)
    text = render_ascii(tree, ascii_only=True)
    assert "(repeated)" in text and "(loop)" not in text


def test_neighborhood_can_change_direction_at_each_step():
    g = graph(edge(), edge(source="C", target="B"))
    assert set(neighborhood(g, "A", depth=2)) == {"A", "B", "C"}
    assert set(neighborhood(g, "A", depth=1)) == {"A", "B"}


def test_paths_and_sign_product_are_separate_operations():
    g = graph(edge(effect="inhibits"), edge(source="B", target="C", effect="inhibits"))
    assert shortest_paths(g, "A", "C", k=2) == [["A", "B", "C"]]
    assert shortest_paths(g, "C", "A") == []
    assert shortest_paths(g, "A", "C", k=0) == []
    assert path_sign(g, ["A", "B", "C"]) == 1
    with pytest.raises(nx.NodeNotFound):
        shortest_paths(g, "missing", "C")


def test_self_loop_is_reported_once():
    assert len(edges_between(graph(edge(target="A")), "A", "A")) == 1


def test_reference_existence_does_not_resolve_invalid_source():
    problems = validate([node("A"), node("B")], [edge()], {1: {"unresolved": True}})
    assert len(problems) == 1 and "reference [1]" in problems[0]


def test_explicit_data_directory_overrides_environment(monkeypatch, tmp_path):
    monkeypatch.setenv("MTOR_GRAPH_DATA", str(tmp_path / "environment"))
    assert data_dir() == tmp_path / "environment"
    assert data_dir(tmp_path / "explicit") == tmp_path / "explicit"


def test_exports_preserve_claims_and_svg_is_repeatable(tmp_path):
    g = graph(edge(), edge(effect="inhibits", context="liver"))
    gexf = export.write_gexf(g, tmp_path / "graph.gexf")
    assert "lastmodifieddate" not in gexf.read_text(encoding="utf-8")
    assert nx.read_gexf(gexf).number_of_edges() == 2
    svg1 = viz.draw_layered(g, tmp_path / "one.svg")
    svg2 = viz.draw_layered(g, tmp_path / "two.svg")
    assert svg1.read_bytes() == svg2.read_bytes()


def test_html_title_is_text(tmp_path):
    path = viz.write_html(graph(edge()), tmp_path / "graph.html", title="<b>Graph</b>")
    assert "<title>&lt;b&gt;Graph&lt;/b&gt;</title>" in path.read_text(encoding="utf-8")
