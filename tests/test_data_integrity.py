"""Integrity and 'golden interaction' tests for the knowledge base."""

import pytest

from mtor_graph import build_graph, load_edges, load_nodes, load_references, validate, edges_between
from mtor_graph.schema import EFFECT_SIGN


@pytest.fixture(scope="module")
def G():
    return build_graph()


def test_validate_clean():
    problems = validate(load_nodes(), load_edges(), load_references())
    assert problems == [], "\n".join(problems)


def test_bibliography_is_numbered_without_gaps():
    refs = load_references()
    assert set(refs) == set(range(1, len(refs) + 1))
    assert all(r.get("citation") for r in refs.values())


CANONICAL_NODES = [
    "MTOR", "MTORC1", "MTORC2", "RPTOR", "RICTOR", "MLST8", "MAPKAP1", "DEPTOR", "AKT1S1",
    "RHEB", "TSC1", "TSC2", "TBC1D7", "AKT", "PDK1", "PI3K", "PTEN", "PIP3", "AMPK", "STK11",
    "RRAGA", "RRAGC", "RAGULATOR", "GATOR1", "GATOR2", "KICSTOR", "FLCN", "SESN2", "CASTOR1", "SAMTOR", "SLC38A9",
    "RPS6KB1", "EIF4EBP1", "EIF4E", "ULK1", "TFEB", "LARP1", "MAF1", "SREBF1", "LPIN1", "HIF1A", "ATF4", "CAD",
    "IRS1", "GRB10", "SGK1", "FOXO1", "BMAL1", "PER2", "RAPAMYCIN", "FKBP12", "EVEROLIMUS", "RMC5552",
]


@pytest.mark.parametrize("node", CANONICAL_NODES)
def test_canonical_node_present(G, node):
    assert node in G, f"{node} missing from the knowledge base"


# (sources, targets, expected sign, required site substring or None)
GOLDEN = [
    (["TSC2", "TSC_COMPLEX"], ["RHEB"], -1, None),
    (["RHEB"], ["MTORC1"], +1, None),
    (["AKT"], ["TSC2", "TSC_COMPLEX"], -1, None),
    (["AKT"], ["AKT1S1"], -1, "246"),
    (["AKT1S1"], ["MTORC1", "RPTOR"], -1, None),
    (["MTORC2"], ["AKT"], +1, "473"),
    (["PDK1"], ["AKT"], +1, "308"),
    (["MTORC1"], ["RPS6KB1"], +1, "389"),
    (["MTORC1"], ["EIF4EBP1"], -1, None),
    (["EIF4EBP1"], ["EIF4E"], -1, None),
    (["MTORC1"], ["ULK1"], -1, "757"),
    (["AMPK"], ["RPTOR", "MTORC1"], -1, "792"),
    (["AMPK"], ["TSC2", "TSC_COMPLEX"], +1, None),
    (["AMPK"], ["ULK1"], +1, None),
    (["MTORC1"], ["TFEB"], -1, "211"),
    (["GATOR1"], ["RRAGA", "RRAGB", "RAG_HETERODIMER"], -1, None),
    (["GATOR2"], ["GATOR1"], -1, None),
    (["SESN2"], ["GATOR2"], -1, None),
    (["CASTOR1"], ["GATOR2"], -1, None),
    (["SAMTOR"], ["GATOR1", "KICSTOR"], +1, None),
    (["DDIT4"], ["TSC2", "TSC_COMPLEX"], +1, None),
    (["RPS6KB1"], ["IRS1"], -1, None),
    (["RPS6KB1"], ["RICTOR", "MTORC2"], -1, None),
    (["MTORC1", "MTORC2"], ["DEPTOR"], -1, None),
    (["RAPAMYCIN", "FKBP12", "FKBP12_RAPAMYCIN"], ["MTORC1", "FRB_DOMAIN", "MTOR"], -1, None),
    (["BMAL1"], ["MTORC1"], -1, None),
    (["PER2"], ["MTORC1", "TSC1"], None, None),
    (["RPS6KB1"], ["BMAL1"], None, None),
    (["MTORC2"], ["SGK1"], +1, "422"),
    (["PTEN"], ["PIP3"], -1, None),
    (["MTORC1"], ["LARP1"], -1, None),
    (["MTORC1"], ["MAF1"], -1, None),
    (["MTORC1"], ["GRB10"], +1, None),
    (["MTORC1", "RPS6KB1"], ["CAD"], +1, "1859"),
    (["RPS6KB1"], ["RPS6"], +1, None),
    (["RRAGA", "RAG_HETERODIMER", "RRAGC"], ["RPTOR", "MTORC1"], +1, None),
    (["FLCN", "FLCN_FNIP"], ["RRAGC", "RRAGD", "RAG_HETERODIMER"], None, None),
    (["RAGULATOR"], ["RRAGA", "RRAGB", "RAG_HETERODIMER"], +1, None),
    (["PIP3"], ["MAPKAP1", "MTORC2"], +1, None),
    (["RPS6KB1"], ["EIF4B"], +1, "422"),
    (["RPS6KB1"], ["PDCD4"], -1, None),
    (["RPS6KB1", "MTORC1"], ["EEF2K"], -1, None),
    (["TRPML1", "MCOLN1", "PPP3", "CALCINEURIN"], ["TFEB"], +1, None),
]


@pytest.mark.parametrize("sources,targets,sign,site", GOLDEN, ids=[f"{s[0]}->{t[0]}" for s, t, _, _ in GOLDEN])
def test_golden_interaction(G, sources, targets, sign, site):
    hits = []
    for s in sources:
        for t in targets:
            if s in G and t in G and G.has_edge(s, t):
                for d in G[s][t].values():
                    ok_sign = sign is None or EFFECT_SIGN.get(d["effect"], 0) == sign
                    ok_site = site is None or any(site in x for x in d.get("sites", []))
                    if ok_sign and ok_site:
                        hits.append(d)
    assert hits, f"no edge {sources} -> {targets} with sign {sign} and site {site}"


def test_every_edge_cites_the_bibliography(G):
    refs = G.graph["references"]
    for u, v, d in G.edges(data=True):
        assert d.get("refs"), f"{u}->{v} has no citation"
        for r in d["refs"]:
            assert r in refs, f"{u}->{v} cites [{r}], which is not in the bibliography"


def test_every_bibliography_entry_is_cited(G):
    cited = {r for _, _, d in G.edges(data=True) for r in d.get("refs", [])}
    assert set(G.graph["references"]) == cited


def test_core_layers(G):
    assert G.nodes["MTORC1"]["layer"] == 4
    assert G.nodes["MTORC2"]["layer"] == 4
    assert G.nodes["RHEB"]["layer"] <= 3
    assert G.nodes["RPS6KB1"]["layer"] >= 5


def test_feedback_loops_exist(G):
    """The pathway is cyclic: S6K1 → IRS-1 ⊣ PI3K → … → mTORC1 → S6K1."""
    import networkx as nx

    H = nx.DiGraph(G)
    assert "MTORC1" in nx.descendants(H, "RPS6KB1") and "RPS6KB1" in nx.descendants(H, "MTORC1")
    assert "IRS1" in nx.descendants(H, "RPS6KB1") and "RPS6KB1" in nx.descendants(H, "IRS1")
