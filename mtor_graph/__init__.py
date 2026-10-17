"""mtor_graph – a curated interaction graph of the mTOR signalling pathway.

The package loads the knowledge base in ``data/`` (nodes, edges, references),
builds a :class:`networkx.MultiDiGraph`, and provides queries, hierarchical
"tree" views, exports (GraphML, GEXF, Cytoscape JSON, DOT, SIF, CSV) and
visualisations (static layered figure, interactive HTML).

Typical use::

    from mtor_graph import build_graph, build_tree, render_ascii
    G = build_graph()
    print(render_ascii(build_tree(G, "MTORC1", direction="downstream", max_depth=3)))
"""

from .loader import build_graph, load_edges, load_nodes, load_references, load_modules, to_simple, data_dir
from .queries import (
    upstream,
    downstream,
    neighborhood,
    shortest_paths,
    by_module,
    filter_graph,
    search,
    edges_between,
    summary,
    hubs,
)
from .tree import build_tree, render_ascii, render_markdown, tree_to_dict
from .schema import Node, Edge, validate

__all__ = [
    "build_graph",
    "load_edges",
    "load_nodes",
    "load_references",
    "load_modules",
    "to_simple",
    "data_dir",
    "upstream",
    "downstream",
    "neighborhood",
    "shortest_paths",
    "by_module",
    "filter_graph",
    "search",
    "edges_between",
    "summary",
    "hubs",
    "build_tree",
    "render_ascii",
    "render_markdown",
    "tree_to_dict",
    "Node",
    "Edge",
    "validate",
]

__version__ = "1.0.0"
