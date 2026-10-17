"""Data model and validation for the mTOR interaction knowledge base.

Every node and edge in ``data/nodes.json`` / ``data/edges.json`` follows the
schema below. Controlled vocabularies are enforced by :func:`validate`.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
import hashlib
import json
from typing import Any, Iterable

INTERACTION_TYPES: frozenset[str] = frozenset(
    {
        "phosphorylation",
        "dephosphorylation",
        "binding",
        "complex_formation",
        "gap_activity",
        "gef_activity",
        "recruitment",
        "sequestration",
        "allosteric_activation",
        "competitive_inhibition",
        "ubiquitination",
        "degradation",
        "stabilisation",
        "transcriptional_regulation",
        "translational_regulation",
        "transport",
        "sensing",
        "lipid_binding",
        "scaffolding",
        "drug_inhibition",
        "drug_activation",
        "metabolic_conversion",
        "genetic_loss_of_function",
        "genetic_gain_of_function",
        "physiological_output",
        "autoinhibition",
        "acetylation",
        "O_GlcNAcylation",
        "membrane_targeting",
        "other",
    }
)

EFFECTS: frozenset[str] = frozenset(
    {"activates", "inhibits", "modulates", "required_for", "component_of", "converts_to"}
)

EVIDENCE_LEVELS: frozenset[str] = frozenset({"established", "emerging", "contested"})

NODE_TYPES: frozenset[str] = frozenset(
    {
        "protein",
        "complex",
        "gene_variant",
        "metabolite",
        "lipid",
        "drug",
        "process",
        "stimulus",
        "phenotype",
        "tissue",
        "disease",
        "organelle",
        "rna",
        "transcription_program",
    }
)

DIRECTIONS: frozenset[str] = frozenset(
    {
        "upstream_of_mTORC1",
        "upstream_of_mTORC2",
        "downstream_of_mTORC1",
        "downstream_of_mTORC2",
        "feedback",
        "core",
        "pharmacological",
        "pathophysiological",
        "circadian",
    }
)

MODULE_KEYS: tuple[str, ...] = (
    "lysosomal_nutrient_sensing",
    "growth_factor_input",
    "energy_stress_input",
    "core_complexes",
    "mtorc1_output_translation_ribosome",
    "mtorc1_output_autophagy_metabolism_feedback",
    "mtorc2_regulation_and_outputs",
    "circadian",
    "tissue_physiology_disease",
    "pharmacology",
)

LAYERS: dict[int, str] = {
    0: "Stimuli, nutrients and drugs",
    1: "Receptors, first messengers and sensors",
    2: "Upstream transducers",
    3: "Proximal switches",
    4: "mTOR complexes and subunits",
    5: "Direct substrates",
    6: "Secondary effectors",
    7: "Cellular processes and programmes",
    8: "Tissue physiology and clinical phenotypes",
}

# Positive / negative sign of an effect, used for colouring and for sign-aware paths.
EFFECT_SIGN: dict[str, int] = {
    "activates": +1,
    "required_for": +1,
    "inhibits": -1,
    "modulates": 0,
    "component_of": 0,
    "converts_to": 0,
}


@dataclass
class Node:
    id: str
    label: str
    type: str
    module: str
    layer: int
    description: str = ""
    gene: str | None = None
    aliases: list[str] = field(default_factory=list)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Node":
        return cls(
            id=str(d["id"]),
            label=str(d.get("label") or d["id"]),
            type=str(d.get("type", "protein")),
            module=str(d.get("module", "")),
            layer=int(d.get("layer", 4)),
            description=str(d.get("description") or ""),
            gene=(d.get("gene") if d.get("gene") not in ("", None) else None),
            aliases=list(d.get("aliases") or []),
        )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class Edge:
    source: str
    target: str
    interaction: str
    effect: str
    mechanism: str = ""
    evidence: str = "established"
    sites: list[str] = field(default_factory=list)
    refs: list[int] = field(default_factory=list)
    context: str | None = None
    direction: str = "core"
    notes: str | None = None
    module: str = ""

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Edge":
        return cls(
            source=str(d["source"]),
            target=str(d["target"]),
            interaction=str(d.get("interaction", "other")),
            effect=str(d.get("effect", "modulates")),
            mechanism=str(d.get("mechanism") or ""),
            evidence=str(d.get("evidence", "established")),
            sites=[str(s) for s in (d.get("sites") or [])],
            refs=[int(r) for r in (d.get("refs") or [])],
            context=(d.get("context") or None),
            direction=str(d.get("direction", "core")),
            notes=(d.get("notes") or None),
            module=str(d.get("module") or ""),
        )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @property
    def key(self) -> str:
        """Stable key for a claim, including its context and evidence grade.

        Citations and notes do not define identity. Distinct mechanisms, sites,
        effects or contexts must not overwrite one another in a multigraph.
        """
        claim = [self.interaction, self.effect, sorted(set(self.sites)), self.context,
                 self.evidence, self.direction, self.module, self.mechanism]
        digest = hashlib.sha256(json.dumps(claim, ensure_ascii=False).encode("utf-8")).hexdigest()[:16]
        return f"{self.interaction}:{digest}"

    @property
    def sign(self) -> int:
        return EFFECT_SIGN.get(self.effect, 0)


def validate(
    nodes: Iterable[Node],
    edges: Iterable[Edge],
    references: dict[int, Any] | None = None,
) -> list[str]:
    """Return a list of human-readable problems (empty list == valid)."""
    nodes = list(nodes)
    edges = list(edges)
    problems: list[str] = []
    ids: dict[str, int] = {}
    for n in nodes:
        ids[n.id] = ids.get(n.id, 0) + 1
        if n.type not in NODE_TYPES:
            problems.append(f"node {n.id}: unknown type {n.type!r}")
        if n.module not in MODULE_KEYS:
            problems.append(f"node {n.id}: unknown module {n.module!r}")
        if n.layer not in LAYERS:
            problems.append(f"node {n.id}: layer {n.layer} outside {min(LAYERS)}..{max(LAYERS)}")
        if not n.id.replace("_", "").replace("-", "").isalnum():
            problems.append(f"node {n.id}: id must be alphanumeric/underscore")
    for nid, c in ids.items():
        if c > 1:
            problems.append(f"node {nid}: duplicated {c} times")

    seen: dict[tuple[str, str, str], int] = {}
    for i, e in enumerate(edges):
        where = f"edge[{i}] {e.source}->{e.target} ({e.interaction})"
        if e.source not in ids:
            problems.append(f"{where}: unknown source node")
        if e.target not in ids:
            problems.append(f"{where}: unknown target node")
        if e.interaction not in INTERACTION_TYPES:
            problems.append(f"{where}: unknown interaction type {e.interaction!r}")
        if e.effect not in EFFECTS:
            problems.append(f"{where}: unknown effect {e.effect!r}")
        if e.evidence not in EVIDENCE_LEVELS:
            problems.append(f"{where}: unknown evidence level {e.evidence!r}")
        if e.direction not in DIRECTIONS:
            problems.append(f"{where}: unknown direction {e.direction!r}")
        if e.module and e.module not in MODULE_KEYS:
            problems.append(f"{where}: unknown module {e.module!r}")
        if not e.refs:
            problems.append(f"{where}: no citation")
        if references is not None:
            for r in e.refs:
                if r not in references:
                    problems.append(f"{where}: reference [{r}] is not in the bibliography")
                elif references[r].get("unresolved"):
                    problems.append(f"{where}: reference [{r}] is marked unresolved and may not be cited")
        k = (e.source, e.target, e.key)
        seen[k] = seen.get(k, 0) + 1
    for k, c in seen.items():
        if c > 1:
            problems.append(f"edge {k[0]}->{k[1]} [{k[2]}]: duplicated {c} times")
    return problems
