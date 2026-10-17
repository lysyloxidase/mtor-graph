# Data model

The knowledge base lives in four JSON files under `data/`:

| File | Content |
|---|---|
| `nodes.json` | every molecular entity, process, tissue, disease, lesion and drug (one object per node) |
| `edges.json` | every directed interaction (one object per record; a pair of nodes may carry several records when they state different claims) |
| `references.json` | the bibliography; `id` is the number used in `refs` |
| `modules.json` | the ten pathway modules (title, summary, colour group) |

## Node

```json
{
  "id": "TSC2",
  "label": "TSC2 (tuberin)",
  "type": "protein",
  "gene": "TSC2",
  "module": "growth_factor_input",
  "layer": 2,
  "description": "GAP for Rheb within the TSC1–TSC2–TBC1D7 complex …",
  "aliases": ["tuberin"]
}
```

* `id` — stable UPPER_SNAKE_CASE identifier. Human proteins use the HGNC
  symbol (`RPTOR`, `MAPKAP1`, `AKT1S1`, `RPS6KB1`, `EIF4EBP1`); complexes,
  metabolites, processes, tissues (`TISSUE_*`), diseases (`DISEASE_*`),
  phenotypes (`PHENO_*`) and gene lesions (`VAR_*`) use descriptive ids.
* `type` — one of `protein`, `complex`, `gene_variant`, `metabolite`, `lipid`,
  `drug`, `process`, `stimulus`, `phenotype`, `tissue`, `disease`, `organelle`,
  `rna`, `transcription_program`.
* `module` — the pathway module the node belongs to (see `modules.json`).
* `layer` — position in the hierarchical layout, 0 (stimuli, nutrients,
  drugs) … 8 (tissue physiology and clinical phenotypes). See
  `mtor_graph.schema.LAYERS`.

## Edge

```json
{
  "source": "AKT",
  "target": "TSC2",
  "interaction": "phosphorylation",
  "effect": "inhibits",
  "sites": ["Ser939", "Thr1462"],
  "mechanism": "AKT phosphorylates TSC2 and promotes displacement of the TSC complex from lysosomes, reducing local restraint of Rheb",
  "evidence": "established",
  "refs": [92, 93, 94],
  "context": "growth-factor stimulation",
  "direction": "upstream_of_mTORC1",
  "notes": null,
  "module": "growth_factor_input"
}
```

* `interaction` — controlled vocabulary (`mtor_graph.schema.INTERACTION_TYPES`):
  phosphorylation, dephosphorylation, binding, complex_formation, gap_activity,
  gef_activity, recruitment, sequestration, allosteric_activation,
  competitive_inhibition, ubiquitination, degradation, stabilisation,
  transcriptional_regulation, translational_regulation, transport, sensing,
  lipid_binding, scaffolding, drug_inhibition, drug_activation,
  metabolic_conversion, genetic_loss_of_function, genetic_gain_of_function,
  physiological_output, autoinhibition, acetylation, O_GlcNAcylation,
  membrane_targeting, other.
* `effect` — `activates`, `inhibits`, `modulates`, `required_for`,
  `component_of`, `converts_to`. The sign used for colouring and for
  sign-aware paths is +1 for activates/required_for, −1 for inhibits, 0
  otherwise.
* `evidence` — a recorded judgement: `established` (reproduced across
  laboratories), `emerging` (one laboratory, or not shown in vivo) or
  `contested` (conflicting results). It is a label, not a certificate.
* `refs` — bibliography entry numbers supporting the record. Every record
  cites at least one. A citation being present does not prove that it
  supports the exact sites, context or mechanism recorded.
* `direction` — the record's place in the pathway logic: `upstream_of_mTORC1`,
  `upstream_of_mTORC2`, `downstream_of_mTORC1`, `downstream_of_mTORC2`,
  `feedback`, `core`, `pharmacological`, `pathophysiological`, `circadian`.
* `context` — tissue, condition or species qualifier when the interaction is
  conditional. This is free text, and path queries do not check whether
  contexts are compatible across steps.
* `notes` — optional qualification of the evidence.

Two records with the same source, target and interaction type are distinct
claims when they differ in effect, sites, context, evidence, mechanism,
direction or module; `mtor_graph.schema.Edge.key` derives a stable key from
exactly those fields, so parallel claims never overwrite one another.

## Bibliography

```json
{"id": 92, "citation": "Manning BD, Tee AR, … Mol Cell. 2002;10(1):151-162.", "year": 2002, "pmid": "12150915", "doi": "10.1016/s1097-2765(02)00568-3"}
```

Entries are numbered from 1 without gaps, sorted by first author and year.
Every entry is cited by at least one record, and all but a few carry a PMID
and a DOI.

## Validation

`python scripts/validate_data.py` (or `mtor-graph validate`) checks the
controlled vocabularies, dangling endpoints, duplicate records, records
without a citation and citation numbers that are not in the bibliography.
`pytest` additionally asserts a set of canonical textbook interactions (see
`tests/test_data_integrity.py`) and the query semantics.

`python scripts/assemble_kb.py` normalises the data files in place after a
manual edit: canonical identifiers, controlled vocabularies, merging of
records that state the same claim, home module and layer per node, and a
stable sort order.
