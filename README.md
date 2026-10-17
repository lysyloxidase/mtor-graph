# mTOR signalling interactome

A curated, machine-readable graph of the mTOR signalling pathway: the
lysosomal nutrient-sensing module, growth-factor, energy and stress inputs,
the architecture of mTORC1 and mTORC2, their substrates and output
programmes, circadian coupling, tissue physiology, disease lesions and the
pharmacology of mTOR inhibitors.

Each interaction carries its direction, effect sign, interaction type,
phosphosites, a one-sentence mechanism, a tissue or condition context, an
evidence label and the publications supporting it. The repository ships the
data, a Python package built on [networkx](https://networkx.org/) to query
and visualise it, a command-line tool, exports in the usual network formats,
pre-rendered figures and an interactive HTML explorer.

<!--STATS-->
**545 nodes · 1128 directed source–target pairs · 1284 claim records** (a pair can carry several records when they differ in interaction type, sign, sites, context or evidence). Every record cites the bibliography (708 distinct publications); recorded evidence grades: 1047 established, 205 emerging, 32 contested. Node types: 249 protein, 75 phenotype, 34 complex, 34 process, 31 disease, 27 drug, 26 gene variant, 22 stimulus, 18 metabolite, 16 tissue, 5 lipid, 4 rna, 3 organelle, 1 transcription program. Claim records per module: Tissue physiology and disease 240, mTORC1 output: autophagy, metabolism and feedback 153, Pharmacology 146, Lysosomal nutrient sensing 129, mTORC2 regulation and outputs 126, mTORC1 output: translation and ribosome biogenesis 111, Core complexes (mTORC1 and mTORC2 architecture) 110, Energy, oxygen and stress input (AMPK, REDD1, p53) 99, Growth-factor input (PI3K–AKT–TSC–Rheb) 97, Circadian coupling 73.
<!--/STATS-->

**Po polsku.** Repozytorium zawiera bazę interakcji szlaku mTOR w formacie
nadającym się do analizy maszynowej oraz pakiet Pythona do budowy grafu
skierowanego, zapytań, drzew hierarchicznych i eksportu (GraphML, GEXF,
Cytoscape, DOT, SIF, CSV), gotowe rysunki i interaktywny eksplorator HTML.
Każdy wpis ma typ interakcji, znak efektu, miejsca modyfikacji, kontekst
tkankowy, etykietę siły dowodów i odwołania bibliograficzne z PMID lub DOI.
Etykieta dowodów jest oceną kuratorską, a nie potwierdzeniem poprawności.
Szybki start: sekcja *Quick start*.

## Contents

| Path | What it is |
|---|---|
| `data/nodes.json`, `data/edges.json` | the knowledge base (see [docs/data-model.md](docs/data-model.md)) |
| `data/references.json` | the bibliography, numbered from 1, every entry cited |
| `data/modules.json` | the ten pathway modules |
| `mtor_graph/` | Python package: loader, queries, tree views, exports, figures, CLI |
| `output/mtor_interactome.html` | interactive explorer: filter by module, evidence and effect, hierarchical or force layout, node details with citations, table view |
| `output/figures/` | layered pathway figures: the whole graph, one per module, and the direct interactors of mTORC1 and mTORC2 |
| `output/trees/` | text and Markdown tree views rooted at mTORC1, mTORC2, leucine, insulin, rapamycin and AMPK |
| `output/exports/` | GraphML, GEXF, Cytoscape JSON, Graphviz DOT, SIF, CSV, JSON |
| `output/stats.md` | descriptive statistics |
| `tests/` | data integrity, canonical interactions and query semantics |
| `scripts/` | `build_all.py` (regenerate `output/`), `validate_data.py`, `assemble_kb.py` (normalise `data/` in place) |

## Quick start

```bash
pip install -e .
mtor-graph stats
mtor-graph tree --root MTORC1 --direction downstream --depth 3
mtor-graph tree --root MTORC1 --direction upstream --depth 3 --refs
mtor-graph node RHEB --refs
mtor-graph path LEUCINE RPS6
mtor-graph search sestrin
mtor-graph module lysosomal_nutrient_sensing
mtor-graph export output/exports
mtor-graph draw output/figures/full.png
mtor-graph html output/mtor_interactome.html
python scripts/build_all.py        # regenerate everything in output/
pytest                             # integrity and query tests
```

`--evidence established,emerging` before any subcommand restricts the graph
to records with those labels; `--modules circadian,pharmacology` restricts it
to modules.

### Python

```python
from mtor_graph import build_graph, build_tree, render_ascii, shortest_paths, upstream, by_module

G = build_graph()                       # networkx.MultiDiGraph; edge attributes = the data model
print(G.number_of_nodes(), G.number_of_edges())

# all regulators within two steps of mTORC1, with distances
print(upstream(G, "MTORC1", depth=2))

# shortest directed paths (topology only; contexts may differ between steps)
for p in shortest_paths(G, "INSULIN", "AUTOPHAGY", k=3):
    print(" -> ".join(p))

# a hierarchical view of what mTORC2 acts on
print(render_ascii(build_tree(G, "MTORC2", direction="downstream", max_depth=3), show_refs=True))

# a sub-graph of one module
C = by_module(G, "circadian")
```

Abbreviated illustration of a downstream tree (not the complete output):

```
mTORC1
├── ──▶ S6K1  (Thr389; phosphorylation)
│   ├── ──┤ IRS-1  (phosphorylation)
│   ├── ──▶ rpS6  (Ser235, Ser236, Ser240, Ser244; phosphorylation)
│   └── ──┤ Rictor  (Thr1135; phosphorylation)
├── ──┤ 4E-BP1  (Thr37, Thr46, Ser65, Thr70; phosphorylation)
│   └── ──┤ eIF4E  (binding, sequestration)
├── ──┤ ULK1  (Ser757; phosphorylation)
└── ──┤ TFEB  (Ser211; phosphorylation)
```

(`──▶` activates, `──┤` inhibits, `──◇` modulates, `──■` component of;
`(repeated)` marks a node already expanded elsewhere, which can be a shared
descendant or a feedback loop.)

## The modules

| Module | Scope |
|---|---|
| `lysosomal_nutrient_sensing` | Rag GTPases, Ragulator, v-ATPase, GATOR1/2, KICSTOR, FLCN–FNIP, Sestrin2, CASTOR1, SAMTOR, SLC38A9, LARS1, SAR1B |
| `growth_factor_input` | RTK → IRS-1 → PI3K → PDK1/AKT ⊣ TSC ⊣ Rheb; PRAS40; ERK/RSK, Wnt–GSK3, IKKβ |
| `energy_stress_input` | LKB1–AMPK (TSC2, Raptor Ser792, ULK1), hypoxia–REDD1, DNA damage–p53, metformin, exercise, caloric restriction |
| `core_complexes` | subunits and domains of mTORC1 and mTORC2, TOS and RAIP motifs, FRB occlusion, mSin1 autoinhibition, distal Akt docking, DEPTOR and PRAS40 brakes |
| `mtorc1_output_translation_ribosome` | S6K1 and its substrates, 4E-BP1–eIF4E, LARP1 and TOP mRNAs, TIF-IA/UBF, MAF1 |
| `mtorc1_output_autophagy_metabolism_feedback` | ULK1/ATG13/TFEB/TFE3, HIF-1α, SREBP/Lipin1, G6PD, CAD, ATF4–MTHFD2, c-Myc, PGC-1α; S6K1–IRS-1 and Grb10 feedback |
| `mtorc2_regulation_and_outputs` | PIP3–mSin1, ribosome association, AKT Ser473/Thr450, SGK1, PKC isoforms, FOXO, ENaC, Tiam1–Rac1 |
| `circadian` | BMAL1 and PER2 ⊣ mTORC1, S6K1 → BMAL1, SCN entrainment, rhythmic autophagy |
| `tissue_physiology_disease` | tissue-specific outputs, gene lesions, mTORopathies, cancer lesion classes, neurodegeneration, immunity, ageing |
| `pharmacology` | rapalogs, ATP-competitive, dual PI3K/mTOR and bi-steric inhibitors: selectivity, feedback, adverse effects, indications |

Nodes carry a hierarchical `layer` (0 stimuli, nutrients and drugs → 4 the
complexes → 8 tissue physiology and clinical phenotypes), which is what the
layered figures and the tree layout use. The visual encoding is documented in
[docs/visual-design.md](docs/visual-design.md).

## Provenance and limitations

The graph is a reading of the published literature, not a measurement.

- Every record cites at least one publication, and all but a few citations
  carry a PMID and a DOI. A citation being present does not prove that it
  supports the exact sites, context or mechanism recorded; check the primary
  source before relying on a record, including one labelled `established`.
- Evidence labels are curation judgements, not independently validated
  scores.
- Coverage is extensive but not exhaustive. This is not a complete
  interactome.
- Paths join evidence from different organisms, tissues and experimental
  conditions. The product of edge signs along a path is not a prediction of a
  treatment's net effect.
- Drugs, doses and adverse effects appear here for scientific synthesis only.
  Nothing in this repository is clinical guidance.

## Licence

Code: MIT. Data files in `data/`: CC BY 4.0. See `LICENSE`.
