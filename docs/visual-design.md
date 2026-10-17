# Visual encoding of the pathway graph

All figures (static PNG/SVG, DOT and the interactive HTML explorer) use one
encoding so that a reader can move between them without relearning anything.

## Node colour = module group (identity, categorical)

Ten modules are displayed in eight groups using a fixed categorical palette:

| Slot | Group | Modules | Light | Dark |
|---|---|---|---|---|
| 1 | Lysosomal nutrient sensing | `lysosomal_nutrient_sensing` | `#2a78d6` | `#3987e5` |
| 2 | Growth-factor, energy and stress inputs | `growth_factor_input`, `energy_stress_input` | `#eb6834` | `#d95926` |
| 3 | Core complexes | `core_complexes` | `#1baf7a` | `#199e70` |
| 4 | mTORC1 outputs | `mtorc1_output_translation_ribosome`, `mtorc1_output_autophagy_metabolism_feedback` | `#eda100` | `#c98500` |
| 5 | mTORC2 regulation and outputs | `mtorc2_regulation_and_outputs` | `#e87ba4` | `#d55181` |
| 6 | Circadian coupling | `circadian` | `#008300` | `#008300` |
| 7 | Tissue physiology and disease | `tissue_physiology_disease` | `#4a3aa7` | `#9085e9` |
| 8 | Pharmacology | `pharmacology` | `#e34948` | `#e66767` |

The repository does not include a colour-vision or contrast validation
report. Text labels identify nodes, and node *shape* encodes node type
(● protein/process, ■ complex/tissue/organelle, ◆ metabolite/lipid, ⬢ drug,
▲ stimulus, ▼ gene variant/RNA, ★ disease/phenotype). A table view of the
visible interactions is one click away in the HTML explorer.

## Edge colour + arrowhead = effect sign (polarity)

| Effect | Arrowhead | Colour (light / dark) |
|---|---|---|
| activates, required_for | ▶ | `#0ca30c` |
| inhibits | ┤ (bar) | `#d03b3b` / `#e66767` |
| modulates | ○ | `#898781` |
| component_of, binding, converts_to | none / ○ | `#898781` |

Green and red encode positive and negative effects and are reinforced by
arrowhead shape where supported by the renderer.

## Edge line style = evidence grade

solid = established · dashed = emerging · dotted = contested. These styles
display the recorded evidence labels; they are not a verification of them.

## Layout = layer (hierarchical "tree")

Nodes carry a `layer` 0–8 (stimuli → sensors → transducers → switches →
complexes → substrates → effectors → processes → physiology). The static
figure and the hierarchical HTML layout place layer 0 at the top and layer 8 at
the bottom, so the picture reads as a tree from inputs to outcomes; feedback
edges run upward.

## Chrome

System sans-serif; chart surface `#fcfcfb` (light) / `#1a1a19` (dark); primary
ink `#0b0b0b` / `#ffffff`; secondary ink `#52514e` / `#c3c2b7`; muted
`#898781`; hairlines `#e1e0d9` / `#2c2c2a`. The HTML explorer follows the
operating-system colour scheme and has a manual theme toggle.
