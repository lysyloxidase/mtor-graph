"""Normalise ``data/nodes.json`` and ``data/edges.json`` in place.

Maps identifiers to their canonical form, enforces the controlled vocabularies of
``docs/data-model.md``, merges records that state the same claim, resolves each node to
its home module and layout layer, drops records without a usable citation, and sorts
both files. Run it after editing the data by hand; it is idempotent.

Usage::

    python scripts/assemble_kb.py
"""

from __future__ import annotations

import argparse
import collections
import json
import re
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from mtor_graph.schema import (  # noqa: E402
    DIRECTIONS, EFFECTS, EVIDENCE_LEVELS, INTERACTION_TYPES, MODULE_KEYS, NODE_TYPES, Edge,
)

MODULE_ORDER = list(MODULE_KEYS)

# --------------------------------------------------------------------------- canonical ids
ID_ALIASES: dict[str, str] = {
    "RAPTOR": "RPTOR", "S6K1": "RPS6KB1", "S6K": "RPS6KB1", "P70S6K": "RPS6KB1", "P70S6K1": "RPS6KB1", "S6K2": "RPS6KB2",
    "4EBP1": "EIF4EBP1", "4E_BP1": "EIF4EBP1", "4EBP2": "EIF4EBP2", "4E_BP2": "EIF4EBP2", "PRAS40": "AKT1S1",
    "MSIN1": "MAPKAP1", "SIN1": "MAPKAP1", "FKBP1A": "FKBP12", "BMT2": "SAMTOR", "ARNTL": "BMAL1", "PRKAA1": "AMPK",
    "PRKAA2": "AMPK", "LKB1": "STK11", "TRPML1": "MCOLN1", "CALCINEURIN": "PPP3", "PPP3CA": "PPP3", "PPP3CB": "PPP3",
    "REDD1": "DDIT4", "TIF_IA": "RRN3", "TIFIA": "RRN3", "TIF1A": "RRN3", "UBF": "UBTF", "GDH": "GLUD1", "PGC1A": "PPARGC1A",
    "PGC_1A": "PPARGC1A", "GLUT1": "SLC2A1", "SREBP1": "SREBF1", "SREBP1C": "SREBF1", "SREBP2": "SREBF2", "SREBP": "SREBF1",
    "LIPIN1": "LPIN1", "ENAC": "SCNN1", "SCNN1A": "SCNN1", "PKB": "AKT", "AKT1": "AKT", "ERK1_2": "ERK", "MAPK1": "ERK",
    "MAPK3": "ERK", "ERK1": "ERK", "ERK2": "ERK", "RPS6KA1": "RSK", "RSK1": "RSK", "P90RSK": "RSK", "IKKB": "IKBKB",
    "IKK_BETA": "IKBKB", "P53": "TP53", "SIROLIMUS": "RAPAMYCIN", "P_GLYCOPROTEIN": "ABCB1", "PGP": "ABCB1", "P_GP": "ABCB1",
    "14_3_3": "YWHA", "YWHAZ": "YWHA", "FIP200": "RB1CC1", "ATG14L": "ATG14", "POLDIP3": "SKAR", "LAMTOR": "RAGULATOR",
    "ATP6V": "V_ATPASE", "VATPASE": "V_ATPASE", "FLCN_FNIP1": "FLCN_FNIP", "FLCN_FNIP2": "FLCN_FNIP", "TSC1_TSC2": "TSC_COMPLEX",
    "TSC": "TSC_COMPLEX", "TSC1_TSC2_TBC1D7": "TSC_COMPLEX", "RRAGA_RRAGC": "RAG_HETERODIMER", "RAGS": "RAG_HETERODIMER",
    "RAG_GTPASES": "RAG_HETERODIMER", "RAG": "RAG_HETERODIMER", "MTORC1_COMPLEX": "MTORC1", "MTORC2_COMPLEX": "MTORC2",
    "PIK3CA": "PI3K", "HIF_1A": "HIF1A", "HIF1": "HIF1A", "C_MYC": "MYC", "CMYC": "MYC", "CYCLIN_D1": "CCND1", "VEGF": "VEGFA",
    "FOXO3A": "FOXO3", "NA_K_ATPASE": "ATP1A1", "NAK_ATPASE": "ATP1A1", "SESTRIN2": "SESN2", "SESTRIN1": "SESN1",
    "SESTRIN3": "SESN3", "CASTOR": "CASTOR1", "RHEB1": "RHEB", "INSULIN_RECEPTOR": "INSR", "IGF_1": "IGF1", "IGF_1R": "IGF1R",
    "EIF4F_COMPLEX": "EIF4F", "PROTOR1": "PRR5", "PROTOR2": "PRR5L", "PROTOR_1": "PRR5", "PROTOR_2": "PRR5L",
    "BETA_TRCP": "BTRC", "BTRCP": "BTRC", "SAM_SENSOR": "SAMTOR", "S_ADENOSYLMETHIONINE": "SAM", "MO25": "CAB39",
    "TORIN_1": "TORIN1", "TORIN_2": "TORIN2", "RAPALINK_1": "RAPALINK1", "RMC_5552": "RMC5552", "KU_0063794": "KU0063794",
    "PP_242": "PP242", "NAB_SIROLIMUS": "NAB_SIROLIMUS", "AMINO_ACID": "AMINO_ACIDS", "EEF2_KINASE": "EEF2K", "EIF4G": "EIF4G1",
    "EIF4A": "EIF4A1", "ULK1_COMPLEX": "ULK1", "LAMP": "LYSOSOME", "LAT1": "SLC7A5", "FKBP12_RAPAMYCIN_COMPLEX": "FKBP12_RAPAMYCIN",
    "FKBP12_RAPA": "FKBP12_RAPAMYCIN", "RAPAMYCIN_FKBP12": "FKBP12_RAPAMYCIN", "HSP90": "HSP90AA1", "TEL2": "TELO2",
    "CLEAR": "CLEAR_NETWORK", "PKC_ISOFORMS": "PKC", "PKCALPHA": "PRKCA", "PKC_ALPHA": "PRKCA", "FOXO": "FOXO1",
    "GSK3": "GSK3B", "ATF4_TARGETS": "ATF4", "AMP_ATP_RATIO": "AMP", "LEU": "LEUCINE", "ARG": "ARGININE",
    # structural sub-domains are modelled as their own nodes without the parent-protein prefix
    "MSIN1_CRIM_DOMAIN": "CRIM_DOMAIN", "MSIN1_PH_DOMAIN": "PH_DOMAIN", "MSIN1_RBD_DOMAIN": "RBD_DOMAIN",
    "MAPKAP1_CRIM_DOMAIN": "CRIM_DOMAIN", "MAPKAP1_PH_DOMAIN": "PH_DOMAIN", "MAPKAP1_RBD_DOMAIN": "RBD_DOMAIN",
    "MTOR_FATC_DOMAIN": "FATC_DOMAIN", "MTOR_FAT_DOMAIN": "FAT_DOMAIN", "MTOR_HEAT_REPEATS": "HEAT_REPEATS",
    "MTOR_KINASE_DOMAIN": "KINASE_DOMAIN", "MTOR_FRB_DOMAIN": "FRB_DOMAIN",
    # XPLN is the exchange factor ARHGEF3 (ARHGEF26, i.e. SGEF, is a different protein)
    "XPLN": "ARHGEF3",
    # one concept, one id
    "PHENO_SASP": "SASP", "PHENO_INSULIN_SECRETION": "INSULIN_SECRETION", "DISEASE_ALZHEIMER": "DISEASE_AD",
    "DISEASE_HUNTINGTON": "DISEASE_HD", "DISEASE_PARKINSON": "DISEASE_PD", "DISEASE_ALZHEIMERS_DISEASE": "DISEASE_AD",
    "PHENO_DELAYED_SLEEP_PHASE": "PHENO_DELAYED_SLEEP_WAKE_PHASE_DISORDER", "DISEASE_NEUROENDOCRINE_TUMOUR": "DISEASE_NET",
    "PHENO_ANTIVIRAL_GENE_EXPRESSION": "ISG_ANTIVIRAL_PROGRAMME", "PHENO_LTP": "LTP_LATE_PHASE", "PHENO_TH2": "TH2_DIFFERENTIATION",
    "PHENO_PD1_EXHAUSTED_T_CELLS": "PHENO_PD1_T_CELLS",
}

# Home module and layer for shared nodes (overrides whatever a module file assigned).
CANON: dict[str, tuple[str, int]] = {}


def _c(module: str, layer: int, *ids: str) -> None:
    for i in ids:
        CANON[i] = (module, layer)


_c("core_complexes", 4, "MTOR", "MTORC1", "MTORC2", "RPTOR", "RICTOR", "MLST8", "MAPKAP1", "PRR5", "PRR5L", "MTORC3", "ETV7", "FRB_DOMAIN", "ATP_SITE", "KINASE_DOMAIN", "FAT_DOMAIN", "FATC_DOMAIN", "HEAT_REPEATS", "CRIM_DOMAIN", "PH_DOMAIN", "RBD_DOMAIN", "TOS_MOTIF", "RAIP_MOTIF")
_c("core_complexes", 3, "DEPTOR", "AKT1S1")
_c("core_complexes", 2, "TTT_COMPLEX", "TELO2", "TTI1", "TTI2", "HSP90AA1", "CDC37")
_c("core_complexes", 5, "BTRC")
_c("core_complexes", 6, "RIBOSOME")
_c("lysosomal_nutrient_sensing", 3, "RRAGA", "RRAGB", "RRAGC", "RRAGD", "RAG_HETERODIMER", "LYSOSOME")
_c("lysosomal_nutrient_sensing", 2, "RAGULATOR", "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5", "GATOR1", "DEPDC5", "NPRL2", "NPRL3", "GATOR2", "WDR24", "WDR59", "MIOS", "SEH1L", "SEC13", "KICSTOR", "KPTN", "ITFG2", "C12ORF66", "SZT2", "FLCN", "FNIP1", "FNIP2", "FLCN_FNIP", "ARF1", "RAB1A", "RNF152", "SKP2", "SH3BP4", "RUVBL1", "RUVBL2")
_c("lysosomal_nutrient_sensing", 1, "V_ATPASE", "SLC38A9", "NPC1", "SESN1", "SESN2", "SESN3", "CASTOR1", "CASTOR2", "SAMTOR", "LARS1", "SAR1B", "SLC7A5", "SLC3A2", "GCN2", "EIF2AK4")
_c("lysosomal_nutrient_sensing", 0, "LEUCINE", "ARGININE", "METHIONINE", "SAM", "CHOLESTEROL", "AMINO_ACIDS", "GLUTAMINE", "ALPHA_KETOGLUTARATE", "AMINO_ACID_STARVATION", "LYSOSOMAL_AMINO_ACIDS")
_c("growth_factor_input", 0, "INSULIN", "IGF1", "GROWTH_FACTORS", "EGF", "WNT", "TNF", "INFLAMMATION", "MECHANICAL_LOAD", "RESISTANCE_EXERCISE", "ESTROGEN", "ANDROGEN", "LEPTIN")
_c("growth_factor_input", 1, "INSR", "IGF1R", "RTK", "EGFR", "IRS1", "IRS2", "PI3K", "PIK3R1", "PTEN", "PIP2", "PIP3", "PLD1", "PLD2", "PHOSPHATIDIC_ACID", "RAS", "NF1", "LEPR")
_c("growth_factor_input", 2, "PDK1", "AKT", "AKT2", "AKT3", "TSC1", "TSC2", "TBC1D7", "TSC_COMPLEX", "RAF", "MEK", "ERK", "RSK", "GSK3B", "IKBKB", "MAPKAPK2", "MK2", "CK1", "CSNK1A1", "PKA")
_c("growth_factor_input", 3, "RHEB")
_c("energy_stress_input", 0, "HYPOXIA", "DNA_DAMAGE", "ENERGY_STRESS", "AMP", "ATP", "OXYGEN", "GLUCOSE", "GLUCOSE_DEPRIVATION", "EXERCISE", "CALORIC_RESTRICTION", "ER_STRESS", "OXIDATIVE_STRESS", "CALCIUM")
_c("energy_stress_input", 1, "STK11", "STRADA", "CAB39", "LKB1_COMPLEX", "CAMKK2", "HK2_GLUCOSE_SENSING")
_c("energy_stress_input", 2, "AMPK", "DDIT4", "DDIT4L", "TP53", "BNIP3", "PML", "ATM", "GCN2_ATF4")
_c("mtorc1_output_translation_ribosome", 5, "RPS6KB1", "RPS6KB2", "EIF4EBP1", "EIF4EBP2", "EIF4EBP3", "LARP1", "MAF1", "RRN3", "UBTF")
_c("mtorc1_output_translation_ribosome", 6, "RPS6", "EIF4E", "EIF4G1", "EIF4F", "EIF4B", "EIF4A1", "PDCD4", "SKAR", "EEF2K", "EEF2", "TOP_MRNA", "POLR1", "POLR3", "CCND1", "MYC", "VEGFA", "MCL1", "HIF1A", "MNK1", "MKNK1", "CBP80", "NCBP1", "TFIIIC", "PRE_RRNA_47S", "TRNA", "RRNA_5S")
_c("mtorc1_output_translation_ribosome", 7, "CAP_DEPENDENT_TRANSLATION", "RIBOSOME_BIOGENESIS", "PROTEIN_SYNTHESIS", "CELL_GROWTH", "TOP_MRNA_TRANSLATION", "ELONGATION")
_c("mtorc1_output_autophagy_metabolism_feedback", 5, "ULK1", "ATG13", "ATG14", "AMBRA1", "RB1CC1", "TFEB", "TFE3", "LPIN1", "GRB10", "UVRAG", "WIPI2", "PACER", "STAT3", "P300", "EP300", "PPP3", "MITF")
_c("mtorc1_output_autophagy_metabolism_feedback", 4, "MCOLN1")
_c("mtorc1_output_autophagy_metabolism_feedback", 6, "YWHA", "SREBF1", "SREBF2", "G6PD", "CAD", "ATF4", "MTHFD2", "GLUD1", "YY1", "PPARGC1A", "PPARG", "SLC2A1", "HK2", "PFKFB3", "PKM2", "LDHA", "SCAP", "INSIG", "IL1A", "MAPKAPK2_TRANSLATION", "CEBPB")
_c("mtorc1_output_autophagy_metabolism_feedback", 7, "CLEAR_NETWORK", "AUTOPHAGY", "LYSOSOME_BIOGENESIS", "LIPOGENESIS", "CHOLESTEROL_SYNTHESIS", "GLYCOLYSIS", "PENTOSE_PHOSPHATE_PATHWAY", "PYRIMIDINE_SYNTHESIS", "PURINE_SYNTHESIS", "GLUTAMINE_ANAPLEROSIS", "MITOCHONDRIAL_BIOGENESIS", "SASP", "ADIPOGENESIS", "ONE_CARBON_METABOLISM", "NUCLEOTIDE_SYNTHESIS", "AUTOPHAGY_DERIVED_AMINO_ACIDS", "PHAGOPHORE_INITIATION")
_c("mtorc2_regulation_and_outputs", 5, "SGK1", "PKC", "PRKCA", "PRKCB", "PRKCD", "PRKCE", "PRKCZ", "XPLN", "ARHGEF3")
_c("mtorc2_regulation_and_outputs", 6, "FOXO1", "FOXO3", "BAD", "MDM2", "NDRG1", "SCNN1", "ATP1A1", "CLCNK", "TIAM1", "RAC1", "AS160", "TBC1D4", "GLUT4", "SLC2A4", "RHOA", "CDC42", "GLUCOKINASE", "GCK", "GLUCONEOGENESIS")
_c("mtorc2_regulation_and_outputs", 7, "ACTIN_CYTOSKELETON", "CELL_SURVIVAL", "ION_TRANSPORT", "LTP_LATE_PHASE", "CELL_MIGRATION", "EPITHELIAL_POLARITY", "GLUCOSE_UPTAKE")
_c("circadian", 6, "BMAL1", "CLOCK", "CLOCK_BMAL1", "CRY1", "CRY2", "PER1", "VIP")
_c("circadian", 2, "PER2")
_c("circadian", 0, "LIGHT")
_c("circadian", 7, "CIRCADIAN_CLOCK", "RHYTHMIC_TRANSLATION", "RHYTHMIC_AUTOPHAGY")
_c("circadian", 8, "SCN")
_c("pharmacology", 0, "RAPAMYCIN", "EVEROLIMUS", "TEMSIROLIMUS", "RIDAFOROLIMUS", "NAB_SIROLIMUS", "SAPANISERTIB", "VISTUSERTIB", "TORIN1", "TORIN2", "PP242", "KU0063794", "DACTOLISIB", "APITOLISIB", "RAPALINK1", "RMC5552", "RTB101", "ALPELISIB", "CAPIVASERTIB", "TACROLIMUS", "CYCLOSPORINE", "METFORMIN", "AZD8055", "OSI027", "MLN0128", "DEXAMETHASONE", "RAPALOGS", "ATP_COMPETITIVE_INHIBITORS", "BI_STERIC_INHIBITORS")
_c("pharmacology", 1, "FKBP12", "FKBP12_RAPAMYCIN", "ABCB1", "CYP3A4")


def canon_id(raw: str) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", str(raw)).strip("_").upper()
    return ID_ALIASES.get(s, s)


EFFECT_ALIASES = {
    "activation": "activates", "activate": "activates", "recruitment": "activates", "recruits": "activates",
    "stabilises": "activates", "stabilizes": "activates", "stabilisation": "activates", "promotes": "activates",
    "enhances": "activates", "positive": "activates", "inhibition": "inhibits", "inhibit": "inhibits", "represses": "inhibits",
    "suppresses": "inhibits", "sequesters": "inhibits", "negative": "inhibits", "degrades": "inhibits", "degradation": "inhibits",
    "destabilises": "inhibits", "destabilizes": "inhibits", "binds": "modulates",
    "binding": "modulates", "modulate": "modulates", "regulates": "modulates", "neutral": "modulates", "required": "required_for",
    "requires": "required_for", "necessary_for": "required_for", "component": "component_of", "part_of": "component_of",
    "subunit_of": "component_of", "converts": "converts_to", "produces": "converts_to", "generates": "converts_to",
}
INTERACTION_ALIASES = {
    "phosphorylates": "phosphorylation", "phospho": "phosphorylation", "dephosphorylates": "dephosphorylation",
    "binds": "binding", "gap": "gap_activity", "gef": "gef_activity", "recruits": "recruitment", "sequesters": "sequestration",
    "ubiquitylation": "ubiquitination", "stabilization": "stabilisation", "transcription": "transcriptional_regulation",
    "translation": "translational_regulation", "drug_binding": "binding", "allosteric": "allosteric_activation",
    "complex": "complex_formation", "localisation": "membrane_targeting", "localization": "membrane_targeting",
    "o_glcnacylation": "O_GlcNAcylation", "glycosylation": "O_GlcNAcylation", "loss_of_function": "genetic_loss_of_function",
    "gain_of_function": "genetic_gain_of_function", "output": "physiological_output", "phenotype": "physiological_output",
    "metabolic": "metabolic_conversion", "proteolysis": "degradation", "scaffold": "scaffolding", "sense": "sensing",
}
TYPE_ALIASES = {
    "small_molecule": "metabolite", "compound": "drug", "gene": "protein", "enzyme": "protein", "kinase": "protein",
    "gtpase": "protein", "transcription_factor": "protein", "cell_type": "tissue", "organ": "tissue", "pathway": "process",
    "program": "transcription_program", "programme": "transcription_program", "variant": "gene_variant", "mutation": "gene_variant",
    "lesion": "gene_variant", "hormone": "stimulus", "nutrient": "metabolite", "condition": "stimulus", "ion": "metabolite",
    "mrna": "rna", "cellular_process": "process", "clinical_phenotype": "phenotype",
}


def norm_type(t: Any, nid: str) -> str:
    t = str(t or "").strip().lower().replace(" ", "_").replace("-", "_")
    t = TYPE_ALIASES.get(t, t)
    if t in NODE_TYPES:
        return t
    if nid.startswith("VAR_"):
        return "gene_variant"
    if nid.startswith("DISEASE_"):
        return "disease"
    if nid.startswith("PHENO_"):
        return "phenotype"
    if nid.startswith("TISSUE_"):
        return "tissue"
    return "protein"


def default_layer(nid: str, ntype: str) -> int:
    if nid in CANON:
        return CANON[nid][1]
    if nid.startswith(("TISSUE_", "DISEASE_", "PHENO_")) or ntype in ("tissue", "disease", "phenotype"):
        return 8
    if nid.startswith("VAR_") or ntype in ("gene_variant", "drug", "stimulus"):
        return 0
    if ntype in ("metabolite", "lipid"):
        return 1
    if ntype in ("process", "transcription_program"):
        return 7
    return 4


# --------------------------------------------------------------------------- input

def load_data(directory: Path) -> tuple[list[dict], list[dict]]:
    """Read the node and edge files."""
    nodes = json.loads((directory / "nodes.json").read_text(encoding="utf-8"))
    edges = json.loads((directory / "edges.json").read_text(encoding="utf-8"))
    return nodes, edges


# --------------------------------------------------------------------------- normalisation

EVIDENCE_RANK = {"established": 0, "emerging": 1, "contested": 2}


def normalise(nodes_raw: list[dict], edges_raw: list[dict], refs: dict[int, Any]) -> tuple[list[dict], list[dict], dict[str, Any]]:
    log: dict[str, Any] = collections.defaultdict(list)
    renames: dict[str, str] = {}

    # ---- edges: vocabularies, ids, refs
    edges: list[dict] = []
    for e in edges_raw:
        src, tgt = canon_id(e.get("source", "")), canon_id(e.get("target", ""))
        for raw, new in ((e.get("source"), src), (e.get("target"), tgt)):
            if raw and str(raw) != new:
                renames[str(raw)] = new
        if not src or not tgt:
            log["dropped_edges"].append({"edge": e, "why": "empty endpoint"})
            continue
        inter = str(e.get("interaction") or "other").strip().replace(" ", "_").replace("-", "_")
        inter_l = inter.lower()
        if inter not in INTERACTION_TYPES:
            fixed = INTERACTION_ALIASES.get(inter_l, inter_l if inter_l in INTERACTION_TYPES else ("O_GlcNAcylation" if inter_l == "o_glcnacylation" else "other"))
            log["interaction_fixes"].append({"edge": f"{src}->{tgt}", "from": inter, "to": fixed})
            inter = fixed
        eff = str(e.get("effect") or "modulates").strip().lower().replace(" ", "_").replace("-", "_")
        if eff not in EFFECTS:
            fixed = EFFECT_ALIASES.get(eff, "modulates")
            log["effect_fixes"].append({"edge": f"{src}->{tgt}", "from": eff, "to": fixed})
            eff = fixed
        ev = str(e.get("evidence") or "emerging").strip().lower()
        if ev not in EVIDENCE_LEVELS:
            fixed = "contested" if "contest" in ev else ("established" if "establ" in ev else "emerging")
            log["evidence_fixes"].append({"edge": f"{src}->{tgt}", "from": ev, "to": fixed})
            ev = fixed
        direction = str(e.get("direction") or "").strip()
        if direction not in DIRECTIONS:
            fixed = "core"
            if "pharm" in direction.lower() or e.get("module") == "pharmacology":
                fixed = "pharmacological"
            elif e.get("module") == "tissue_physiology_disease":
                fixed = "pathophysiological"
            elif e.get("module") == "circadian":
                fixed = "circadian"
            elif e.get("module") in ("mtorc1_output_translation_ribosome", "mtorc1_output_autophagy_metabolism_feedback"):
                fixed = "downstream_of_mTORC1"
            elif e.get("module") == "mtorc2_regulation_and_outputs":
                fixed = "downstream_of_mTORC2"
            elif e.get("module") in ("lysosomal_nutrient_sensing", "growth_factor_input", "energy_stress_input"):
                fixed = "upstream_of_mTORC1"
            log["direction_fixes"].append({"edge": f"{src}->{tgt}", "from": direction, "to": fixed})
            direction = fixed
        good_refs, bad_refs = [], []
        for r in e.get("refs") or []:
            try:
                ri = int(r)
            except (TypeError, ValueError):
                bad_refs.append(r)
                continue
            (good_refs if ri in refs else bad_refs).append(ri)
        if bad_refs:
            log["bad_refs_dropped"].append({"edge": f"{src}->{tgt}", "refs": bad_refs})
        if not good_refs:
            log["dropped_edges"].append({"edge": f"{src}->{tgt} [{inter}]", "why": "no citation in the bibliography"})
            continue
        if src == tgt and eff in ("component_of",):
            log["dropped_edges"].append({"edge": f"{src}->{tgt} [{inter}]", "why": "self-referential component_of"})
            continue
        module = e.get("module") if e.get("module") in MODULE_KEYS else e.get("_origin")
        if module not in MODULE_KEYS:
            module = "core_complexes"
        sites = sorted({str(s).strip() for s in (e.get("sites") or []) if str(s).strip()})
        edges.append({
            "source": src, "target": tgt, "interaction": inter, "effect": eff, "sites": sites,
            "mechanism": str(e.get("mechanism") or "").strip(), "evidence": ev, "refs": sorted(set(good_refs)),
            "context": (str(e["context"]).strip() if e.get("context") else None), "direction": direction,
            "notes": (str(e["notes"]).strip() if e.get("notes") else None), "module": module,
        })

    # Merge only the same claim; preserve context, sites, polarity and evidence.
    merged: dict[tuple, dict] = {}
    dup = 0
    for e in edges:
        k = (e["source"], e["target"], Edge.from_dict(e).key)
        if k not in merged:
            merged[k] = e
            continue
        dup += 1
        m = merged[k]
        m["refs"] = sorted(set(m["refs"]) | set(e["refs"]))
        if e["notes"] and e["notes"] != m["notes"]:
            m["notes"] = (m["notes"] + " | " + e["notes"]) if m["notes"] else e["notes"]
    edges = list(merged.values())
    log["merged_duplicate_edges"] = dup

    # ---- nodes
    occ: dict[str, list[dict]] = collections.defaultdict(list)
    for n in nodes_raw:
        if not n or not n.get("id"):
            continue
        nid = canon_id(n["id"])
        if str(n["id"]) != nid:
            renames[str(n["id"])] = nid
        occ[nid].append(n)
    touch: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    for e in edges:
        touch[e["source"]][e["module"]] += 1
        touch[e["target"]][e["module"]] += 1
    used = {e["source"] for e in edges} | {e["target"] for e in edges}
    nodes: list[dict] = []
    for nid in sorted(used | set(occ)):
        if nid not in used:
            log["unused_nodes_dropped"].append(nid)
            continue
        cands = occ.get(nid, [])
        if not cands:
            log["nodes_without_record"].append(nid)
        ntype = norm_type(next((c.get("type") for c in cands if c.get("type")), None), nid)
        if nid in CANON:
            module, layer = CANON[nid]
        else:
            module = touch[nid].most_common(1)[0][0] if touch[nid] else next((c.get("module") for c in cands if c.get("module") in MODULE_KEYS), "core_complexes")
            layers = [int(c["layer"]) for c in cands if isinstance(c.get("layer"), (int, float))]
            layer = collections.Counter(layers).most_common(1)[0][0] if layers else default_layer(nid, ntype)
            layer = max(0, min(8, layer))
        home = [c for c in cands if c.get("module") == module] or cands
        label = next((c.get("label") for c in home if c.get("label")), None) or next((c.get("label") for c in cands if c.get("label")), None) or nid.replace("_", " ").title()
        desc = max((str(c.get("description") or "") for c in cands), key=len, default="")
        gene = next((c.get("gene") for c in cands if c.get("gene")), None)
        if isinstance(gene, list):
            gene = "/".join(str(g) for g in gene)
        aliases = sorted({str(a) for c in cands for a in (c.get("aliases") or []) if str(a) and str(a) != label})
        nodes.append({"id": nid, "label": str(label), "type": ntype, "gene": gene, "module": module, "layer": layer,
                      "description": desc, "aliases": aliases})
    log["id_renames"] = dict(sorted(renames.items()))
    return nodes, edges, dict(log)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, default=ROOT / "data")
    args = ap.parse_args()
    refs = {int(r["id"]): r for r in json.loads((args.data / "references.json").read_text(encoding="utf-8"))}
    nodes_raw, edges_raw = load_data(args.data)
    nodes, edges, log = normalise(nodes_raw, edges_raw, refs)
    edges.sort(key=lambda e: (MODULE_ORDER.index(e["module"]), e["source"], e["target"], Edge.from_dict(e).key))
    nodes.sort(key=lambda n: (n["layer"], MODULE_ORDER.index(n["module"]), n["id"]))
    (args.data / "nodes.json").write_text(json.dumps(nodes, ensure_ascii=False, indent=1), encoding="utf-8")
    (args.data / "edges.json").write_text(json.dumps(edges, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"{len(nodes_raw)} nodes, {len(edges_raw)} records in -> {len(nodes)} nodes, {len(edges)} records out")
    for k, v in log.items():
        if v:
            print(f"  {k}: {len(v) if hasattr(v, '__len__') else v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
