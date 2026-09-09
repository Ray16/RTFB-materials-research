"""Electrolyte / solvation and electrochemical-referencing parameters.

Facts live in config/project.json (the machine-readable, documented single source of truth);
this module LOADS them and exposes the names the pipeline imports (SOLVENT, REFERENCE,
WINDOW_V_VS_FC, FARADAY, ...). Edit values in project.json, not here.

The DFT+SMD step (Stage 2) reads SOLVENT; the redox-potential calculation reads REFERENCE to
convert absolute potentials to a reported scale. UMA/OMol (Stage 1) is gas-phase and does NOT
use these — solvation is applied only at the DFT single point.
"""
import json
from pathlib import Path

_FACTS = json.loads((Path(__file__).with_name("project.json")).read_text())

# --- Implicit solvation (Stage-2 DFT) + solvent dielectrics (outer-sphere reorg) ---
_s = _FACTS["solvent"]
SOLVENT = dict(
    name=_s["name"],
    model=_s["model"],            # PySCF: mf = mf.SMD(); mf.with_solvent.solvent = name
    eps=_s["eps_r"],              # static permittivity (back-compat key)
    abbrev=_s["abbrev"],
    eps_r=_s["eps_r"],            # static (relative) permittivity
    eps_optical=_s["eps_optical"],  # optical permittivity = refractive_index**2
    refractive_index=_s["refractive_index"],
    solvent_probe_radius_A=_s["solvent_probe_radius_A"],  # SASA probe for r_hyd estimate
)

# Supporting-electrolyte counterion (context; enters via ion-pairing / activity, not the
# continuum by default). Model explicitly only if ion-pairing is being studied.
COUNTERION = _FACTS["counterion"]

# --- Electrochemical stability window (V vs Fc/Fc+) + anolyte/catholyte divider ---
WINDOW_V_VS_FC = tuple(_FACTS["electrochemical_window_V_vs_Fc"])
ANOLYTE_CATHOLYTE_DIVIDER_V = _FACTS["anolyte_catholyte_divider_V"]

# --- Electrochemical referencing (absolute -> reported scale) ---
# E_abs = -dG / (n F).  E_vs_ref = E_abs - E_ref_abs.
FARADAY = _FACTS["constants"]["faraday_C_per_mol"]        # C/mol
_r = _FACTS["reference"]
REFERENCE = dict(
    she_abs_V=_r["she_abs_V"],
    fc_vs_she_V=_r["fc_vs_she_V"],
    default_scale=_r["default_scale"],
)

# Level-matched ferrocene reference (PREFERRED): compute Fc/Fc+ at the SAME level and subtract
# it, cancelling the systematic error in the absolute reference. redox.py prefers this over the
# thermodynamic she_abs + fc_vs_she constant.
FC_ABS_COMPUTED_V = _r["fc_abs_computed_V"]
# Published level-matched Fc references in MeCN (OROP SI) — cross-check for our computed value.
FC_ABS_REF_MeCN_V = dict(_r["fc_abs_ref_MeCN_V"])
