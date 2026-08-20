"""Electrolyte / solvation and electrochemical-referencing parameters.

Single source of truth for the environment. The DFT+SMD step (Stage 2) reads SOLVENT;
the redox-potential calculation reads REFERENCE to convert absolute potentials to a
reported scale. UMA/OMol (Stage 1) is gas-phase and does NOT use these — solvation is
applied only at the DFT single-point.
"""

# --- Implicit solvation (Stage-2 DFT single point) ---
SOLVENT = dict(
    name="acetonitrile",
    model="SMD",          # PySCF: mf = mf.SMD(); mf.with_solvent.solvent = "acetonitrile"
    eps=37.5,             # static dielectric constant
    abbrev="MeCN",
)

# Supporting electrolyte counterion (context; enters via ion-pairing / activity, not the
# continuum by default). Model explicitly only if ion-pairing is being studied.
COUNTERION = "PF6-"

# --- Electrochemical referencing (absolute -> reported scale) ---
# E_abs = -dG / (n F).  E_vs_ref = E_abs - E_ref_abs.
FARADAY = 96485.33212            # C/mol
REFERENCE = dict(
    # Absolute potential of SHE (IUPAC ~4.44 V; common computational value 4.28-4.44 V).
    she_abs_V=4.44,
    # Fc/Fc+ vs SHE in acetonitrile (~+0.40 V); internal reference commonly used in MeCN.
    fc_vs_she_V=0.40,
    default_scale="Fc/Fc+",      # report potentials vs ferrocene by default
)

# --- Level-matched ferrocene reference (PREFERRED, physics not fitting) ---
# The rigorous way to put computed potentials on the Fc/Fc+ scale is to compute the
# ferrocene reference at the SAME functional/basis/solvent and subtract it:
#     E_vs_Fc = E_abs(analyte) - E_abs(Fc/Fc+)_same_level
# This cancels the systematic error in the absolute reference that otherwise dominates
# implicit-solvation redox potentials (the OROP 313-system MeCN benchmark shows raw
# implicit DFT carries a ~+0.3 V systematic bias, MAE ~0.5 V; referencing to Fc at the
# same level removes most of the constant part). We compute our own value in the
# validation run; until then FC_ABS_COMPUTED_V is None and redox.py falls back to the
# thermodynamic she_abs + fc_vs_she constant above.
FC_ABS_COMPUTED_V = 4.434749378960987         # filled from calcs/dft/ferrocene once computed at our level

# Published level-matched ferrocene reference potentials in MeCN (OROP SI,
# raw_ferrocene-ref-values.txt) — for cross-checking our own computed value.
FC_ABS_REF_MeCN_V = {
    "b3lyp":     4.63210987710688,   # b3lyp (no D3) / 6-31G*
    "b3lyp-d3":  4.66224854035885,   # b3lyp-D3 / 6-31G*
    "wb97x-d3":  4.618265132412174,  # wB97X-D3 / 6-31G*
}
