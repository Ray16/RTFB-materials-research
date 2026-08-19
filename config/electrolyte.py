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
