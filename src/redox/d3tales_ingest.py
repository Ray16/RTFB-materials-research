"""D3TaLES potential-column HYGIENE — keep D3TaLES's own DFT potentials off our vs-Fc axis.

WHY: The D3TaLES public dump reports `solv_oxidation_potential` / `solv_reduction_potential`
on an ABSOLUTE scale, NOT vs ferrocene. From their source (RedoxPotentialCalc):

    potential = -dG_solv / n + 4.42          # 4.42 V = absolute SHE reference, hard-coded

So e.g. a benzoquinone reduction shows up as ~+8.3 "V" — which is NOT a physical
electrochemical reduction potential on the usual scale, it is that absolute quantity. Left
unlabeled next to our own `E_vs_Fc_V` (redox.py) it invites a category error in candidate
selection. This module RELABELS the raw columns and provides a clearly-marked, approximate
vs-Fc ESTIMATE so the two never get confused.

The two columns are on DIFFERENT effective references (empirically confirmed against
anchors with known MeCN potentials vs Fc/Fc+):
  - oxidation:  E_vs_Fc ~= D3_ox  - OX_OFFSET   (anchor: TEMPO 0.856 -> +0.24)
  - reduction:  E_vs_Fc ~= D3_red - RED_OFFSET  (anchor: anthraquinone 7.832 -> -1.28)
RED_OFFSET is also consistent with the source formula (strip the +4.42 SHE constant, then
subtract an Fc absolute ~4.4-4.8 V => ~8.9-9.2). These offsets are calibrated on only a
handful of in-D3TaLES anchors, so the estimate is good to ~0.3 V at best (the implicit-DFT
floor anyway) and is for RANKING/placement, never a precise absolute number. Add anchors and
refit (calibrate_offsets) to tighten them.

Usage:
    from redox.d3tales_ingest import relabel_potentials, d3tales_red_to_vs_fc
    df = relabel_potentials(df)   # renames raw cols -> *_D3TaLES_abs_V, adds *_vs_Fc_est_V
"""
from __future__ import annotations

# Anchor-calibrated offsets (V). E_vs_Fc_est = D3TaLES_value - OFFSET.
# Refine with calibrate_offsets() as more overlap anchors become available.
OX_OFFSET = 0.62    # oxidation column  (n=1 anchor: TEMPO; ~ Fc-vs-SHE-ish)
RED_OFFSET = 9.11   # reduction column  (n=1 anchor: anthraquinone; ~ 2*SHE + Fc, see source)

# Raw D3TaLES column -> (absolute-relabel, vs-Fc-estimate) column names.
_POTENTIAL_COLS = {
    "solv_oxidation_potential": ("solv_oxidation_potential_D3TaLES_abs_V",
                                 "E_ox_vs_Fc_est_V", OX_OFFSET),
    "oxidation_potential":      ("oxidation_potential_D3TaLES_abs_V",
                                 "E_ox_vs_Fc_est_V", OX_OFFSET),
    "solv_reduction_potential": ("solv_reduction_potential_D3TaLES_abs_V",
                                 "E_red_vs_Fc_est_V", RED_OFFSET),
    "reduction_potential":      ("reduction_potential_D3TaLES_abs_V",
                                 "E_red_vs_Fc_est_V", RED_OFFSET),
}


def d3tales_ox_to_vs_fc(v):
    """D3TaLES oxidation potential (absolute) -> approximate V vs Fc/Fc+ (None-safe)."""
    try:
        return round(float(v) - OX_OFFSET, 3)
    except (TypeError, ValueError):
        return None


def d3tales_red_to_vs_fc(v):
    """D3TaLES reduction potential (absolute) -> approximate V vs Fc/Fc+ (None-safe)."""
    try:
        return round(float(v) - RED_OFFSET, 3)
    except (TypeError, ValueError):
        return None


def relabel_potentials(df):
    """Return a copy of `df` with raw D3TaLES potential columns renamed to explicit
    *_D3TaLES_abs_V and an approximate *_vs_Fc_est_V column added alongside. Idempotent-ish:
    only acts on columns present. Keeps every other column untouched."""
    df = df.copy()
    for raw, (abs_name, est_name, offset) in _POTENTIAL_COLS.items():
        if raw not in df.columns:
            continue
        import pandas as pd
        vals = pd.to_numeric(df[raw], errors="coerce")
        df[est_name] = (vals - offset).round(3)
        df = df.rename(columns={raw: abs_name})
    return df


def calibrate_offsets(anchors):
    """Given anchors = list of dicts {d3tales_value, exp_vs_Fc, column: 'ox'|'red'}, return
    the mean empirical offset per column. Utility for refitting OX_OFFSET/RED_OFFSET when more
    overlap molecules (D3TaLES ∩ known-vs-Fc) are available. Does not mutate the constants."""
    out = {}
    for col in ("ox", "red"):
        diffs = [a["d3tales_value"] - a["exp_vs_Fc"]
                 for a in anchors if a.get("column") == col]
        if diffs:
            out[col] = round(sum(diffs) / len(diffs), 3)
    return out
