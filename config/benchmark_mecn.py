"""Experimental MeCN benchmark for E1, E2 and the wave spacing (E1 - E2), built from the
hand-curated literature table data/raw/validation/two_wave_mecn/two_wave_mecn.csv (every row
carries its source, DOI, table/page and original reference electrode; see that folder's
README.md). Entries use the config/validation.py schema so they run through the identical
production pipeline (build -> UMA -> DFT+SMD opt -> active-protocol SP -> thermal).

Which experimental value is usable is decided HERE, per wave, from what the source reports
(no value is corrected or fitted):
  tier "A"  : E deg' or E1/2 of a (quasi)reversible wave           -> accuracy statistics
  tier "B"  : half-peak potential Ep/2 (Schimanofsky 2022 AQs; the source notes Ep/2 is offset
              from E1/2 for a reversible wave)                     -> reported separately
  excluded  : cathodic peak Epc only (not an E deg; wave may be irreversible), unspecified
              stereocenters (CLAUDE.md stereo rule), or inferred compound identity.
Molecules already in the library under another id are mapped to it (REUSE) instead of being
recomputed.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

_CSV = (Path(__file__).resolve().parents[1] / "data" / "raw" / "validation" / "two_wave_mecn"
        / "two_wave_mecn.csv")

# Existing library ids for the same compound (same SMILES) — reuse their calculations, with
# their existing state names: experimental couple label -> that id's event name.
REUSE = {
    "9,10-anthraquinone": ("anthraquinone_parent", {"0/-1": "neu->red1", "-1/-2": "red1->red2"}),
    "methyl viologen (1,1'-dimethyl-4,4'-bipyridinium)":
        ("methyl_viologen", {"+2/+1": "ox2->ox1", "+1/0": "ox1->neu"}),
    "2,2,6,6-tetramethylpiperidine-1-oxyl (TEMPO)": ("tempo_parent", {"+1/0": "ox->rad"}),
    "10H-phenothiazine": ("phenothiazine_parent", {"+1/0": "ox->neu"}),
}

# Rows excluded entirely, with the reason (kept visible in the benchmark output).
EXCLUDE = {
    "N,N'-bis(2-ethylhexyl)pyromellitic diimide":
        "2 unspecified stereocenters (2-ethylhexyl); CLAUDE.md requires fully specified stereo",
    "N,N'-bis(2-ethylhexyl)-1,4,5,8-naphthalene diimide":
        "compound identity inferred by the curator (source labels it only 'NDI'); stereocenters",
}

_EPC_ONLY = ("5-hydroxy-1,4-naphthoquinone", "5,8-dihydroxy-1,4-naphthoquinone")
_EP2 = ("9,10-anthraquinone", "1-hydroxy-9,10-anthraquinone", "1,4-dihydroxy-9,10-anthraquinone",
        "1-amino-9,10-anthraquinone", "1,4-diamino-9,10-anthraquinone",
        "1-amino-4-hydroxy-9,10-anthraquinone")


def _ref_is_primary(conversion: str) -> bool:
    """Absolute value is grounded only if reported directly vs Fc, or converted with the
    authors' OWN measured calibration. A literature conversion constant we could not verify
    (Pavlishchuk & Addison via a secondary transcription) is not."""
    c = conversion.lower()
    return c.startswith("none") or "authors' own" in c


def _tier(r: dict, wave: int) -> tuple[str, str]:
    name = r["name"].strip()
    if name.startswith(_EPC_ONLY):
        return "excluded", "cathodic peak potential at 1000 mV/s (not an E1/2)"
    if wave == 2 and r.get("wave2_reversible", "").strip().lower() != "true":
        return "excluded", "second wave not reported as (quasi)reversible with an E1/2"
    if name.startswith(_EP2):
        return "B", "half-peak potential Ep/2 (source: Ep/2 offset from E1/2)"
    if "h2o" in r["electrolyte"].lower() or "water" in r["electrolyte"].lower():
        return "B", "electrolyte contains added water (not dry MeCN)"
    if r["confidence"].strip().lower() == "low":
        return "B", "curator confidence low (secondary / extrapolated value)"
    if not _ref_is_primary(r["conversion"]):
        return "B", "absolute value relies on an unverified literature conversion constant"
    return "A", "E deg' / E1/2 vs Fc (direct or authors' own calibration)"


def _slug(name: str) -> str:
    s = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
    return "bm_" + s


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _n_electrons(smiles: str) -> int:
    from rdkit import Chem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return sum(a.GetAtomicNum() for a in m.GetAtoms()) - Chem.GetFormalCharge(m)


def _layout(r: dict):
    """(states, [(couple, event, value, wave)]) for a new (non-reused) molecule."""
    smi, couple = r["smiles_neutral"], (r.get("couple") or "").strip()
    e1, e2 = _f(r["E1_vs_Fc_V"]), _f(r["E2_vs_Fc_V"])
    if "[n+]" in smi:      # viologen: SMILES is the dication; waves 2+/1+ and 1+/0
        st = [("ox2", 2, 1, 0), ("ox1", 1, 2, -1)] + ([("neu", 0, 1, -2)] if e2 is not None else [])
        return st, [("+2/+1", "ox2->ox1", e1, 1), ("+1/0", "ox1->neu", e2, 2)]
    if couple.startswith("+1/0"):   # oxidation of a neutral; multiplicities from e- parity
        odd = _n_electrons(smi) % 2
        st = [("neu", 0, 2 if odd else 1, 0), ("ox", 1, 1 if odd else 2, +1)]
        return st, [("+1/0", "ox->neu", e1, 1)]
    # neutral quinone / imide: 0/-1 and -1/-2 (dianion assumed closed-shell singlet)
    st = [("neu", 0, 1, 0), ("red1", -1, 2, -1)] + ([("red2", -2, 1, -2)] if e2 is not None else [])
    return st, [("0/-1", "neu->red1", e1, 1), ("-1/-2", "red1->red2", e2, 2)]


def load() -> list[dict]:
    """One entry per MOLECULE (rows from several sources are merged): id, name, smiles, states,
    events (one per source x wave, each with tier + reason + citation) and reference-free
    wave-spacing events (E1 - E2) where both waves of one source are usable."""
    by_name: dict[str, dict] = {}
    for r in csv.DictReader(_CSV.open()):
        name = r["name"].strip()
        reuse = REUSE.get(name)
        states, waves = _layout(r)
        ent = by_name.setdefault(name, dict(
            id=reuse[0] if reuse else _slug(name), name=name, smiles=r["smiles_neutral"],
            family="validation", bench_family=r["family"], states=states, events=[],
            reuse=bool(reuse), excluded=EXCLUDE.get(name)))
        cite = dict(source=r["source"], doi=r["doi"], where=r["page_or_table"],
                    confidence=r["confidence"])
        tiers = {}
        for couple, event, val, wave in waves:
            if val is None:
                continue
            if reuse:
                event = reuse[1][couple]
            tier, why = _tier(r, wave)
            tiers[wave] = tier
            ent["events"].append(dict(event=event, kind="E", exp_V_vs_Fc=val, tier=tier,
                                      tier_reason=why, **cite))
        if 1 in tiers and 2 in tiers and "excluded" not in tiers.values():
            ev1 = [e for e in ent["events"] if e["kind"] == "E" and e["source"] == r["source"]]
            sp_tier = "B" if "half-peak" in ev1[0]["tier_reason"] else "A"
            ent["events"].append(dict(event="dE12", kind="spacing",
                                      exp_V=_f(r["E1_vs_Fc_V"]) - _f(r["E2_vs_Fc_V"]),
                                      tier=sp_tier, tier_reason="E1 - E2, reference-free", **cite))
    return list(by_name.values())


BENCHMARK = load()
