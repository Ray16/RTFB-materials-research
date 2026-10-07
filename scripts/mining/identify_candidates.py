#!/usr/bin/env python
"""Systematic candidate identification over the D3TaLES dump (redox.screening.discovery).

Writes results/discovery/:
  rediscovery_check.csv     the rules re-generate every registered grafted candidate
  current_candidates.csv    registered candidates with descriptors computed by the SAME code
  d3tales_pool_status_counts.csv  molecules per (family, rule outcome) — every outcome counted
  d3tales_pool.csv          every in-scope D3TaLES molecule x graft site, with the
                            rule outcome (`status`), nominal capacity, grafted SA, and the
                            D3TaLES-reported values (THEIR level of theory — pre-filter only)
  d3tales_prefilter_front.csv  status=='ok' rows that are non-dominated, per family, on
                            (capacity max, SA min, D3TaLES electron lambda min), using only
                            physical D3TaLES lambda (0-1.5 eV). No scalar score.

  PYTHONPATH=src python scripts/mining/identify_candidates.py [--pubchem]
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from redox.core.common import load_config  # noqa: E402
from redox.screening.discovery import TETHER, describe, identify, rediscover  # noqa: E402

OUT = ROOT / "results" / "discovery"
D3 = ROOT / "data" / "raw" / "validation" / "D3TaLES" / "d3tales_public.csv"
CONFIGS = ("starting_candidates", "merrifield_multielectron", "redox_groups",
           "discovered_candidates")


def _front(df, cols_max=("capacity_nominal_mAh_g",), cols_min=("SA_grafted", "d3_electron_reorg_eV")):
    """Non-dominated rows (all objectives present) — no weights, no scalar."""
    d = df.dropna(subset=list(cols_max) + list(cols_min)).copy()
    v = [(-d[c]).values for c in cols_max] + [d[c].values for c in cols_min]
    import numpy as np
    X = np.vstack(v).T
    keep = []
    for i in range(len(X)):
        dom = ((X <= X[i]).all(1) & (X < X[i]).any(1)).any()
        keep.append(not dom)
    return d[keep]


def pubchem_names(smiles):
    """[(CID, IUPAC name, title)] from PubChem PUG-REST by structure (public SMILES only; no
    personal data sent). Missing -> (None, None, 'not found')."""
    import json
    import time
    import urllib.parse
    import urllib.request
    out = []
    for s in smiles:
        url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
               + urllib.parse.quote(s, safe="") + "/property/IUPACName,Title/JSON")
        try:
            d = json.load(urllib.request.urlopen(url, timeout=30))["PropertyTable"]["Properties"][0]
            out.append((d.get("CID"), d.get("IUPACName"), d.get("Title")))
        except Exception as exc:
            out.append((None, None, f"not found ({type(exc).__name__})"))
        time.sleep(0.3)          # PubChem asks for <= 5 requests/s
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--pubchem", action="store_true",
                    help="name the pre-filter front via PubChem (sends parent SMILES only)")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    F = load_config("electrolyte").FARADAY
    groups = [dict(g, config=c) for c in CONFIGS for g in load_config(c).GROUPS if "frag" in g]
    sc = {r["id"] for r in csv.DictReader(open(ROOT / "results" / "scorecard.csv"))}

    redisc = pd.DataFrame(rediscover(groups, F))
    redisc.to_csv(OUT / "rediscovery_check.csv", index=False)
    print(f"rediscovered {int(redisc.rediscovered.sum())}/{len(redisc)} registered candidates")

    cur = []
    for g, r in zip(groups, redisc.to_dict("records")):
        n_e = 2 if r["multi_e"] else 1
        cur.append(dict(id=g["id"], name=g["name"], config=g["config"], families=r["families"],
                        multi_e=r["multi_e"], grafted_smiles=r["target"], in_scorecard=g["id"] in sc,
                        **describe(r["target"], n_e, F)))
    cur = pd.DataFrame(cur)
    cur.to_csv(OUT / "current_candidates.csv", index=False)

    known = set()
    for r in csv.DictReader(open(ROOT / "library" / "manifest.csv")):
        m = Chem.MolFromSmiles(r["smiles"])
        if m:
            known.add(Chem.MolToInchiKey(m)[:14])
    d3 = pd.read_csv(D3, low_memory=False,
                     usecols=["_id", "smiles", "groundState_charge", "sa_score",
                              "electron_reorganization_energy", "hole_reorganization_energy",
                              "solv_reduction_potential", "solv_oxidation_potential"])
    d3 = d3[d3.groundState_charge == 0]
    recs = [dict(source="D3TaLES", source_id=r["_id"], smiles=r["smiles"],
                 d3_sa_score=r["sa_score"], d3_electron_reorg_eV=r["electron_reorganization_energy"],
                 d3_hole_reorg_eV=r["hole_reorganization_energy"],
                 d3_solv_reduction_potential_raw=r["solv_reduction_potential"],
                 d3_solv_oxidation_potential_raw=r["solv_oxidation_potential"])
            for r in d3.to_dict("records")]
    pool = pd.DataFrame(identify(recs, F, frozenset(known)))
    pool.to_csv(OUT / "d3tales_pool.csv", index=False)
    print(f"{pool.source_id.nunique()} D3TaLES molecules in a family; {len(pool)} molecule x site rows")
    print(pool.groupby(["families", "status"]).source_id.nunique().unstack(fill_value=0).to_string())

    # D3TaLES lambda is a pre-filter only, and only where physical: a 4-point lambda is >= 0
    # by construction (redox.properties.reorg), and D3TaLES's >1.5 eV tail is dominated by
    # failed optimizations (FINDINGS #10). Outside [0, 1.5] eV it is treated as unavailable.
    lam = pool.d3_electron_reorg_eV
    pool["d3_lambda_valid"] = lam.between(0.0, 1.5)
    # every outcome is counted; the per-row table omits only the out-of-scope ring-N grafts
    # (N-heterocycles that would become 1e- aziniums: ~13k rows) to keep the artifact small
    pool.groupby(["families", "status"]).source_id.nunique().rename("n_molecules"
        ).reset_index().to_csv(OUT / "d3tales_pool_status_counts.csv", index=False)
    pool[pool.status != "ring-N graft not forming a viologen"].to_csv(
        OUT / "d3tales_pool.csv", index=False)
    ok = pool[(pool.status == "ok") & ~pool.parent_in_library.astype(bool)]
    print(f"status ok, new: {ok.source_id.nunique()} molecules; D3TaLES lambda usable for "
          f"{ok[ok.d3_lambda_valid].source_id.nunique()} (rest outside [0, 1.5] eV or missing)")
    okv = ok[ok.d3_lambda_valid]
    fronts = pd.concat([_front(g).assign(front_family=f) for f, g in okv.groupby("families")]
                       ) if len(okv) else okv
    if a.pubchem and len(fronts):
        nm = pubchem_names(list(fronts.parent_smiles))
        fronts = fronts.assign(pubchem_cid=[n[0] for n in nm], iupac_name=[n[1] for n in nm],
                               pubchem_title=[n[2] for n in nm])
    fronts.to_csv(OUT / "d3tales_prefilter_front.csv", index=False)
    print(f"pre-filter front (per family, non-dominated on capacity/SA/D3TaLES lambda): "
          f"{len(fronts)} rows, {fronts.source_id.nunique() if len(fronts) else 0} molecules")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
