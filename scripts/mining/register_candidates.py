#!/usr/bin/env python
"""Append newly identified leads to the growing candidate registry
(config/discovered_candidates.py). APPEND-ONLY: entries already registered — in ANY candidate
config, matched by grafted-structure InChIKey — are never added twice or edited; existing
entries are never rewritten. Each new entry records its provenance and the run protocol.

  PYTHONPATH=src python scripts/mining/register_candidates.py \
      --from results/discovery/d3tales_prefilter_front.csv --names names.csv [--dry-run]

`--names` (CSV: source_id,name,flag) supplies the grafted species' real name — derived from
the PubChem parent name, never invented shorthand — and an optional chemistry flag.
"""
from __future__ import annotations

import argparse
import datetime
import pprint
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from redox.core.common import load_config  # noqa: E402
from redox.screening.discovery import TETHER  # noqa: E402

REG = ROOT / "config" / "discovered_candidates.py"
CONFIGS = ("starting_candidates", "merrifield_multielectron", "redox_groups",
           "discovered_candidates")
_TETHER_CORE = Chem.MolFromSmarts("[CH3]c1ccc([CH2])cc1")
QUINONE_STATES = [("neu", 0, 1, 0), ("red1", -1, 2, -1), ("red2", -2, 1, -2)]


def frag_from_grafted(smi: str) -> str:
    """Grafted SMILES -> config fragment: the 4-methylbenzyl tether becomes [*:1]."""
    m = Chem.MolFromSmiles(smi)
    hits = m.GetSubstructMatches(_TETHER_CORE)
    if len(hits) != 1:
        raise ValueError(f"expected exactly one tether in {smi}, found {len(hits)}")
    tether = set(hits[0])
    ch2 = hits[0][5]                     # pattern atom order: CH3, c x4, CH2, c x2
    anchor = next(n.GetIdx() for n in m.GetAtomWithIdx(ch2).GetNeighbors()
                  if n.GetIdx() not in tether)
    rw = Chem.RWMol(m)
    d = rw.AddAtom(Chem.Atom(0)); rw.GetAtomWithIdx(d).SetAtomMapNum(1)
    rw.AddBond(anchor, d, Chem.BondType.SINGLE)
    for i in sorted(tether, reverse=True):
        rw.RemoveAtom(i)
    frag = Chem.MolToSmiles(rw.GetMol())
    back = Chem.MolToSmiles(Chem.molzip(Chem.MolFromSmiles(TETHER + "." + frag)))
    if back != Chem.MolToSmiles(m):
        raise ValueError(f"fragment round-trip failed: {smi} -> {frag} -> {back}")
    return frag


def registered_keys() -> set:
    keys = set()
    for c in CONFIGS:
        try:
            groups = load_config(c).GROUPS
        except (FileNotFoundError, AttributeError):
            continue
        for g in groups:
            if "frag" in g:
                m = Chem.molzip(Chem.MolFromSmiles(TETHER + "." + g["frag"]))
                keys.add(Chem.MolToInchiKey(m)[:14])
    return keys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="src", required=True)
    ap.add_argument("--names", required=True)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    leads = pd.read_csv(a.src)
    names = pd.read_csv(a.names).set_index("source_id")
    have = registered_keys()
    new = []
    for r in leads.to_dict("records"):
        key = Chem.MolToInchiKey(Chem.MolFromSmiles(r["grafted_smiles"]))[:14]
        if key in have:
            print(f"[skip] {r['source_id']}: already registered")
            continue
        if r["source_id"] not in names.index:
            raise SystemExit(f"{r['source_id']}: no name in {a.names} (real names required)")
        if not str(r["families"]).endswith("quinone"):
            raise SystemExit(f"{r['source_id']}: states are only defined here for quinones")
        nm = names.loc[r["source_id"]]
        g = dict(id=f"d3_{r['source_id'].lower()}", name=nm["name"], family="quinone (n-type)",
                 frag=frag_from_grafted(r["grafted_smiles"]), states=QUINONE_STATES,
                 provenance=dict(source=r["source"], source_id=r["source_id"],
                                 parent_smiles=r["parent_smiles"],
                                 pubchem_cid=None if pd.isna(r.get("pubchem_cid")) else int(r["pubchem_cid"]),
                                 pubchem_title=r.get("pubchem_title"), graft=r["graft"],
                                 discovery_rule="redox.screening.discovery",
                                 registered=datetime.date.today().isoformat()))
        if isinstance(nm.get("flag"), str) and nm["flag"].strip():
            g["flag"] = nm["flag"].strip()
        new.append(g); have.add(key)
    print(f"{len(new)} new candidate(s)")
    if a.dry_run or not new:
        for g in new:
            print(f"  {g['id']}: {g['name']} | {g['frag']}")
        return
    text = REG.read_text()
    block = "".join("    " + pprint.pformat(g, width=96, sort_dicts=False).replace("\n", "\n    ")
                    + ",\n" for g in new)
    marker = "# >>> append new entries above this line (register_candidates.py) <<<"
    if marker not in text:
        raise SystemExit(f"registry marker missing in {REG}")
    REG.write_text(text.replace(marker, block + marker))
    print(f"appended to {REG}")


if __name__ == "__main__":
    main()
