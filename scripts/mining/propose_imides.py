#!/usr/bin/env python
"""Enumerate graftable PMDI/NDI designs; source-inspired proposals, not property predictions."""
import argparse
import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from rdkit import Chem
from redox.build.contracts import validate_group
from redox.core.common import load_config, read_manifest
from redox.screening.discovery import describe, precursor_from_frag, identify

SCAFFOLD = "Cc1ccc(C[*:1])cc1"
STATES = [("neu", 0, 1, 0), ("red1", -1, 2, -1), ("red2", -2, 1, -2)]
CORES = {
    "pmdi": ("[*:1]n1c(=O)c2cc3c(=O)n({cap})c(=O)c3cc2c1=O",
             "pyromellitic diimide", "https://doi.org/10.1021/acsaem.1c01490"),
    "ndi": ("[*:1]N1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N({cap})C3=O",
            "naphthalene diimide", "https://doi.org/10.1002/cssc.202501011"),
}
CAPS = {
    "methyl": ("C", "methyl"), "ethyl": ("CC", "ethyl"),
    "methoxyethyl": ("CCOC", "(2-methoxyethyl)"),
    "methoxyethoxyethyl": ("CCOCCOC", "[2-(2-methoxyethoxy)ethyl]"),
    "isopropyl": ("C(C)C", "isopropyl"),
}


def proposals():
    known = {Chem.MolToInchiKey(Chem.MolFromSmiles(r["smiles"])) for r in read_manifest()}
    faraday = load_config("electrolyte").FARADAY
    groups, rows = [], []
    for core, (template, core_name, source) in CORES.items():
        for capid, (cap, cap_name) in CAPS.items():
            frag = template.format(cap=cap)
            group = dict(id=f"design_{core}_{capid}",
                         name=f"N-(4-methylbenzyl)-N'-{cap_name} {core_name}",
                         family="imide (n-type)", frag=frag, states=STATES,
                         provenance=dict(source=source, type="source-inspired grafted design",
                                         measured_for_this_structure=False))
            mol = validate_group(SCAFFOLD, group)
            if Chem.MolToInchiKey(mol) in known:
                continue
            smi = Chem.MolToSmiles(mol)
            precursor = precursor_from_frag(frag)
            identified = identify([dict(source_id=group["id"], smiles=precursor)], faraday)
            if not any(r.get("grafted_smiles") == smi and r.get("status") == "ok" for r in identified):
                raise ValueError(f"{group['id']}: grafting rules do not regenerate the proposed bis-imide")
            groups.append(group)
            rows.append(dict(id=group["id"], name=group["name"], precursor_smiles=precursor,
                             grafted_smiles=smi, frag=frag, source=source,
                             status="proposed; no production QM", **describe(smi, 2, faraday)))
    return groups, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    if a.out.exists():
        ap.error("choose a new directory; prior proposals are preserved")
    groups, rows = proposals()
    a.out.mkdir(parents=True)
    (a.out / "proposals.json").write_text(json.dumps(groups, indent=2) + "\n")
    with (a.out / "imide_proposals.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    (a.out / "proposals.py").write_text(
        '"""Uncomputed design proposals; validate before adopting into custom_candidates.py."""\n'
        + f"SCAFFOLD = {SCAFFOLD!r}\nGROUPS = {groups!r}\n")
    print(f"{len(groups)} new graftable imide proposals; no new QM submitted. {a.out}")


if __name__ == "__main__":
    main()
