#!/usr/bin/env python
"""Rebuild the exact capped-graft SMILES of every scored candidate (CPU only, RDKit).

Each candidate is defined in a `config/*.py` candidate file as a fragment that is grafted
onto that file's 4-methylbenzyl scaffold. This script grafts every definition with the
package's own validator (`redox.structures.contracts.validate_group`, the same call the
library builder uses) and keeps, for each scored candidate, the definition whose average
molecular weight matches the `MW` column of the scored table to within 0.15 g/mol. A
candidate with zero or more than one distinct matching structure is an error.

  PYTHONPATH=src python scripts/manuscript/pool_structures.py
    -> data/manuscript/pool_structures.csv   (id, smiles, config)
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from rdkit import Chem, RDLogger  # noqa: E402
from rdkit.Chem.Descriptors import MolWt  # noqa: E402

from redox.core.common import load_python_config  # noqa: E402
from redox.structures.contracts import validate_group  # noqa: E402

SCORED = ROOT / "data/manuscript/full_pool_scored_legacyG_MeCN_TBAPF6.csv"
OUT = ROOT / "data/manuscript/pool_structures.csv"
MW_TOL = 0.15


def definitions():
    found = {}
    for path in sorted((ROOT / "config").glob("*.py")):
        try:
            module = load_python_config(path)
        except Exception:
            continue
        if not (hasattr(module, "SCAFFOLD") and hasattr(module, "GROUPS")):
            continue
        for group in module.GROUPS:
            if not isinstance(group, dict) or "id" not in group:
                continue
            try:
                smi = Chem.MolToSmiles(validate_group(module.SCAFFOLD, group))
            except Exception:
                continue
            found.setdefault(group["id"], []).append((path.relative_to(ROOT).as_posix(), smi))
    return found


def main() -> int:
    RDLogger.DisableLog("rdApp.*")
    defs = definitions()
    rows, errors = [], []
    for r in csv.DictReader(open(SCORED)):
        match = [(cfg, smi) for cfg, smi in defs.get(r["id"], [])
                 if abs(MolWt(Chem.MolFromSmiles(smi)) - float(r["MW"])) < MW_TOL]
        distinct = {Chem.CanonSmiles(smi) for _, smi in match}
        if len(distinct) != 1:
            errors.append(f"{r['id']}: {len(distinct)} structures match MW {r['MW']}")
            continue
        rows.append(dict(id=r["id"], smiles=distinct.pop(), config=match[0][0]))
    if errors:
        raise SystemExit("\n".join(errors))
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "smiles", "config"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} structures -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
