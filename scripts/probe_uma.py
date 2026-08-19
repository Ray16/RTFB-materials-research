"""Probe the installed fairchem/UMA API: which model names load, and how charge+spin
are passed. Runs single points on the pyridinium ox (cation, singlet) and red (neutral
radical, doublet) seeds to confirm charge/spin conditioning changes the energy."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

from fairchem.core import pretrained_mlip, FAIRChemCalculator
from ase.io import read

# 1. find a UMA model this fairchem version can load
CANDIDATES = ["uma-s-1p2p1", "uma-s-1p2", "uma-s-1p1", "uma-s-1"]
try:
    avail = pretrained_mlip.available_models
    print("available_models:", avail)
except Exception as e:
    print("available_models: (n/a)", e)

model = None
predictor = None
for name in CANDIDATES:
    try:
        predictor = pretrained_mlip.get_predict_unit(name, device="cuda")
        model = name
        print(f"LOADED model: {name}")
        break
    except Exception as e:
        print(f"  cannot load {name}: {type(e).__name__}: {str(e)[:160]}")
if predictor is None:
    sys.exit("no UMA model could be loaded")

calc = FAIRChemCalculator(predictor, task_name="omol")

# 2. single points on pyridinium ox (q=+1,m=1) and red (q=0,m=2)
for state, q, m in [("ox", 1, 1), ("red", 0, 2)]:
    atoms = read(str(ROOT / "library" / "pyridinium" / f"{state}.xyz"))
    atoms.info["charge"] = q
    atoms.info["spin"] = m
    atoms.calc = calc
    e = atoms.get_potential_energy()
    print(f"pyridinium/{state}: q={q} m={m}  E={e:.6f} eV  natoms={len(atoms)}")

print(f"\nMODEL_USED={model}")
