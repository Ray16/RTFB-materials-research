#!/usr/bin/env python
"""Per-couple catholyte anchor: compute the EXACT BMEPZ parent (5,10-bis(2-methoxyethyl)-
5,10-dihydrophenazine) at the same DFT+SMD / legacy-xTB-RRHO level as the pool, and compare
its two oxidation potentials to the measured -0.29 / +0.50 V vs Fc (ScienceDirect S2451929419303158).
This turns the catholyte bias from 'corroborated' into a measured per-couple error for the
EXACT same chemistry as the proposed DHP catholytes.

Working geometries go to LOCAL scratch (TMPDIR). The result is persisted as a record
(manuscript/evidence/<date>_gpu_checks/bmepz_anchor.json + the three optimized geometries):
per-state energies, thermal corrections, protocol, SCF diagnostics and timings, the Fc
reference and the sha256 of the ferrocene records it came from, E1/E2 and their errors.
Run on a GPU via gpu_reserve.  Measured couples from the BMEPZ parent paper
(Kwon et al., Chem 2019, 5, 2642; DOI 10.1016/j.chempr.2019.07.006).
"""
import hashlib, json, os, shutil, socket, sys, tempfile, time
from pathlib import Path
os.environ.setdefault('XTB_BIN', '/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/xtb')
from ase.io import read as aseread, write as asewrite
from rdkit import Chem
from rdkit.Chem import AllChem
from redox.backends.dft import dft_smd, _thermal_correction
from redox.core.common import read_result, DFT, ROOT
from redox.core.io import atomic_write_json
from redox.core.protocol import software_versions

EVIDENCE = ROOT / "manuscript" / "evidence" / f"{time.strftime('%Y-%m-%d')}_gpu_checks"

SCRATCH = Path(os.environ.get("TMPDIR", "/tmp")) / "bmepz_anchor"
SCRATCH.mkdir(parents=True, exist_ok=True)
BMEPZ = "COCCN1c2ccccc2N(CCOC)c2ccccc21"   # 5,10-bis(2-methoxyethyl)-5,10-dihydrophenazine
STATES = [("neu", 0, 1), ("ox1", 1, 2), ("ox2", 2, 1)]
MEAS = {"E1": -0.29, "E2": 0.50}  # V vs Fc, MeCN

def start_xyz():
    m = Chem.AddHs(Chem.MolFromSmiles(BMEPZ))
    AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE); AllChem.MMFFOptimizeMolecule(m)
    p = SCRATCH / "bmepz_start.xyz"
    conf = m.GetConformer()
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    with open(p, "w") as f:
        f.write(f"{m.GetNumAtoms()}\n\n")
        for s, i in zip(syms, range(m.GetNumAtoms())):
            pos = conf.GetAtomPosition(i); f.write(f"{s} {pos.x} {pos.y} {pos.z}\n")
    return p

def _sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def fc_ref(prov=None):
    fg = {}
    for s in ("neu", "ox"):
        r = read_result("ferrocene", s)
        at = aseread(str(DFT / "ferrocene" / s / "opt.xyz"))
        th = _thermal_correction(at, int(r["charge"]), int(r["mult"]))["g_thermal_eV"]
        fg[s] = r["e_smd_eV"] + th
        if prov is not None:
            prov[s] = dict(e_smd_eV=r["e_smd_eV"], g_thermal_eV=th, G_eV=fg[s],
                           geometry_sha256=_sha256(DFT / "ferrocene" / s / "opt.xyz"),
                           result_sha256=_sha256(DFT / "ferrocene" / s / "result.json"))
    return fg["ox"] - fg["neu"]

def main():
    t_start = time.perf_counter()
    start = start_xyz()
    fc_prov = {}
    fcref = fc_ref(fc_prov)
    print(f"ferrocene legacy fcref = {fcref:.4f} eV", flush=True)
    G = {}
    states = {}
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    for name, q, mult in STATES:
        opt_out = SCRATCH / f"{name}_opt.xyz"
        t0 = time.perf_counter()
        res = dft_smd(start, q, mult, backend="gpu", do_opt=True, do_freq=False, opt_out=opt_out)
        e_smd = res["e_smd_eV"]
        at = aseread(str(opt_out))
        th = _thermal_correction(at, q, mult)["g_thermal_eV"]
        G[name] = e_smd + th
        kept = EVIDENCE / f"bmepz_{name}_opt.xyz"
        shutil.copy(opt_out, kept)
        states[name] = dict(
            charge=q, mult=mult, e_smd_eV=e_smd, e_gas_eV=res.get("e_gas_eV"),
            g_thermal_eV=th, G_eV=G[name], thermal_method="gfn2-xtb RRHO (legacy tier)",
            converged=dict(geometry=res.get("geometry_optimization_converged"),
                           smd=res.get("converged_smd"), gas=res.get("converged_gas")),
            scf_guesses=res.get("scf_guesses"),
            smd_scf_guess_spread_eV=res.get("smd_scf_guess_spread_eV"),
            smd_scf_escalated=res.get("smd_scf_escalated"),
            smd_s_squared=res.get("smd_s_squared"),
            opt_level=f"{res.get('opt_xc')}-{res.get('opt_disp')}/{res.get('opt_basis')} SMD",
            energy_protocol=res.get("energy_protocol"),
            n_atoms=res.get("n_atoms"), n_ao_sp=res.get("n_ao_sp"),
            timing_s=dict(res.get("timing_s") or {}, state_wall=round(time.perf_counter() - t0, 1)),
            geometry=str(kept.relative_to(ROOT)), geometry_sha256=_sha256(kept))
        print(f"[{name}] q={q} mult={mult}  e_smd={e_smd:.3f}  g_thermal={th:.3f}  G={G[name]:.3f} eV", flush=True)
    E1 = (G["ox1"] - G["neu"]) - fcref
    E2 = (G["ox2"] - G["ox1"]) - fcref
    print("\n=== BMEPZ parent, computed vs measured (V vs Fc) ===")
    print(f"  E1: calc {E1:+.3f}  meas {MEAS['E1']:+.3f}  err {E1-MEAS['E1']:+.3f}")
    print(f"  E2: calc {E2:+.3f}  meas {MEAS['E2']:+.3f}  err {E2-MEAS['E2']:+.3f}")
    mae = (abs(E1-MEAS['E1']) + abs(E2-MEAS['E2'])) / 2
    bias = ((E1-MEAS['E1']) + (E2-MEAS['E2'])) / 2
    print(f"  per-couple catholyte MAE = {mae:.3f} V   bias = {bias:+.3f} V")
    record = dict(
        purpose="family-specific calibration check for the DMPZ catholytes (manuscript Validation)",
        molecule="5,10-bis(2-methoxyethyl)-5,10-dihydrophenazine (BMEPZ)", smiles=BMEPZ,
        route="RDKit ETKDG(seed 0xC0FFEE)+MMFF start for every state -> DFT+SMD opt -> "
              "active-protocol SP -> GFN2-xTB RRHO (legacy ranking tier); single seed, "
              "no UMA/xTB-ALPB pre-optimization",
        measured=dict(MEAS, units="V vs Fc, MeCN", doi="10.1016/j.chempr.2019.07.006"),
        fc_reference_eV=fcref, fc_reference_records=fc_prov, states=states,
        computed=dict(E1_V=E1, E2_V=E2), error_V=dict(E1=E1 - MEAS["E1"], E2=E2 - MEAS["E2"]),
        MAE_V=mae, bias_V=bias,
        host=socket.gethostname(), wall_s=round(time.perf_counter() - t_start, 1),
        software=software_versions(), date=time.strftime("%Y-%m-%d"))
    out = EVIDENCE / "bmepz_anchor.json"
    atomic_write_json(out, record)
    print(f"wrote {out}", flush=True)

if __name__ == "__main__":
    main()
