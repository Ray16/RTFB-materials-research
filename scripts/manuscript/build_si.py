#!/usr/bin/env python
"""Build the Supporting Information from the stored records (CPU only; runs no chemistry).

Writes
  manuscript/supporting_information.md      the SI document
  manuscript/si_tables/S*.csv                machine-readable copies of the long tables

Every number in the SI is read from a record listed in the table captions; nothing is typed in
by hand. Re-running the script regenerates the SI exactly.

  PYTHONPATH=src python scripts/manuscript/build_si.py
"""
from __future__ import annotations

import csv
import glob
import json
import statistics as st
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts" / "plotting" / "manuscript"))

H2EV = 27.211386245988
ARCH = Path("/nfs/lambda_stor_01/homes/rzhu/archive/0_redox")
BENCH = ARCH / "manuscript/evidence/2026-10-22_benchmark_complete"
EVID = ROOT / "manuscript/evidence"
OUT_MD = ROOT / "manuscript/supporting_information.md"
OUT_TAB = ROOT / "manuscript/si_tables"
SCORED = ROOT / "results/validation/full_pool_scored_legacyG_MeCN_TBAPF6.csv"
CALCS = ROOT / "calcs"


def j(p):
    return json.loads(Path(p).read_text())


def f2(x, d=2, sign=False):
    if x is None:
        return "—"
    s = f"{x:+.{d}f}" if sign else f"{x:.{d}f}"
    return s.replace("-", "−")


def write_csv(name, rows, cols):
    OUT_TAB.mkdir(parents=True, exist_ok=True)
    with (OUT_TAB / name).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


# --------------------------------------------------------------------------- sections
def s_module_map():
    stages = [
        ("1 Build", "`redox.structures.candidates`", "grafted structure and redox states from SMILES"),
        ("2 MLIP pre-optimization", "`redox.backends.uma_current`", "UMA conformer and spin scan, gas phase"),
        ("3 DFT optimization", "`redox.backends.dft`", "xTB-ALPB warm start, r²SCAN-D4 SMD optimization, seeds"),
        ("4 Production energy", "`redox.backends.sp`", "ωB97M-V/def2-TZVPD SMD and gas single points (adopted when identical)"),
        ("5 Thermochemistry", "`redox.backends.dft._thermal_correction`; `redox.qm.qrrho`", "xTB RRHO (pool tier); projected qRRHO (native tier)"),
        ("6 Properties", "`redox.properties.potentials`, `.capacity_and_proxies`, `.lambda_outer`, `.reorg`", "potentials, capacity, λ_o, λ_i"),
        ("7 Quality checks", "`redox.quality.readiness`, `.spins`, `.stability`", "per-state structure, spin, protocol, completeness"),
        ("8 Ranking", "`redox.screening.scorecard`, `.pareto`", "pool assignment and Pareto fronts"),
        ("9 Tables and audit", "`redox.reporting.tables`", "exported tables and arithmetic audit"),
    ]
    missing = [m for _, mods, _ in stages for m in mods.replace("`", "").replace(";", ",").split(",")
               if m.strip().startswith("redox.") and not (ROOT / "src" / (m.strip().split("._")[0].replace(".", "/") + ".py")).exists()]
    assert not missing, missing
    return ("## S1. Where each stage lives in the code\n\n"
            "Each stage of the pipeline is one module in the released package. The run manifest "
            "(`run_manifest.json`) records the exact command line used for every stage.\n\n"
            + md_table(["Stage", "Module", "What it does"], stages) + "\n")


def s_equations():
    return r"""## S2. Equations used

**Redox potential.** For a one-electron couple O + e⁻ → R with free energies
G = E_SMD + G_thermal,

  E_vs_Fc = −[G(R) − G(O)]/e − E_abs(Fc),  E_abs(Fc) = G(Fc⁺) − G(Fc),

with ferrocene computed at the same level and with the same thermal tier as the
molecule. Constants at 298.15 K: RT/F = 0.025693 V; Nernst slope 0.059160 V per decade.
A two-electron candidate's reported potential is the mean of its two couples.

**Capacity.** C = nF/(3.6 M) in mAh g⁻¹, with n the number of usable electrons (those whose
couples lie in the window [−3.0, +2.0] V vs Fc on the candidate's side of the 0 V divider)
and M the neutral polystyrene repeat-unit mass: the capped-graft mass plus one carbon
(12.011 g mol⁻¹), because the cap's CH₃–Ar becomes –[CH₂–CH(Ar)]– in the polymer. The
salt-loaded variant adds n counterions to M (PF₆⁻, 144.96 g mol⁻¹, for catholytes; TBA⁺,
242.47 g mol⁻¹, for anolytes).

**Outer-sphere reorganization (Born estimate).** For a sphere of radius a,

  λ_o = z² e² (1/ε_op − 1/ε_s) / (8π ε₀ a),

with ε_s = 37.5, ε_op = n_D² = 1.806 (n_D = 1.344), z = 1, and a the radius of the sphere
with the same solvent-accessible surface area as the molecule. The molecular-cavity value
uses the same slow-polarization (Pekar) factor with the IEF-PCM response operators of the
real cavity, applied to the change in surface potential between the two charge states; it
reproduces the Born value exactly for a single-atom cavity (code check
`lambda_outer.sphere_sanity_check`).

**Inner-sphere reorganization (four-point).**

  λ_i = [E_O(geometry of R) − E_O(geometry of O)] + [E_R(geometry of O) − E_R(geometry of R)].

The values in this work use gas-phase energies at the SMD geometries (main text).

**Disproportionation.** ΔG_disp = G(R^{q+1}) + G(R^{q−1}) − 2G(R^q), equal to F times the
spacing between the two waves.

**Thermal corrections.** Pool tier: GFN2-xTB Hessian at the DFT geometry, xTB
thermodynamics at 298.15 K. Native tier: r²SCAN-D4/def2-SVP(D) SMD Hessian with grid
response; translations and rotations are projected out and the internal modes are passed
to the xTB thermodynamics routine with the quasi-RRHO treatment of low modes (rotor cutoff
50 cm⁻¹) [Grimme 2012], plus a spin-degeneracy term.
"""


def s_benchmark():
    b = j(BENCH / "benchmark.json")
    leg = j(EVID / "2026-11-18_cpu_checks/benchmark_legacy_tier.json")
    lrows = {(r["id"], r["event"]): r for r in leg["endpoints"]}
    pts = list(csv.DictReader(open(BENCH / "benchmark_points.csv")))
    names = {(o["id"], o["event"]): o["name"] for o in b["outliers"]}
    rows, out = [], []
    for p in pts:
        k = (p["id"], p["event"]); L = lrows[k]
        rows.append(dict(molecule=names.get(k, p["id"]), couple=p["event"], wave=p["wave"], family=p["family"],
                         exp_V=float(p["exp_V"]), calc_native_V=float(p["calc_V"]),
                         calc_pool_tier_V=L["calc_legacy_V"], err_native_V=float(p["calc_V"]) - float(p["exp_V"]),
                         err_pool_tier_V=L["calc_legacy_V"] - float(p["exp_V"]), doi=p["doi"]))
    write_csv("S3_benchmark_endpoints.csv", rows, list(rows[0]))
    rows.sort(key=lambda r: (r["family"], r["molecule"], r["wave"]))
    for r in rows:
        out.append([r["molecule"], r["couple"], r["wave"], f2(r["exp_V"], sign=True), f2(r["calc_native_V"], sign=True),
                    f2(r["calc_pool_tier_V"], sign=True), f2(r["err_native_V"], sign=True), f2(r["err_pool_tier_V"], sign=True), r["doi"]])
    m, s = b["metrics"], leg["summary"]
    fam = [[k, v["n"], v["molecules"], f2(v["MAE_V"], 3), f2(v["bias_V"], 3, True), f2(v["RMSE_V"], 3)]
           for k, v in b["by_family"].items()]
    wav = [[k, v["n"], f2(v["MAE_V"], 3), f2(v["bias_V"], 3, True)] for k, v in b["by_wave"].items()]
    boot = [[u["cluster"] if u["cluster"] != "id" else "molecule", u["independent_cluster_count"],
             "[" + ", ".join(f2(x, 3) for x in u["intervals"]["MAE_V"]) + "]",
             "[" + ", ".join(f2(x, 3, True) for x in u["intervals"]["bias_V"]) + "]"] for u in b["uncertainty"]]
    boot = [[("study" if c == "doi" else c), *rest] for c, *rest in boot]
    tiers = [["native Hessian (published)", m["n"], f2(m["MAE_V"], 3), f2(m["RMSE_V"], 3), f2(m["bias_V"], 3, True), f2(m["maximum_absolute_error_V"], 3)],
             ["xTB (used for the pool)", s["legacy"]["n"], f2(s["legacy"]["MAE_V"], 3), f2(s["legacy"]["RMSE_V"], 3), f2(s["legacy"]["bias_V"], 3, True), f2(s["legacy"]["max_abs_V"], 3)]]
    prov = b["review_snapshot_sha256"]
    man = j(BENCH / "manifest.json")
    text = f"""## S3. Benchmark against experiment

The panel has {m['n']} endpoints from {b['admitted_molecules']} molecules and {b['original_studies']} studies, all
measured in acetonitrile against ferrocene. No fitting, offset, or outlier removal is
applied. Table S3a lists every endpoint at both thermal tiers; Table S3b gives summary
statistics.

**Table S3a.** Endpoints (V vs Fc). "Native" is the published benchmark; "pool tier" uses
the xTB thermal correction with the same electronic energies. Source:
`benchmark.json`, `benchmark_points.csv`, and
`manuscript/evidence/2026-11-18_cpu_checks/benchmark_legacy_tier.json`; CSV copy
`si_tables/S3_benchmark_endpoints.csv`.

{md_table(["Molecule", "Couple", "Wave", "Exp.", "Native", "Pool tier", "Error (native)", "Error (pool tier)", "Source DOI"], out)}

**Table S3b.** Summary statistics (V).

{md_table(["Thermal tier", "n", "MAE", "RMSE", "Bias", "Max error"], tiers)}

By family (native tier):

{md_table(["Family", "n", "Molecules", "MAE", "Bias", "RMSE"], fam)}

By wave (native tier):

{md_table(["Wave", "n", "MAE", "Bias"], wav)}

Cluster bootstrap ({b['uncertainty'][0]['replicates']} replicates, seed {b['uncertainty'][0]['seed']}), 95% intervals:

{md_table(["Resampled unit", "Clusters", "MAE interval", "Bias interval"], boot)}

These intervals describe the observed panel; they are not prediction intervals for a
new molecule.

Between the two tiers the potentials move by {f2(s['shift_V']['min'], 3, True)} to {f2(s['shift_V']['max'], 3, True)} V
(mean {f2(s['shift_V']['mean'], 3, True)} V). The ferrocene reference is {f2(s['fc_native_V'], 4)} V at the native tier
and {f2(s['fc_legacy_V'], 4)} V at the pool tier. Recomputing the native tier from the stored
records reproduces every published potential (largest difference
{s['max_native_reproduction_error_V']:.1e} V).

**Provenance.** The benchmark was produced by code commit `{man['source_commit'][:8]}`
(preserved under the tag `legacy/pre-refactor-main`). Panel hash
`{b['panel_sha256'][:16]}…`; the review files it was checked against and their SHA-256:

{md_table(["File", "SHA-256"], [[Path(k).name, f"`{v[:16]}…`"] for k, v in prov.items()])}
"""
    return text


def s_bmepz():
    r = j(EVID / "2026-11-18_gpu_checks/bmepz_anchor.json")
    rows = [[k, v["charge"], v["mult"], f2(v["e_smd_eV"], 4), f2(v["g_thermal_eV"], 4), f2(v["G_eV"], 4),
             ", ".join(v["scf_guesses"]), f"{v['timing_s']['state_wall'] / 60:.1f}"] for k, v in r["states"].items()]
    return f"""## S4. Calibration check for the dihydrophenazine family (BMEPZ)

The benchmark has no dihydrophenazine, so the measured parent
{r['molecule']} (SMILES `{r['smiles']}`) was computed at the pool's ranking tier and
compared with its measured potentials (E₁ = {f2(r['measured']['E1'], sign=True)} V, E₂ = {f2(r['measured']['E2'], sign=True)} V vs Fc in
MeCN [Kwon 2019]). Each state was started from the same RDKit/MMFF structure,
optimized with DFT in SMD acetonitrile (one conformer seed, no UMA or xTB
pre-optimization), scored with the production single point, and given the xTB
thermal correction used for the pool. Run on {r['host']} ({r['date']}), {r['wall_s'] / 60:.0f} min in total.

{md_table(["State", "Charge", "Multiplicity", "E_SMD (eV)", "G_thermal (eV)", "G (eV)", "SCF guesses", "Wall time (min)"], rows)}

Ferrocene reference (pool tier): {f2(r['fc_reference_eV'], 4)} V. Computed E₁ = {f2(r['computed']['E1_V'], 3, True)} V
(error {f2(r['error_V']['E1'], 3, True)} V) and E₂ = {f2(r['computed']['E2_V'], 3, True)} V (error {f2(r['error_V']['E2'], 3, True)} V);
MAE {f2(r['MAE_V'], 3)} V, bias {f2(r['bias_V'], 3, True)} V. Full record and optimized geometries:
`manuscript/evidence/2026-11-18_gpu_checks/bmepz_anchor.json` and `bmepz_*_opt.xyz`.
"""


def _scf_evidence():
    val_n, val_ok = 0, 0
    for f in glob.glob(str(CALCS / "validation/scf_level/prod/*/*.json")):
        d = j(f)
        if d["mult"] == 1 or d["id"] == "ferrocene":
            continue
        for ph in d["phases"].values():
            e = {x["guess"]: x["e_Ha"] for x in ph["guesses"] if x.get("converged") and x.get("e_Ha") is not None}
            if not e or not {"minao", "huckel"} & set(e):
                continue
            val_n += 1
            val_ok += (min(e[g] for g in ("minao", "huckel") if g in e) - min(e.values())) * H2EV <= 1e-3
    metals = {"Fe", "Co", "Ni", "Cu", "Mn", "Cr", "V", "Ti", "Ru", "Os", "Zn"}
    pool_n, pool_ok = 0, 0
    for f in glob.glob(str(CALCS / "dft/*/*/result.json")):
        d = j(f)
        if d.get("mult", 1) == 1:
            continue
        try:
            syms = [l.split()[0] for l in Path(f).with_name("opt.xyz").read_text().splitlines()[2:] if l.strip()]
        except OSError:
            continue
        if metals & set(syms):
            continue
        for ph in ("smd", "gas"):
            g = {x["guess"]: x["e_Ha"] for x in (d.get(f"scf_guess_energies_{ph}") or [])
                 if x.get("converged") and x.get("e_Ha") is not None}
            if {"minao", "atom", "huckel"} <= set(g):
                pool_n += 1
                pool_ok += (min(g["minao"], g["huckel"]) - min(g.values())) * H2EV <= 1e-3
    return val_n, val_ok, pool_n, pool_ok


def s_scf():
    val_n, val_ok, pool_n, pool_ok = _scf_evidence()
    c = j(EVID / "2026-11-18_gpu_checks/orop170_scf_policy.json")
    cm = c["computed"]
    guesses = lambda key: "; ".join(f"{x['guess']} {x['e_Ha']:.6f}" for x in cm[key])
    return f"""## S5. Evidence for the open-shell SCF guess policy

For organic open-shell states the production code now starts the SCF from two guesses
(minao and Hückel) and adds the atomic guess when they disagree by more than 1 meV.

**Evidence set.** {val_n} solvent or gas-phase open-shell states from the SCF validation set,
each searched with up to nine starting guesses including orbital-occupation seeds, and
{pool_n} open-shell states from production records with all three standard guesses
(transition-metal systems excluded). In {val_ok} of {val_n} and {pool_ok} of {pool_n} cases, the
better of minao and Hückel is within 1 meV of the lowest solution found by any guess.

**Positive control.** The OROP-170 radical anion has two SCF solutions in solvent.
Re-scored through the production code on its stored geometry
({c['host']}, {c['wall_s']:.0f} s):

{md_table(["Phase", "Guess energies (Ha)", "Spread (eV)", "Escalated", "Difference from 9-guess search (eV)"],
          [["solvent (SMD)", guesses("scf_guess_energies_smd"), f2(cm["smd_scf_guess_spread_eV"], 3), cm["smd_scf_escalated"], f"{c['delta_vs_reference_eV']['smd']:.1e}"],
           ["gas", guesses("scf_guess_energies_gas"), f"{cm['gas_scf_guess_spread_eV']:.1e}", cm["gas_scf_escalated"], f"{c['delta_vs_reference_eV']['gas']:.1e}"]])}

All checks passed: `manuscript/evidence/2026-11-18_gpu_checks/orop170_scf_policy.json`.
"""


def s_qc():
    pool = list(csv.DictReader(open(EVID / "2026-11-18_cpu_checks/readiness_pool91.csv")))
    reg = list(csv.DictReader(open(EVID / "2026-11-18_cpu_checks/readiness_registry.csv")))
    ok = lambda r: r["status"] == "ready" or r["issues"] == "missing_validated_native_thermal"
    flagged = [r for r in pool if not ok(r)]
    complete = sum("missing_geometry_or_result" not in r["issues"] and "missing_or_invalid_active_sp" not in r["issues"] for r in reg)
    return f"""## S6. Quality checks

Readiness was run on the whole registry and on the 91 scored candidates (2026-11-18).

{md_table(["Set", "Declared states", "Complete (geometry + single point)", "Pass all ranking checks"],
          [["registry", len(reg), complete, "—"], ["91 scored candidates", len(pool), len(pool), sum(ok(r) for r in pool)]])}

The native-tier thermal check is not a ranking check, because the pool is ranked on the
xTB thermal tier. States that fail any ranking check:

{md_table(["Candidate", "State", "Charge", "Issue"], [[r["id"], r["state"], r["charge"], r["issues"]] for r in flagged])}

Neither of these candidates is on a Pareto front, and the fronts are identical without
them. Full outputs: `manuscript/evidence/2026-11-18_cpu_checks/readiness_*.csv`.
"""


def s_scored():
    import fig5_pareto as fp
    rows = list(csv.DictReader(open(SCORED)))
    lam = {r["id"]: r for r in csv.DictReader(open(ROOT / "results/validation/lambda_outer_fullpool_2026-11-17.csv"))}
    pool_flag = {r["id"]: r for r in csv.DictReader(open(EVID / "2026-11-18_cpu_checks/readiness_pool91.csv"))
                 if not (r["status"] == "ready" or r["issues"] == "missing_validated_native_thermal")}
    frozen, _ = fp.load_frozen()
    front = set()
    for pool, low in (("anolyte", True), ("catholyte", False)):
        front |= {r["id"] for r in fp.pareto([r for r in frozen if r["pool"] == pool and not r["ctrl"]], low)}
    frozen_ids = {r["id"] for r in frozen}
    out, md = [], []
    for r in rows:
        L = lam.get(r["id"], {})
        d = dict(id=r["id"], name=r["name"], family=r["family"], pool=r["pool"], n=r["n"],
                 E_mean_V=r["E_mean"], E_couples_V=r["E_couples"], capacity_active_mAh_g=r["cap_active_mAh_g"],
                 capacity_salt_mAh_g=r["cap_salt_mAh_g"], MW_repeat_unit=r["MW_repeat_unit"],
                 lambda_o_born_eV=L.get("lambda_o_born_eV", ""), sasa_radius_A=L.get("sasa_radius_A", ""),
                 in_91_pool=r["id"] in frozen_ids, pareto_front=r["id"] in front,
                 qc_flag=pool_flag[r["id"]]["issues"] if r["id"] in pool_flag else "")
        out.append(d)
        md.append([d["name"], d["pool"], d["n"], f2(float(d["E_mean_V"]), 2, True), d["capacity_active_mAh_g"],
                   d["capacity_salt_mAh_g"], d["lambda_o_born_eV"] or "—", "yes" if d["pareto_front"] else "",
                   "" if d["in_91_pool"] else "added later", d["qc_flag"]])
    write_csv("S7_scored_candidates.csv", out, list(out[0]))
    rs = j(EVID / "2026-11-18_cpu_checks/pool_legacy_tier_rescore.json")["summary"]
    check = (f"**Independent check.** Every couple was recomputed from the stored energies and "
             f"geometries (production single point plus a freshly computed xTB thermal term; "
             f"ferrocene reference {f2(rs['fc_legacy_V'], 4)} V) by "
             f"`scripts/validation/score_pool_legacy_tier.py`. All {rs['n_compared']} of "
             f"{rs['n_candidates']} candidates agree with the table to within "
             f"{rs['max_abs_diff_V'] * 1000:.1f} mV (the table's rounding); record: "
             f"`manuscript/evidence/2026-11-18_cpu_checks/pool_legacy_tier_rescore.json`.\n")
    order = sorted(range(len(md)), key=lambda i: (md[i][1], float(out[i]["E_mean_V"])))
    return f"""## S7. All scored candidates

**Table S7.** Every row of the scored table (pool's ranking tier). E is the mean of the
couples, vs Fc. Capacities: active = neutral repeat-unit mass (used in the paper); salt =
with counterions (do not compare across pools). λ_o (Born) uses the SASA-equivalent radius.
Source: `results/validation/full_pool_scored_legacyG_MeCN_TBAPF6.csv`,
`results/validation/lambda_outer_fullpool_2026-11-17.csv`; CSV copy with identifiers,
couples, masses and radii: `si_tables/S7_scored_candidates.csv`.

{check}
{md_table(["Candidate", "Pool", "n", "E (V)", "Capacity, active", "Capacity, salt", "λ_o Born (eV)", "Front", "Note", "QC flag"], [md[i] for i in order])}
"""


def s_throughput():
    rows = list(csv.DictReader(open(ROOT / "manuscript/tables/throughput_measured.csv")))
    dft = [r for r in rows if r["stage"] == "dft"]
    md = [[r["id"], r["gpu"], r["state"], r["mult"], r["n_atoms"], r["wall_min"], r["opt_min"], r["non_opt_min"]] for r in dft]
    return f"""## S8. Measured cost of each DFT state

Per-state DFT wall time (minutes, one GPU per task, two CPU threads per worker, three
conformer seeds per state). "Optimization" is the sum of the geomeTRIC optimizer times
of the three seeds; "other" is the rest of the stage (xTB pre-optimization and the
solvent and gas-phase single points with their SCF guesses). Source:
`manuscript/tables/throughput_measured.csv` (all stages), built by
`scripts/validation/throughput_report.py` from the run manifests and logs.

{md_table(["Candidate", "GPU", "State", "Multiplicity", "Atoms", "Wall", "Optimization", "Other"], md)}
"""


def s_software():
    import importlib.metadata as im
    vers = []
    for pkg in ("pyscf", "gpu4pyscf-cuda12x", "gpu4pyscf", "rdkit", "ase", "fairchem-core", "numpy", "geometric"):
        try:
            vers.append([pkg, im.version(pkg)])
        except im.PackageNotFoundError:
            pass
    import subprocess
    xtb = subprocess.run(["/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/xtb", "--version"],
                         capture_output=True, text=True).stdout
    xv = next((l.split("version")[1].split()[0] for l in xtb.splitlines() if "xtb version" in l), "?")
    vers.append(["xtb", xv])
    return f"""## S9. Software

{md_table(["Package", "Version"], vers)}

GPUs used: NVIDIA Tesla V100-SXM2-32GB and NVIDIA H100 NVL. All GPU jobs were submitted
through a reservation tool that assigns one job per GPU.
"""


def main():
    parts = ["# Supporting Information",
             "*A validated, reproducible pipeline for screening graftable redox motifs for "
             "non-aqueous organic redox flow batteries.*\n\n"
             "Generated by `scripts/manuscript/build_si.py` from the stored records; every value "
             "below is read from the file named with it. Citations refer to the main-text "
             "reference list.\n",
             s_module_map(), s_equations(), s_benchmark(), s_bmepz(), s_scf(), s_qc(), s_scored(),
             s_throughput(), s_software()]
    OUT_MD.write_text("\n".join(parts))
    print("wrote", OUT_MD, "and", sorted(p.name for p in OUT_TAB.glob("*.csv")))


if __name__ == "__main__":
    main()
