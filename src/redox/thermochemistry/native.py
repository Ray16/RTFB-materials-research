"""Addressed, derivative-checked solution DFT thermochemistry.

Uses the geometry-level r2SCAN-D4/SVP(D)/SMD surface, full grid response and
the same nuclear qRRHO definition for every state and Fc. Legacy optimizer
JSON/geometry files are preserved. This module is opt-in until a matched cohort
and Fc have been computed; readers never mix this definition with legacy GFN2.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import os
import socket
import subprocess
import time
from importlib import metadata
from pathlib import Path

import numpy as np

from redox.core.common import DFT, ROOT, read_manifest, read_result
from redox.core.protocol import load_record, software_versions
from redox.qm.xtb_runtime import configured_binary
from redox.thermochemistry.vibrations import parse_xtb_thermal_eV, write_xtb_hessian
from redox.thermochemistry.xtb_diagnostic import spin_degeneracy_eV

PRODUCER_SOURCE_BYTES = Path(__file__).read_bytes()
PRODUCER_SOURCE_SHA256 = hashlib.sha256(PRODUCER_SOURCE_BYTES).hexdigest()

NATIVE_THERMAL = dict(
    schema=1,
    model="solution geometry-level DFT Hessian",
    xc="r2scan",
    disp="d4",
    basis="def2-svp",
    anion_basis="def2-svpd",
    solvent="acetonitrile",
    solvent_model="SMD",
    density_fit=True,
    grid_level=3,
    cphf_grid_level=3,
    grid_response=True,
    conv_tol=1e-12,
    conv_tol_grad=1e-8,
    temperature_K=298.15,
    rotor_cutoff_cm=50,
    spin_partition="unsplit ground spin multiplet",
    gradient_max_acceptance_Ha_Bohr=0.00045,
    min_frequency_acceptance_cm=-5,
    finite_difference_step_Bohr=0.0005,
    finite_difference_relative_acceptance=0.01,
    max_spin_contamination=0.1,
    df_batch_cap_GiB=6,
    validated_pyscf_version="2.14.0",
    validated_gpu4pyscf_version="1.8.1",
)
XTB = configured_binary()
THERMAL_DEFINITION = "native_solution_dft_qrrho_spin"


def installed_versions():
    """Inspect distribution metadata without importing PySCF/CUDA libraries."""
    result = {}
    packages = (
        ("pyscf", ("pyscf",)),
        ("gpu4pyscf", ("gpu4pyscf", "gpu4pyscf-cuda12x", "gpu4pyscf-cuda11x")),
    )
    for name, candidates in packages:
        found = set()
        for distribution in candidates:
            try:
                found.add(metadata.version(distribution))
            except metadata.PackageNotFoundError:
                pass
        result[name] = next(iter(found)) if len(found) == 1 else None
    return result


def select_policy(sd, definition=THERMAL_DEFINITION):
    if definition not in (THERMAL_DEFINITION, "legacy_gfn2_gas_rrho"):
        raise ValueError("Unknown thermal definition")
    p = Path(sd) / "thermal_policy.json"
    policy = dict(
        schema=1,
        thermal_definition=definition,
        selected_at=time.time(),
        native_protocol=NATIVE_THERMAL if definition == THERMAL_DEFINITION else None,
    )
    with (Path(sd) / ".thermal_policy.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            existing = json.loads(p.read_text())
            if (
                existing.get("schema") == policy["schema"]
                and existing.get("thermal_definition") == policy["thermal_definition"]
                and existing.get("native_protocol") == policy["native_protocol"]
            ):
                return
        except (OSError, ValueError, TypeError, AttributeError):
            pass
        tmp = p.with_name(f".{p.name}.{os.getpid()}.{time.time_ns()}.tmp")
        try:
            tmp.write_text(json.dumps(policy, indent=2))
            tmp.replace(p)
        finally:
            tmp.unlink(missing_ok=True)


def identity(sd, q, m, protocol=NATIVE_THERMAL):
    return dict(
        geometry_sha256=hashlib.sha256((Path(sd) / "opt.xyz").read_bytes()).hexdigest(),
        charge=int(q),
        mult=int(m),
        protocol=dict(protocol),
    )


def cache_path(sd, q, m, protocol=NATIVE_THERMAL):
    ident = identity(sd, q, m, protocol)
    key = hashlib.sha256(json.dumps(ident, sort_keys=True).encode()).hexdigest()[:16]
    return Path(sd) / "native_thermal" / f"{key}.json", ident


def _is_finite_number(value) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(value)


def _probe_atom_index(symbols: list[str]) -> int:
    return next(
        (i for i, sym in enumerate(symbols) if sym == "O"),
        next((i for i, sym in enumerate(symbols) if sym == "Fe"), len(symbols) - 1),
    )


def _derivative_checks_ok(record: dict, protocol: dict, symbols: list[str]) -> bool:
    checks = record.get("derivative_checks", [])
    tol = protocol["finite_difference_relative_acceptance"]
    if any(not _is_finite_number(c.get("relative_error")) or not 0 <= c["relative_error"] <= tol for c in checks):
        return False
    probe = _probe_atom_index(symbols)
    required = {(ia, axis, step) for ia, axis in ((0, 1), (probe, 2)) for step in (0.0005, 0.001)}
    observed = {(c.get("atom_index"), c.get("axis"), c.get("step_Bohr")) for c in checks}
    return required <= observed


def _load_native_surface_record(sd, q, m, protocol=NATIVE_THERMAL, *, controls_matcher=None):
    """Admit the backing derivative record; its legacy partition is not primary."""
    try:
        p, ident = cache_path(sd, q, m, protocol)
        r = json.loads(p.read_text())
        if r.get("identity") != ident or r.get("status") != "derivative_checked_solution_thermal":
            return None
        from redox.backends.hessian_memory import record_budget_ok, record_derivatives_ok, record_pcm_ok

        if not record_budget_ok(r, protocol) or not record_pcm_ok(r, protocol):
            return None
        if r.get("stationary_minimum") is not True or r.get("converged") is not True or r.get("grid_response") is not True:
            return None
        g = r.get("g_thermal_total_eV")
        freqs = r.get("frequencies_cm")
        if not _is_finite_number(g):
            return None
        if not isinstance(freqs, list) or not freqs or not all(_is_finite_number(x) for x in freqs):
            return None
        if min(freqs) < protocol["min_frequency_acceptance_cm"]:
            return None
        lines = (Path(sd) / "opt.xyz").read_text().splitlines()
        nat = int(lines[0])
        if not record_derivatives_ok(r, protocol, nat):
            return None
        symbols = [line.split()[0] for line in lines[2 : 2 + nat]]
        coords = np.array([[float(x) for x in line.split()[1:4]] for line in lines[2 : 2 + nat]])
        linear = np.linalg.matrix_rank(coords - coords.mean(axis=0), tol=1e-7) <= 1
        if len(freqs) != 3 * nat - (5 if linear else 6):
            return None
        asym = r.get("maximum_asymmetry", math.inf)
        if not _is_finite_number(asym) or not 0 <= asym <= 1e-6:
            return None
        nuclear = r.get("g_nuclear_qrrho_eV")
        electronic_spin = r.get("g_spin_partition_eV")
        if not (_is_finite_number(nuclear) and _is_finite_number(electronic_spin)):
            return None
        if abs(electronic_spin - spin_degeneracy_eV(m, protocol["temperature_K"])) > 1e-8:
            return None
        if abs(g - nuclear - electronic_spin) > 1e-8:
            return None
        if abs((r.get("nuclear_by_cutoff_eV") or {}).get("50", math.inf) - nuclear) > 1e-8:
            return None
        expected_sw = {
            "pyscf": protocol["validated_pyscf_version"],
            "gpu4pyscf": protocol["validated_gpu4pyscf_version"],
        }
        if r.get("software") != expected_sw:
            return None
        gradient = r.get("gradient_max", math.inf)
        if not _is_finite_number(gradient) or not 0 <= gradient <= protocol["gradient_max_acceptance_Ha_Bohr"]:
            return None
        if not _derivative_checks_ok(r, protocol, symbols):
            return None
        spin_diag = r.get("smd_spin_diagnostics", {})
        spin = spin_diag.get("spin_contam", math.inf)
        if not _is_finite_number(spin) or abs(spin) > protocol["max_spin_contamination"]:
            return None
        if any(not _is_finite_number(spin_diag.get(name)) for name in ("homo_eV", "s_squared")):
            return None
        expected_controls = implementation_controls()
        stored_controls = r.get("implementation_controls")
        if controls_matcher is None:
            if stored_controls != expected_controls:
                return None
        elif not controls_matcher(stored_controls, expected_controls):
            return None
        return r
    except (OSError, ValueError, TypeError, KeyError, RuntimeError):
        return None


def load_native_thermal(sd, q, m, protocol=NATIVE_THERMAL):
    """Admit native thermochemistry only with its corrected internal partition."""
    raw = _load_native_surface_record(sd, q, m, protocol)
    if raw is None:
        return None
    from redox.thermochemistry.projected import load_projected_thermal, sha

    path, _ = cache_path(sd, q, m, protocol)
    partition = load_projected_thermal(path, Path(sd) / "opt.xyz")
    if partition is None:
        return None
    result = dict(raw)
    result.update(
        raw_native_g_thermal_total_eV=raw["g_thermal_total_eV"],
        nuclear_by_cutoff_eV=partition["nuclear_by_cutoff_eV"],
        frequencies_cm=partition["frequencies_cm"],
        g_nuclear_qrrho_eV=partition["g_nuclear_qrrho_eV"],
        g_spin_partition_eV=partition["g_spin_partition_eV"],
        g_thermal_total_eV=partition["g_thermal_total_eV"],
        thermal_partition_protocol=partition["identity"]["protocol"],
        thermal_postprocessing_path=partition["postprocessing_path"],
        thermal_postprocessing_sha256=sha(partition["postprocessing_path"]),
        thermal_partition_address=path.stem + "-" + Path(partition["postprocessing_path"]).stem,
    )
    return result


def free_energy_native(gid, state):
    """Return a native-definition free energy only when both addressed records exist."""
    raw = read_result(gid, state, raw=True)
    if raw is None:
        return None, "missing_geometry"
    selected = read_result(gid, state)
    if selected.get("geometry_spin_qc", "ok") != "ok":
        return None, selected["geometry_spin_qc"]
    sd = DFT / gid / state
    energy = load_record(sd, raw["charge"], raw["mult"])
    thermal = load_native_thermal(sd, raw["charge"], raw["mult"])
    if energy is None:
        return None, "missing_active_energy"
    if thermal is None:
        return None, "missing_validated_native_thermal"
    return energy["e_smd_eV"] + thermal["g_thermal_total_eV"], "ok"


def implementation_controls():
    certificate = json.loads((ROOT / "config/native_thermal_validation.json").read_text())
    digest = certificate.pop("certificate_content_sha256", None)
    payload = json.dumps(certificate, sort_keys=True, separators=(",", ":")).encode()
    if hashlib.sha256(payload).hexdigest() != digest:
        raise RuntimeError("Native implementation certificate integrity failed")
    if certificate.get("schema") != 1 or certificate.get("native_protocol") != NATIVE_THERMAL:
        raise RuntimeError("Native implementation certificate protocol mismatch")
    versions = installed_versions()
    if versions != certificate.get("software"):
        raise RuntimeError("Native Hessian library version changed; repeat the implementation controls")
    expected = [
        "analytic_full_response_verified",
        "uks_full_response_verified",
        "fine_response_implementation_verified",
    ]
    controls = certificate.get("controls", [])
    if [c.get("admission_status") for c in controls] != expected:
        raise RuntimeError("Incomplete implementation certificate")
    proof = []
    for control in controls:
        p = ROOT / control["relative_path"]
        error = control["relative_full_matrix_error"]
        if control.get("software") != versions or not _is_finite_number(error) or not 0 <= error < 0.01:
            raise RuntimeError("Invalid certified derivative control")
        if p.exists():
            if hashlib.sha256(p.read_bytes()).hexdigest() != control["source_sha256"]:
                raise RuntimeError("Control evidence differs from immutable certificate")
            r = json.loads(p.read_text())
            if r.get("status") != control["source_status"]:
                raise RuntimeError("Control status mismatch")
        item = dict(path=str(p), sha256=control["source_sha256"], status=control["admission_status"])
        if control["admission_status"] == "fine_response_implementation_verified":
            if control.get("fine_response_implementation_verified") is not True:
                raise RuntimeError("Fine response control not certified")
            item["relative_full_matrix_error"] = error
        proof.append(item)
    return proof


def ensure_native_thermal(gid, state, backend="gpu", state_dir=None, memory_cap_gib=None, qualification=False):
    if backend != "gpu":
        raise ValueError("This analytical implementation was validated on GPU only")
    from ase.io import read
    from pyscf import gto
    from pyscf.hessian import thermo
    from redox.backends.dft import _build_smd_mf, _ecp_for, _free_gpu_pool, _opt_protocol_mismatch, _orbital_diag
    from redox.backends.hessian_memory import execution_budget, stage_budget_bytes

    proof = implementation_controls()
    sd = Path(state_dir) if state_dir is not None else DFT / gid / state
    if qualification and (state_dir is None or DFT.resolve() in sd.resolve().parents):
        raise ValueError("Memory qualification must use an isolated directory outside production")
    memory = execution_budget(NATIVE_THERMAL, memory_cap_gib, qualification)
    raw = json.loads((sd / "result.json").read_text()) if (sd / "result.json").exists() else None
    if raw is None or _opt_protocol_mismatch(raw, raw["charge"]):
        raise RuntimeError("Source geometry missing or optimized at another level")
    q, m = raw["charge"], raw["mult"]
    p, ident = cache_path(sd, q, m)
    existing = _load_native_surface_record(sd, q, m)
    if existing:
        from redox.thermochemistry.projected import ensure_projected_thermal

        ensure_projected_thermal(p, sd / "opt.xyz")
        corrected = load_native_thermal(sd, q, m)
        if corrected is None:
            raise RuntimeError("Corrected native partition failed admission")
        return "present", corrected
    if p.exists():
        raise RuntimeError(f"Invalid addressed thermal record must be preserved/inspected: {p}")
    manifest = next(r for r in read_manifest() if r["id"] == gid and r["state"] == state)
    from redox.core.common import smiles_to_xyz_map

    if smiles_to_xyz_map(manifest["smiles"], sd / "opt.xyz") is None:
        raise RuntimeError("Source has changed declared connectivity")
    atoms = read(sd / "opt.xyz")
    symbols = atoms.get_chemical_symbols()
    coordinates = atoms.positions / 0.52917721067
    basis = NATIVE_THERMAL["anion_basis"] if q < 0 else NATIVE_THERMAL["basis"]
    work = p.with_suffix("")
    if work.exists():
        raise RuntimeError(f"Partial native thermal outputs must be inspected/preserved before rerun: {work}")
    work.mkdir(parents=True)
    (work / "producer_source.py").write_bytes(PRODUCER_SOURCE_BYTES)
    baseline = None

    def cpu(x):
        return np.asarray(x.get() if hasattr(x, "get") else x)

    import gpu4pyscf.df.hessian.rhf as dr
    import gpu4pyscf.df.hessian.uhf as du

    for module in (dr, du):
        if hasattr(module, "get_avail_mem"):
            previous = getattr(module, "_redox_native_original_avail_mem", module.get_avail_mem)
            module._redox_native_original_avail_mem = previous

            def budgeted_available(*args, _old=previous, _module=module.__name__, **kwargs):
                import sys

                return min(_old(*args, **kwargs), stage_budget_bytes(memory, _module, sys._getframe(1).f_code.co_name))

            module.get_avail_mem = budgeted_available

    def solve(coords, warm=True):
        mol = gto.M(
            atom=list(zip(symbols, coords.tolist())),
            unit="Bohr",
            charge=q,
            spin=m - 1,
            basis=basis,
            ecp=_ecp_for(symbols, basis),
            verbose=0,
            max_memory=6000,
        )
        mf = _build_smd_mf(mol, "r2scan", "acetonitrile", backend, disp="d4")
        mf.grids.level = 3
        mf.conv_tol = 1e-12
        mf.conv_tol_grad = 1e-8
        mf.max_cycle = 200
        mf.kernel(dm0=baseline if warm and baseline is not None else None)
        if not mf.converged:
            raise RuntimeError("Native thermal SCF failed")
        return mol, mf

    started = time.time()
    mol, mf = solve(coordinates)
    diag = _orbital_diag(mf, m - 1)
    if abs(diag.get("spin_contam", math.inf)) > 0.1:
        raise RuntimeError("SMD spin contamination failed")
    grad = mf.nuc_grad_method()
    grad.grid_response = True
    g = cpu(grad.kernel())
    if abs(g).max() > 0.00045:
        raise RuntimeError("Source is not stationary under the conservative solution gradient; inspect before replacing geometry")
    baseline = cpu(mf.make_rdm1())
    mf.cphf_grids = mf.grids.copy()
    hobj = mf.Hessian()
    hobj.grid_response = True
    hobj.max_cycle = 100
    H = cpu(hobj.kernel())
    mat = H.transpose(0, 2, 1, 3).reshape(3 * len(atoms), 3 * len(atoms))
    asym = float(abs(mat - mat.T).max())
    np.save(work / "hessian_cartesian.npy", mat)
    (work / "hessian_diagnostics.json").write_text(
        json.dumps(dict(maximum_asymmetry=asym, all_finite=bool(np.isfinite(mat).all()), execution_memory=memory), indent=2)
    )
    if not np.isfinite(mat).all() or asym > 1e-6:
        raise RuntimeError("Native Hessian symmetry/finite-value failure")
    # The gradient also owns its base SCF object. Release it before displaced
    # solves, otherwise the original GPU density-fitting tensors remain live.
    del mf, hobj, grad
    _free_gpu_pool()

    probe = _probe_atom_index(symbols)
    checks = []
    for ia, axis in ((0, 1), (probe, 2)):
        for step in (0.0005, 0.001):
            gradients = []
            for sign in (-1, 1):
                coords = coordinates.copy()
                coords[ia, axis] += sign * step
                _, fm = solve(coords, warm=False)
                gm = fm.nuc_grad_method()
                gm.grid_response = True
                gradients.append(cpu(gm.kernel()))
                del fm, gm
                _free_gpu_pool()
            numerical = (gradients[1] - gradients[0]) / (2 * step)
            analytic = H[ia, :, axis, :]
            error = float(np.linalg.norm(numerical - analytic) / np.linalg.norm(numerical))
            checks.append(dict(atom_index=ia, axis=axis, step_Bohr=step, relative_error=error))
            if error > 0.01:
                raise RuntimeError("Independent finite-gradient Hessian check failed")

    data = thermo.harmonic_analysis(mol, H, exclude_trans=True, exclude_rot=True, imaginary_freq=False)
    frequencies = np.asarray(data["freq_wavenumber"])
    if min(frequencies) < -5:
        raise RuntimeError("Source DFT geometry is not a minimum of its own solution surface")
    write_xtb_hessian(work / "hessian", (mat + mat.T) / 2)
    (work / "mol.xyz").write_bytes((sd / "opt.xyz").read_bytes())
    thermal = {}
    for cutoff in (25, 50, 100):
        log = work / f"thermo_{cutoff}.log"
        with log.open("w") as f:
            proc = subprocess.run(
                [str(XTB), "thermo", "--sthr", str(cutoff), "--temp", "298.15", "mol.xyz", "hessian"],
                cwd=work,
                stdout=f,
                stderr=subprocess.STDOUT,
                timeout=120,
            )
        if proc.returncode:
            raise RuntimeError("Native Hessian qRRHO failed")
        thermal[str(cutoff)] = parse_xtb_thermal_eV(log.read_text())
    if identity(sd, q, m) != ident:
        raise RuntimeError("Source changed during thermochemistry")
    spin = spin_degeneracy_eV(m)
    rec = dict(
        id=gid,
        state=state,
        status="memory_budget_qualification" if qualification else "derivative_checked_solution_thermal",
        identity=ident,
        converged=True,
        stationary_minimum=True,
        grid_response=True,
        smd_spin_diagnostics=diag,
        gradient_max=float(abs(g).max()),
        frequencies_cm=frequencies.tolist(),
        derivative_checks=checks,
        maximum_asymmetry=asym,
        g_nuclear_qrrho_eV=thermal["50"],
        g_spin_partition_eV=spin,
        g_thermal_total_eV=thermal["50"] + spin,
        nuclear_by_cutoff_eV=thermal,
        implementation_controls=proof,
        software=software_versions(),
        host=socket.gethostname(),
        wall_s=time.time() - started,
        source_producer_sha256=PRODUCER_SOURCE_SHA256,
        source_producer_snapshot=str(work / "producer_source.py"),
        producer_hash_convention="source snapshot at module import",
        execution_memory=memory,
        experimental_accuracy_validated=False,
        legacy_production_modified=False,
    )
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(rec, indent=2))
    tmp.replace(p)
    if not qualification:
        from redox.thermochemistry.projected import ensure_projected_thermal

        ensure_projected_thermal(p, sd / "opt.xyz")
        corrected = load_native_thermal(sd, q, m)
        if corrected is None:
            raise RuntimeError("Corrected native partition failed admission")
        return "computed", corrected
    return "computed", rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", required=True, help="comma separated IDs or id:state")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument(
        "--select-policy",
        action="store_true",
        help="durably select native definition; missing records never fall back",
    )
    ap.add_argument(
        "--repair-minimum",
        action="store_true",
        help="bounded signed-mode repair after actual DFT stationarity failure",
    )
    a = ap.parse_args()
    if len(os.sched_getaffinity(0)) != a.threads or os.getloadavg()[0] + a.threads > (os.cpu_count() or 1):
        ap.error("pinned CPU allocation/capacity required")
    if not os.getenv("CUDA_VISIBLE_DEVICES"):
        ap.error("run through gpu_reserve, which selects the GPU")
    selected = set(a.only.split(","))
    rows = [r for r in read_manifest() if r["id"] in selected or r["id"] + ":" + r["state"] in selected]
    available = {name for r in read_manifest() for name in (r["id"], r["id"] + ":" + r["state"])}
    if selected - available:
        ap.error(f"unknown declared inputs: {sorted(selected - available)}")
    if not rows:
        ap.error("no declared states selected")
    if a.select_policy:
        for row in rows:
            select_policy(DFT / row["id"] / row["state"])
    failures = 0
    total_failures = 0
    for row in rows:
        with (ROOT / "logs/pipeline_locks" / f"{row['id']}_native_thermal.lock").open("a+") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            try:
                try:
                    kind, record = ensure_native_thermal(row["id"], row["state"])
                except RuntimeError as exc:
                    repairable = "not stationary" in str(exc) or "not a minimum" in str(exc)
                    sd = DFT / row["id"] / row["state"]
                    if "Partial native thermal outputs" in str(exc):
                        gh = hashlib.sha256((sd / "opt.xyz").read_bytes()).hexdigest()
                        history = (
                            ROOT
                            / "results/validation/physical_thermal_spin_2026-10-10/native_minimum_repairs"
                            / row["id"]
                            / row["state"]
                        )
                        for p in history.glob("*/status.json"):
                            evidence = json.loads(p.read_text())
                            reason = evidence.get("reason", "")
                            if evidence.get("source_sha256") == gh and (
                                "not stationary" in reason or "not a minimum" in reason
                            ):
                                repairable = True
                    if not a.repair_minimum or not repairable:
                        raise
                    from redox.thermochemistry.native_repair import repair_minimum

                    kind, record = repair_minimum(row["id"], row["state"], str(exc))
                print(kind, row["id"], row["state"], record["g_thermal_total_eV"], flush=True)
                failures = 0
            except Exception as exc:
                print("FAILED", row["id"], row["state"], repr(exc), flush=True)
                failures += 1
                total_failures += 1
                if failures >= 3:
                    raise SystemExit("Three consecutive actual failures; native thermal circuit breaker")
    if total_failures:
        raise SystemExit(f"{total_failures} states failed native thermal validation")


if __name__ == "__main__":
    main()
