"""Recover an isolated native calculation from a source-bound saved Hessian.

Only an independently replayed additive CDS replacement is supported. Standard
SCF, stationarity, finite-gradient, minimum and thermodynamic gates still run.
"""

import json
from pathlib import Path

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def reconstructed_matrix(spec):
    required = {
        spec[name]
        for name in (
            "source_hessian",
            "coarse_CDS",
            "fine_CDS",
            "refined_hessian",
            "reconstruction_replay",
        )
    }
    if not required <= set(spec["input_sha256"]):
        raise ValueError("Reconstruction input hash coverage is incomplete")
    api.verify_sources(spec["input_sha256"])
    replay = json.loads(Path(spec["reconstruction_replay"]).read_text())
    if replay["status"] != "two_full_recomputation_replays_passed" or {
        r["mult"] for r in replay["controls"]
    } != {1, 2}:
        raise ValueError("Both full recomputation controls are required")
    for row in replay["controls"]:
        if not row.get("input_sha256"):
            raise ValueError("Control replay input evidence is empty")
        api.verify_sources(row["input_sha256"])
        if any(
            not 0 <= row[k] <= 1e-6
            for k in (
                "maximum_matrix_difference",
                "relative_matrix_error",
                "maximum_thermal_difference_eV",
            )
        ):
            raise ValueError("Reconstruction control exceeds its fixed tolerance")
    source = np.load(spec["source_hessian"])
    coarse = np.load(spec["coarse_CDS"])
    fine = np.load(spec["fine_CDS"])
    matrix = source - coarse + fine
    if (
        any(a.shape != source.shape for a in (coarse, fine))
        or not np.isfinite(matrix).all()
    ):
        raise ValueError("CDS reconstruction has invalid shape/values")
    if np.max(abs(matrix - matrix.T)) > 1e-6:
        raise ValueError("Reconstructed Hessian fails the unchanged symmetry gate")
    if not np.array_equal(matrix, np.load(spec["refined_hessian"])):
        raise ValueError("Reconstructed matrix differs from the diagnosed replacement")
    return matrix, list(
        map(int, np.unravel_index(abs(coarse - coarse.T).argmax(), coarse.shape))
    )


def verify_reconstruction_record(record, geometry):
    settings = record["execution_hessian_reconstruction"]
    if settings["module_sha256"] != api.digest(__file__) or settings[
        "source_geometry_sha256"
    ] != api.digest(geometry):
        raise ValueError("Saved-Hessian reconstruction source/geometry differs")
    matrix, indices = reconstructed_matrix(settings["spec"])
    if api.digest(settings["final_hessian"]) != settings[
        "final_hessian_sha256"
    ] or not np.array_equal(np.load(settings["final_hessian"]), matrix):
        raise ValueError("Final Hessian is not the checked additive reconstruction")
    derivatives = record["execution_derivatives"]
    required = {(column, step) for column in indices for step in (0.0005, 0.001)}
    observed = {
        (c["column"], c["step_Bohr"]) for c in derivatives["worst_contact_checks"]
    }
    if not required <= observed:
        raise ValueError("Worst-column validation coverage is incomplete")
    for check in derivatives["worst_contact_checks"]:
        column, step = check["column"], check["step_Bohr"]
        arrays = []
        for sign in (-1, 1):
            matches = [
                Path(p)
                for p in derivatives["target_evidence_sha256"]
                if p.endswith(f"gradient_{column}_{step}_{sign}.npy")
            ]
            if len(matches) != 1:
                raise ValueError("Worst-column derivative evidence missing")
            arrays.append(np.load(matches[0]))
        numerical = (arrays[1] - arrays[0]) / (2 * step)
        error = float(
            np.linalg.norm(numerical - matrix[column]) / np.linalg.norm(numerical)
        )
        if (
            not np.isfinite(error)
            or error > 0.01
            or abs(error - check["relative_error"]) > 1e-10
        ):
            raise ValueError("Saved worst-column derivative check does not reproduce")
    return matrix


def worst_column_checks(state, matrix, indices, protocol):
    from ase.io import read
    from pyscf import gto
    from redox.qm.dft import _build_smd_mf, _free_gpu_pool
    from redox.qm.stable_geometry_scf import stabilize

    atoms = read(state / "opt.xyz")
    coordinates = atoms.positions / 0.52917721067
    directory = state / "CDS_worst_column_checks"
    directory.mkdir(exist_ok=False)
    checks = []
    evidence = {}
    for column in sorted(set(indices)):
        atom, axis = divmod(column, 3)
        for step in (0.0005, 0.001):
            gradients = []
            for sign in (-1, 1):
                shifted = coordinates.copy()
                shifted[atom, axis] += sign * step
                mol = gto.M(
                    atom=list(zip(atoms.get_chemical_symbols(), shifted.tolist())),
                    unit="Bohr",
                    basis="def2-svp",
                    charge=0,
                    spin=0,
                    verbose=0,
                    max_memory=6000,
                )
                mf = _build_smd_mf(mol, "r2scan", "acetonitrile", "gpu", disp="d4")
                mf.grids.level = protocol["grid_level"]
                mf.conv_tol = protocol["conv_tol"]
                mf.conv_tol_grad = protocol["conv_tol_grad"]
                mf.max_cycle = 300
                stabilize(mf, directory / "SCF_recovery.jsonl")
                mf.kernel()
                if not mf.converged:
                    raise RuntimeError("Worst-column SCF failed")
                method = mf.nuc_grad_method()
                method.grid_response = True
                gradient = method.kernel()
                gradient = (
                    np.asarray(gradient.get() if hasattr(gradient, "get") else gradient)
                    .ravel()
                    .copy()
                )
                path = directory / f"gradient_{column}_{step}_{sign}.npy"
                np.save(path, gradient)
                evidence[str(path)] = api.digest(path)
                gradients.append(gradient)
                atomic_write_json(
                    directory / "progress.json",
                    dict(completed_evaluations=len(evidence), output_sha256=evidence),
                )
                del mf, method, mol
                _free_gpu_pool()
            numerical = (gradients[1] - gradients[0]) / (2 * step)
            error = float(
                np.linalg.norm(numerical - matrix[column]) / np.linalg.norm(numerical)
            )
            checks.append(dict(column=column, step_Bohr=step, relative_error=error))
            if not np.isfinite(error) or error > 0.01:
                raise RuntimeError("Worst-contact finite-gradient check failed")
    return checks, evidence


def annotate_finer_hessian(path, geometry):
    """Attach qualified step provenance and independent worst-column evidence."""
    from redox.qm import native_thermal as native
    from redox.qm.hessian_memory import derivative_certificate

    path = Path(path)
    state = path.parent.parent
    helper = api.load_module(
        "fine_CDS_target_component",
        ROOT / "scripts/validation/diagnose_saved_hessian.py",
    )
    coarse = helper.cds_for_state(state, 1e-4)
    component = state / "coarse_CDS_for_validation.npy"
    np.save(component, coarse)
    indices = list(
        map(int, np.unravel_index(abs(coarse - coarse.T).argmax(), coarse.shape))
    )
    matrix = np.load(path.with_suffix("") / "hessian_cartesian.npy")
    checks, evidence = worst_column_checks(
        state, matrix, indices, native.NATIVE_THERMAL
    )
    evidence[str(component)] = api.digest(component)
    record = json.loads(path.read_text())
    record["execution_derivatives"] = dict(
        cds_step_Bohr=1e-5,
        qualification_only=False,
        scientific_model_changed=False,
        qualification_sha256=derivative_certificate(native.NATIVE_THERMAL, 1e-5),
        source_geometry_sha256=api.digest(geometry),
        worst_contact_indices=indices,
        worst_contact_checks=checks,
        target_evidence_sha256=evidence,
        execution_helper_sha256=api.digest(__file__),
        finite_CDS_module_sha256=api.digest(ROOT / "src/redox/qm/cds_finite.py"),
    )
    atomic_write_json(path, record)


def resume_native(gid, state_name, state, spec, manifest_row):
    from redox.workflows.resources import verify_gate_owner

    verify_gate_owner()
    from redox.qm import native_thermal as native, projected_thermal as projected
    from redox.qm.hessian_memory import derivative_certificate
    from gpu4pyscf.solvent.hessian import smd as smd_hessian
    from ase.io import read

    state = Path(state)
    matrix, indices = reconstructed_matrix(spec)
    if api.digest(state / "opt.xyz") != spec["source_geometry_sha256"]:
        raise ValueError("Cached Hessian cannot be transferred to another geometry")
    coordinates = read(state / "opt.xyz").positions / 0.52917721067
    n = len(coordinates)
    certificate = derivative_certificate(native.NATIVE_THERMAL, 1e-5)
    previous_kernel = smd_hessian.WithSolventHess.kernel
    previous_projected = projected.ensure_projected_thermal
    previous_manifest = native.read_manifest
    calls = []

    def cached(hessian, *args, **kwargs):
        mf = hessian.base
        if (
            mf.mol.charge != 0
            or mf.mol.spin != 0
            or mf.xc != "r2scan"
            or mf.disp != "d4"
            or mf.grids.level != 3
            or hessian.grid_response is not True
            or mf.with_solvent.solvent != "acetonitrile"
            or str(mf.mol.basis).lower() != "def2-svp"
            or mf.conv_tol != 1e-12
            or mf.conv_tol_grad != 1e-8
            or np.max(abs(mf.mol.atom_coords(unit="Bohr") - coordinates)) > 1e-9
        ):
            raise ValueError(
                "Cached Hessian requested under different physical settings"
            )
        calls.append(True)
        hessian.de = matrix.reshape(n, 3, n, 3).transpose(0, 2, 1, 3).copy()
        return hessian.de

    def seal(path, geometry):
        if len(calls) != 1:
            raise ValueError("Saved Hessian was not used exactly once")
        checks, evidence = worst_column_checks(
            state, matrix, indices, native.NATIVE_THERMAL
        )
        raw = json.loads(Path(path).read_text())
        final = Path(path).with_suffix("") / "hessian_cartesian.npy"
        if not np.array_equal(np.load(final), matrix):
            raise ValueError("Native saved Hessian changed the cached matrix")
        raw["execution_hessian_reconstruction"] = dict(
            spec=spec,
            module_sha256=api.digest(__file__),
            source_geometry_sha256=api.digest(geometry),
            final_hessian=str(final),
            final_hessian_sha256=api.digest(final),
            original_electronic_PCM_hessian_reused=True,
            CDS_only_replaced=True,
        )
        raw["execution_derivatives"] = dict(
            cds_step_Bohr=1e-5,
            qualification_only=False,
            scientific_model_changed=False,
            qualification_sha256=certificate,
            source_geometry_sha256=api.digest(geometry),
            worst_contact_indices=indices,
            worst_contact_checks=checks,
            target_evidence_sha256=evidence,
        )
        raw["execution_bounded_pcm"] = dict(
            adapter_sha256=api.digest(ROOT / "src/redox/qm/pcm_bounded.py"),
            driver_sha256=api.digest(
                ROOT / "scripts/pipeline/resume_bounded_polymer_thermal.py"
            ),
            replay_evidence_sha256=spec["memory_replay_evidence_sha256"],
            source_geometry_sha256=api.digest(geometry),
            target_evidence_sha256=spec["input_sha256"],
            scientific_model_changed=False,
            surface_grid_changed=False,
            mode="cached_fused_PCM_Hessian_with_qualified_CDS_replacement",
        )
        atomic_write_json(path, raw)
        return previous_projected(path, geometry)

    smd_hessian.WithSolventHess.kernel = cached
    projected.ensure_projected_thermal = seal
    native.read_manifest = lambda: [manifest_row]
    try:
        return native.ensure_native_thermal(
            gid, state_name, state_dir=state, memory_cap_gib=24
        )
    finally:
        smd_hessian.WithSolventHess.kernel = previous_kernel
        projected.ensure_projected_thermal = previous_projected
        native.read_manifest = previous_manifest
