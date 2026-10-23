#!/usr/bin/env python
"""Export primary paper records and source without running or relabeling QM.

Objects retain their original bytes. The portable example verifies every G
before reconstructing benchmark, route and finite-conformer arithmetic.
"""
import argparse
import hashlib
import io
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import time

from redox.core.common import ROOT, DFT, read_manifest
from redox.core.protocol import load_record, record_path
from redox.qm.native_thermal import cache_path, load_native_thermal
from redox.workflows.checkpoints import verify_state_checkpoint
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def extract_source(archive, destination):
    """Git source exports may contain only ordinary files/directories."""
    with tarfile.open(fileobj=io.BytesIO(archive)) as stream:
        for member in stream.getmembers():
            path = Path(member.name)
            if path.is_absolute() or ".." in path.parts or not (member.isfile() or member.isdir()):
                raise ValueError("Unsafe source archive member: " + member.name)
        stream.extractall(destination, filter="data")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--code-ref", default="HEAD")
    parser.add_argument("--max-GiB", type=float, default=2.)
    args = parser.parse_args()
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    if args.out.exists():
        raise ValueError("Preserve the existing package; choose a new output directory")
    code_ref = subprocess.check_output(["git", "rev-parse", args.code_ref], cwd=ROOT, text=True).strip()
    producer_paths = [ROOT / path for path in (
        "scripts/release/prepare_review_package.py", "scripts/release/reproduce_review_package.py",
        "scripts/release/review_package_hessians.py", "scripts/release/review_package_auxiliary.py",
        "scripts/validation/review_release_package.py", "tests/test_review_package_auxiliary.py")]
    producer_hashes = {str(p): digest(p) for p in producer_paths}
    pointer_path = ROOT / "manuscript/evidence/current.json"
    pointer = json.loads(pointer_path.read_text())
    epoch = (ROOT / pointer["manifest"]).parent
    frozen = json.loads((epoch / "manifest.json").read_text())
    if any(digest(ROOT / p) != h for p, h in frozen["snapshot_sha256"].items()):
        raise ValueError("Frozen evidence differs from its manifest")
    sources = {str(pointer_path): digest(pointer_path), str(epoch / "manifest.json"): digest(epoch / "manifest.json")}
    load = lambda name: json.loads((epoch / name).read_text())
    catalog, aliases = {}, {}
    def state(alias, directory, smiles, expected, alternative_surface=None):
        directory = Path(directory).resolve()
        key = "s_" + hashlib.sha256(str(directory).encode()).hexdigest()[:20]
        if key in catalog:
            if abs(catalog[key]["expected_G_eV"] - expected) > 1e-8:
                raise ValueError("Conflicting checked G values for one source geometry")
        else:
            catalog[key] = dict(id=key, source=directory, smiles=smiles, expected_G_eV=expected, aliases=[], alternative_surface=alternative_surface)
        if alias in aliases and aliases[alias] != key:
            raise ValueError("Conflicting state alias")
        aliases[alias] = key
        catalog[key]["aliases"].append(alias)
        return key
    def receipt(alias, path):
        report = verify_state_checkpoint(path)
        sources[str(path)] = digest(path)
        row = report["task"]["manifest_row"]
        return state(alias, report["state_dir"], row["smiles"], report["G_eV"])
    manifest = {(r["id"], r["state"]): r for r in read_manifest()}
    benchmark = load("benchmark.json")
    review_states = next(Path(p) for p in benchmark["review_snapshot_sha256"] if Path(p).name == "states.json")
    if digest(review_states) != benchmark["review_snapshot_sha256"][str(review_states)]:
        raise ValueError("Benchmark state snapshot changed")
    sources[str(review_states)] = digest(review_states)
    for label, row in json.loads(review_states.read_text()).items():
        gid, redox_state = label.split(":")
        if row["status"] != "ready":
            raise ValueError("Unready benchmark state in primary panel")
        state("benchmark:" + label, DFT / gid / redox_state, manifest[gid, redox_state]["smiles"], row["G_eV"])
    reference = ROOT / "results/grafting/campaigns/2026-10-19_mecn_reference"
    for label, path in load("grafting_status.json")["components"].items():
        receipt("component:" + label, Path(path))
    bases = ROOT / "results/grafting/campaigns/2026-10-19_base_assisted_pilot"
    for task in json.loads((bases / "plan.json").read_text())["tasks"]:
        receipt("base:" + task["key"], bases / "done" / (task["key"] + ".json"))
    grafting_path = ROOT / "results/grafting/current.json"
    grafting = json.loads(grafting_path.read_text())["candidates"]
    sources[str(grafting_path)] = digest(grafting_path)
    reactions = []
    complete = [r for r in grafting if r["deltaG_reference_eV"] is not None]
    if len(complete) != load("grafting_status.json")["completed_grafting_deltaG"]:
        raise ValueError("Completed grafting coverage differs from the frozen paper")
    for row in complete:
        raw = json.loads((Path(row["product_state_path"]) / "result.json").read_text())
        from redox.core.common import free_energy, read_result
        directory = Path(row["product_state_path"])
        if directory.is_relative_to(DFT):
            expected, quality = free_energy(read_result(directory.parent.name, directory.name))
            if quality != "ok" or expected is None:
                raise ValueError("Unqualified completed graft product: " + row["id"] + "; " + quality)
        else:
            product_receipt = reference / "done" / (raw["id"] + ".json")
            report = verify_state_checkpoint(product_receipt)
            if Path(report["state_dir"]).resolve() != directory.resolve():
                raise ValueError("Separate-workspace graft product belongs to another receipt")
            sources[str(product_receipt)] = digest(product_receipt)
            expected = report["G_eV"]
        product = state("product:" + row["id"], row["product_state_path"], row["product_smiles"], expected)
        precursor = aliases["component:" + row["precursor_id"]]
        site = aliases["component:merrifield_local_chloride"]
        byproduct = aliases["component:" + ("hydrogen_chloride" if row["reaction_class"] == "NH_OH_substitution" else "chloride")]
        reactions.append(dict(name="reference:" + row["id"], coefficients={product: 1, byproduct: 1, precursor: -1, site: -1},
            delta_n_solutes=0, source_kJ_mol=row["deltaG_reference_kJ_mol"]))
    for row in load("base_assisted.json")["candidates"]:
        if row["reference_deltaG_kJ_mol"] is None:
            continue
        product, precursor, site = (aliases["product:" + row["id"]], aliases["component:" + row["reaction"]["precursor_id"]], aliases["component:merrifield_local_chloride"])
        for base, modes in row["base_models"].items():
            for kind, expected in modes.items():
                if expected is None:
                    continue
                coefficients = {product: 1, precursor: -1, site: -1, aliases["base:" + base + "_neutral"]: -1}
                if kind == "free_ions_kJ_mol":
                    coefficients.update({aliases["base:" + base + "_protonated"]: 1, aliases["component:chloride"]: 1})
                elif kind == "contact_pair_kJ_mol":
                    coefficients[aliases["base:" + base + "_chloride_pair"]] = 1
                else:
                    raise ValueError("Unknown base representation")
                reactions.append(dict(name=f"{base}:{kind}:{row['id']}", coefficients=coefficients,
                    delta_n_solutes=sum(coefficients.values()), source_kJ_mol=expected))
    for row in load("exact_case_join.json")["cases"]:
        gid = row["candidate_id"]
        for redox_state, g in row["state_G_eV"].items():
            state(f"exact:{gid}:{redox_state}", DFT / gid / redox_state, manifest[gid, redox_state]["smiles"], g)
    finite = load("finite_ensemble.json")
    cases = {row["case"]["key"]: row["case"] for row in load("finite_ensemble_states.json")["cases"]}
    for row in finite["states"]:
        graph = row["key"].split("__pose")[0]
        state("polymer:" + row["key"], row["state_dir"], cases[graph]["smiles"], row["G_eV"])
    for row in load("proton_transfer.json")["states"]:
        receipt("relative_basicity:" + row["id"], Path(row["source_receipt"]))
    spec = importlib.util.spec_from_file_location("auxiliary_export", Path(__file__).with_name("review_package_auxiliary.py"))
    auxiliary = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(auxiliary)
    auxiliary_data = auxiliary.add_auxiliary_states(ROOT, load, state, sources, digest)
    # First inventory and verify every source; no bulk output exists on failure.
    files, exported = {}, []
    def artifact(path):
        path = Path(path)
        h = digest(path)
        sources[str(path)] = h
        relative = "objects/" + h[:2] + "/" + h + path.suffix
        files[relative] = dict(source=path, sha256=h, bytes=path.stat().st_size)
        return relative
    for item in catalog.values():
        sd = item["source"]
        raw = json.loads((sd / "result.json").read_text())
        q, mult = raw["charge"], raw["mult"]
        active = load_record(sd, q, mult)
        if active is None:
            raise ValueError("Unqualified electronic record: " + str(sd))
        paths = dict(geometry=sd / "opt.xyz", result=sd / "result.json", energy=record_path(sd, q, mult))
        if (sd / "monatomic_thermal.json").exists():
            partition = json.loads((sd / "monatomic_thermal.json").read_text())
            thermal, kind = partition["G_thermal_eV"], "monatomic"
            paths["partition"] = sd / "monatomic_thermal.json"
        elif item["alternative_surface"]:
            from redox.workflows import polymer_campaign as api
            independent = api.load_module("basis_export_partition", ROOT / "scripts/validation/audit/native_partition.py")
            native_path = Path(item["alternative_surface"])
            auxiliary.validate_alternative_surface(json.loads(native_path.read_text()))
            partition = independent.read_partition(native_path, sd / "opt.xyz")
            if partition is None:
                raise ValueError("Alternative basis partition failed independent reconstruction")
            paths.update(native=native_path, partition=Path(partition["path"]),
                hessian=native_path.with_suffix("") / "hessian_cartesian.npy")
            thermal, kind = partition["g_thermal_total_eV"], "projected_DFT_alternative_basis_control"
        else:
            native = load_native_thermal(sd, q, mult)
            if native is None:
                raise ValueError("Unqualified native partition: " + str(sd))
            native_path, _ = cache_path(sd, q, mult)
            paths.update(native=native_path, partition=Path(native["thermal_postprocessing_path"]),
                         hessian=native_path.with_suffix("") / "hessian_cartesian.npy")
            diagnostics = native_path.with_suffix("") / "hessian_diagnostics.json"
            if diagnostics.exists():
                paths["hessian_diagnostics"] = diagnostics
            thermal, kind = native["g_thermal_total_eV"], "projected_DFT"
        total = active["e_smd_eV"] + thermal
        if abs(total - item["expected_G_eV"]) > 1e-8:
            raise ValueError("Current raw G differs from frozen evidence: " + str(sd))
        roles = {}
        for role, path in paths.items():
            roles[role] = artifact(path)
        exported.append(dict(id=item["id"], aliases=item["aliases"], smiles=item["smiles"], charge=q, mult=mult,
            partition_kind=kind, files=roles, checked_G_eV=total, original_state_directory=str(sd)))
    auxiliary_data["proof_files"] = [artifact(path) for path in auxiliary_data.pop("proof_sources")]
    auxiliary_data["electronic"] = auxiliary.electronic_inventory(ROOT, load, artifact, digest)
    size = sum(r["bytes"] for r in files.values())
    if size > args.max_GiB * 2**30 or shutil.disk_usage(next(p for p in args.out.parents if p.exists())).free < 2 * size + 2**30:
        raise ValueError("Package inventory exceeds the declared disk budget")
    if any(digest(p) != h for p, h in sources.items()):
        raise ValueError("A source changed during inventory")
    args.out.mkdir(parents=True)
    (args.out / "data").mkdir()
    for relative, row in files.items():
        target = args.out / "data" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(row["source"], target)
        if digest(target) != row["sha256"]:
            raise ValueError("Copied object differs")
    src_paths = ["src", "scripts", "tests", "config", "docs", "environments", "examples", "LICENSE", "pyproject.toml", "setup_env.sh",
        "manuscript/build.py", "manuscript/check_claims.py", "manuscript/freeze_evidence.py", "manuscript/manuscript.tex",
        "manuscript/supporting_information.tex", "manuscript/references.bib", "manuscript/evidence/current.json",
        str(epoch.relative_to(ROOT)), pointer["figures"]]
    archive = subprocess.check_output(["git", "archive", code_ref, *src_paths], cwd=ROOT)
    extract_source(archive, args.out / "code")
    overlays = {}
    for source in producer_paths:
        target = args.out / "code" / source.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        overlays[str(target.relative_to(args.out))] = digest(target)
    potential = dict(oxidized=aliases["product:design_pmdi_ethyl"], reduced=aliases["exact:design_pmdi_ethyl:red1"],
        reference_oxidized=aliases["benchmark:ferrocene:ox"], reference_reduced=aliases["benchmark:ferrocene:neu"])
    data_manifest = dict(schema=1, temperature_K=298.15, states=exported, potential=potential, reactions=reactions,
        source_energy_conversion_kJ_mol_per_eV=96.485, SI_energy_conversion_kJ_mol_per_eV=96.48533212331002,
        aliases=aliases, auxiliary=auxiliary_data, files_sha256={p:r["sha256"] for p,r in files.items()})
    (args.out / "data/manifest.json").write_text(json.dumps(data_manifest, indent=2) + "\n")
    (args.out / "README.md").write_text(
        "# Primary redox-polymer review package\n\n"
        "Run from any directory with standard Python 3.10 or newer:\n\n"
        "`python3 -I -S PACKAGE/code/scripts/release/reproduce_review_package.py PACKAGE`\n\n"
        "The offline replay checks all packaged bytes, state/geometry identities, "
        "atom/charge-balanced sums, 39 benchmark endpoints, 67 reference reactions, "
        "20 base-assisted values, two relative-basicity controls and the twelve-state "
        "finite polymer partition. Additional records support five alternative-basis "
        "controls, two interrupted/resumed comparisons, 26 stationary electronic models "
        "and six frozen-fragment energies. It also replays frozen manuscript claim checks. "
        "Original absolute paths are descriptive provenance and are never opened.\n\n"
        "With NumPy and ASE installed, an additional CPU check reconstructs modes and "
        "nuclear sums from the raw Hessians:\n\n"
        "`python3 PACKAGE/code/scripts/release/review_package_hessians.py PACKAGE`\n\n"
        "Set OMP_NUM_THREADS, OMP_THREAD_LIMIT, MKL_NUM_THREADS, OPENBLAS_NUM_THREADS, "
        "NUMEXPR_NUM_THREADS and XTB_THREADS to 1 for that check.\n\n"
        "`code/` records a Git source commit plus the explicitly hashed release-tool "
        "overlay. `data/objects/` contains original raw records, geometries and Hessians "
        "addressed by content hash. `data/manifest.json` joins aliases and defines sums. "
        "`package_manifest.json` binds every packaged file. Recorded conda/pip versions "
        "are in `code/environments/`; the full CUDA environment is not freshly qualified.\n\n"
        "This is a prepared scientific-record review package, not a DOI deposition. "
        "Execution/numerical diagnostic artifact closure and a fresh GPU restoration "
        "remain separate checks. Broader conformer convergence, exact PMDI occupancy and material "
        "performance remain unresolved. Code uses the existing MIT license; a deposition "
        "license for data and author metadata require confirmation before public release.\n")
    if any(digest(p) != h for p, h in producer_hashes.items()):
        raise ValueError("Release producer changed during preparation")
    info = dict(schema=1, status="prepared_pending_offline_replay", prepared=time.time(), code_commit=code_ref, original_checkout_root=str(ROOT),
        code_overlays_sha256=overlays, evidence_epoch="code/" + str(epoch.relative_to(ROOT)),
        raw_states=len(exported), benchmark_endpoints=39, reference_reactions=len(complete), base_assisted_values=len(reactions)-len(complete),
        finite_native_states=12, alternative_basis_states=len(auxiliary_data["basis"]), restart_states=len(auxiliary_data["restart"]),
        electronic_only_records=len(auxiliary_data["electronic"]["records"]), raw_objects=len(files), raw_object_bytes=size,
        original_sources_sha256=sources, producer_sha256=producer_hashes, code_license="MIT", data_deposition_license=None, public_archive_deposited=False,
        scope="Portable primary and auxiliary scientific records; GPU rerun/environment qualification and execution/numerical diagnostic artifact closure remain separate.",
        files_sha256={str(p.relative_to(args.out)):digest(p) for p in args.out.rglob("*") if p.is_file()})
    (args.out / "package_manifest.json").write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps({k:info[k] for k in ("status","raw_states","reference_reactions","base_assisted_values","raw_objects","raw_object_bytes")}))


if __name__ == "__main__":
    main()
