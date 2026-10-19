"""Matched-comparison candidate expansion with exact-graph registration.

All proposed electrochemical behavior and exact graft routes remain hypotheses.
Preparation writes definitions and bookkeeping descriptors, never inherited parent
measurements. Registration appends validated new definitions without rewriting old ones.
"""
from __future__ import annotations

import argparse
import csv
from datetime import date
import fcntl
import hashlib
import json
from pathlib import Path
import pprint
import runpy

from rdkit import Chem

from redox.build.contracts import validate_group
from redox.core.common import ROOT, load_config, read_manifest
from redox.core.io import atomic_text_writer, atomic_write_json
from redox.screening.discovery import describe, precursor_from_frag

SCAFFOLD = "Cc1ccc(C[*:1])cc1"
REGISTRIES = ("redox_groups", "starting_candidates", "merrifield_multielectron",
              "discovered_candidates", "custom_candidates")
REDUCTIONS = [("neu", 0, 1, 0), ("red1", -1, 2, -1), ("red2", -2, 1, -2)]
OXIDATIONS = [("neu", 0, 1, 0), ("ox1", 1, 2, 1), ("ox2", 2, 1, 2)]
SOURCE_SCOPE = {
    "imide": dict(url="https://pubs.acs.org/doi/10.1021/acsaem.1c01490",
                  scope="Family rationale for imide multiredox screening; no measurements or synthesis proof for these grafts."),
    "ptz": dict(url="https://pubs.rsc.org/en/content/articlelanding/2017/ta/c7ta05883g",
                scope="Methoxy substitution stabilizes a particular parent dication. Methyl/fluoro grafts are contrasting hypotheses, not established two-electron donors."),
}


def registry_definitions():
    return [(module.SCAFFOLD, group) for name in REGISTRIES
            for module in [load_config(name)] for group in module.GROUPS]


def designs():
    groups = []
    templates = {"pmdi": "[*:1]n1c(=O)c2cc3c(=O)n({cap})c(=O)c3cc2c1=O",
                 "ndi": "[*:1]N1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N({cap})C3=O"}
    for core, template in templates.items():
        for cap, fragment, baseline, hypothesis in (
            ("cyclopropyl", "C5CC5", f"design_{core}_isopropyl",
             "Test cap topology and conformer sensitivity against isopropyl; cyclization also removes two H atoms, so this is not a pure rigidity perturbation."),
            ("2fluoroethyl", "CCF", f"design_{core}_ethyl",
             "Test replacing one terminal ethyl H by F against the ethyl cap; evaluate E1/E2, wave spacing and mass cost without assuming improved solubility."),
        ):
            groups.append(dict(id=f"hyp_{core}_{cap}", name=f"N-(4-methylbenzyl)-N'-{cap} {core.upper()}",
                               family="imide (n-type)", frag=template.format(cap=fragment), states=REDUCTIONS,
                               comparison_baseline=baseline, hypothesis=hypothesis, source_key="imide"))
    for name, substituent in (("dimethyl", "C"), ("difluoro", "F")):
        groups.append(dict(id=f"hyp_ptz_37{name}", name=f"N-(4-methylbenzyl)-3,7-{name}phenothiazine",
                           family="amine (p-type)",
                           frag=f"[*:1]N1c2ccc({substituent})cc2Sc2cc({substituent})ccc21", states=OXIDATIONS,
                           comparison_baseline="exp_ptz_37dimethoxy", source_key="ptz",
                           hypothesis="Test ring substitution against the dimethoxy graft. The second oxidation is an explicit stress test: accept two accessible electrons only if both waves, integrity and physical QC pass."))
    for group in groups:
        source = SOURCE_SCOPE[group.pop("source_key")]
        group["provenance"] = dict(type="matched-comparison computational hypothesis", source=source["url"],
                                   evidence_scope=source["scope"], exact_grafted_experiment=False,
                                   production_qm=False, direct_merrifield_route_verified=False,
                                   registration_stage="exact structure/state definition; properties pending")
        group["chemistry_notes"] = (group["hypothesis"] + " Exact graft/precursor synthesis is unverified; "
                                     "parent measurements are not graft properties. Imide mono-N-H precursor isolation "
                                     "or selective phenothiazine N-alkylation requires separate route review.")
    return groups


def prepare(output):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    known = {Chem.MolToInchiKey(validate_group(scaffold, group)) for scaffold, group in registry_definitions()}
    known.update(Chem.MolToInchiKey(Chem.MolFromSmiles(row["smiles"])) for row in read_manifest())
    groups, rows, seen = [], [], set()
    faraday = load_config("electrolyte").FARADAY
    for group in designs():
        molecule = validate_group(SCAFFOLD, group)
        key = Chem.MolToInchiKey(molecule)
        if key in known:
            continue
        if key in seen:
            raise RuntimeError("Duplicate proposed exact graft")
        seen.add(key); groups.append(group)
        smiles = Chem.MolToSmiles(molecule)
        rows.append(dict(id=group["id"], baseline=group["comparison_baseline"], inchikey=key,
                         grafted_smiles=smiles, precursor_smiles=precursor_from_frag(group["frag"]),
                         states=json.dumps(group["states"]), hypothesis=group["hypothesis"],
                         exact_grafted_experiment=False, properties_status="not computed",
                         **describe(smiles, 2, faraday)))
    atomic_write_json(output / "proposals.json", groups)
    with atomic_text_writer(output / "designs.csv", newline="") as stream:
        columns = list(rows[0]) if rows else ["id", "baseline", "properties_status"]
        writer = csv.DictWriter(stream, fieldnames=columns); writer.writeheader(); writer.writerows(rows)
    atomic_write_json(output / "preparation.json", dict(proposals=len(groups), states=sum(len(g["states"]) for g in groups),
                      proposals_sha256=hashlib.sha256((output / "proposals.json").read_bytes()).hexdigest(),
                      source_scopes=SOURCE_SCOPE, quantum_jobs_submitted=0))
    return groups


def register(output, config=None):
    output = Path(output).resolve()
    config = Path(config) if config is not None else ROOT / "config/custom_candidates.py"
    content = (output / "proposals.json").read_bytes()
    summary = json.loads((output / "preparation.json").read_text())
    if hashlib.sha256(content).hexdigest() != summary["proposals_sha256"]:
        raise RuntimeError("Prepared proposals changed before registration")
    proposals = json.loads(content)
    locks = ROOT / "logs/config_locks"; locks.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256(str(config.resolve()).encode()).hexdigest()[:16]
    with (locks / f"{key}.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        original = config.read_bytes()
        current = runpy.run_path(str(config))
        existing = registry_definitions() + [(current["SCAFFOLD"], group) for group in current["GROUPS"]]
        ids = {group["id"]: Chem.MolToInchiKey(validate_group(scaffold, group)) for scaffold, group in existing}
        known = set(ids.values())
        pending = []
        for group in proposals:
            identity = Chem.MolToInchiKey(validate_group(SCAFFOLD, group))
            if group["id"] in ids and ids[group["id"]] != identity:
                raise ValueError(f"{group['id']}: ID belongs to another exact structure")
            if identity in known:
                continue
            pending.append(group); ids[group["id"]] = identity; known.add(identity)
        if pending:
            before = hashlib.sha256(original).hexdigest()
            archive = output / "config_before" / f"{before}.py"
            archive.parent.mkdir(exist_ok=True); archive.write_bytes(original)
            addition = f"\n\n# Matched-comparison hypotheses registered {date.today().isoformat()}; exact-graft properties remain pending.\nGROUPS += " + pprint.pformat(pending, sort_dicts=False, width=100) + "\n"
            if config.read_bytes() != original:
                raise RuntimeError("Candidate registry changed during preparation")
            with atomic_text_writer(config) as stream:
                stream.write(original.decode() + addition)
        result = dict(registered=[group["id"] for group in pending], already_present=len(proposals)-len(pending),
                      config=str(config.resolve()), config_sha256=hashlib.sha256(config.read_bytes()).hexdigest(),
                      quantum_jobs_submitted=0)
        atomic_write_json(output / "registration.json", result)
        return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--register", action="store_true", help="append prepared validated hypotheses to custom_candidates.py")
    args = parser.parse_args(argv)
    if not args.out.exists():
        groups = prepare(args.out)
        print(f"Prepared {len(groups)} exact-graph proposals")
    elif not args.register:
        parser.error("output exists; use another directory or --register the existing preparation")
    if args.register:
        print(json.dumps(register(args.out)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
