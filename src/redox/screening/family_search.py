"""Enumerate a declared family search space without registering or running QM.

Enumerated graphs and all symmetry-equivalent recipes are retained. Cheap
descriptors are bookkeeping, not predicted solution performance or a ranking.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

from rdkit import Chem

from redox.build.contracts import validate_group
from redox.core.common import ROOT, load_config, read_manifest
from redox.core.io import atomic_text_writer, atomic_write_json
from redox.screening.discovery import describe, precursor_from_frag
from redox.screening.expansion import registry_definitions


def ring_variants(fragment, substituents):
    molecule = Chem.MolFromSmiles(fragment)
    if molecule is None:
        raise ValueError("Invalid family template")
    sites = [a.GetIdx() for a in molecule.GetAtoms()
             if a.GetAtomicNum() == 6 and a.GetIsAromatic() and a.GetTotalNumHs() == 1]
    for label, substituent in substituents.items():
        if substituent is None:
            yield label, None, Chem.MolToSmiles(molecule)
            continue
        branch = Chem.MolFromSmiles(substituent)
        if branch is None:
            raise ValueError("Invalid ring substituent")
        for site in sites:
            product = Chem.RWMol(Chem.CombineMols(molecule, branch))
            product.AddBond(site, molecule.GetNumAtoms(), Chem.BondType.SINGLE)
            Chem.SanitizeMol(product)
            yield label, site, Chem.MolToSmiles(product)


def enumerate_space(config=None):
    config = load_config("family_search") if config is None else config
    unique = {}
    for core, spec in config.CORES.items():
        for anchor, anchor_fragment in config.ANCHORS.items():
            for cap, cap_fragment in config.CAPS.items():
                fragment = spec["template"].format(anchor=anchor_fragment, cap=cap_fragment)
                for ring, site, variant in ring_variants(fragment, config.RING_SUBSTITUENTS):
                    recipe = dict(core=core, anchor=anchor, cap=cap, ring=ring, template_atom=site)
                    group = dict(id="pending", name=f"{core} / {anchor} / {cap} / {ring}",
                                 family=spec["family"], frag=variant, states=spec["states"])
                    molecule = validate_group(config.SCAFFOLD, group)
                    smiles = Chem.MolToSmiles(molecule, isomericSmiles=True)
                    key = Chem.MolToInchiKey(molecule)
                    if key in unique:
                        if unique[key]["grafted_smiles"] != smiles:
                            raise RuntimeError("Different exact graphs share a search identity")
                        unique[key]["recipes"].append(recipe)
                        continue
                    group["id"] = "search_" + core + "_" + hashlib.sha256(smiles.encode()).hexdigest()[:12]
                    unique[key] = dict(group=group, inchikey=key, grafted_smiles=smiles,
                                       recipes=[recipe], core=core)
    return list(unique.values())


def prepare(output):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    config = load_config("family_search")
    known = {Chem.MolToInchiKey(validate_group(scaffold, group)) for scaffold, group in registry_definitions()}
    known.update(Chem.MolToInchiKey(Chem.MolFromSmiles(row["smiles"])) for row in read_manifest())
    faraday = load_config("electrolyte").FARADAY
    records, rows = enumerate_space(config), []
    for record in records:
        group = record["group"]
        group["provenance"] = dict(type="systematic family enumeration", exact_grafted_experiment=False,
                                    production_qm=False, route_status=config.ROUTE_STATUS)
        row = dict(id=group["id"], core=record["core"], inchikey=record["inchikey"],
                   grafted_smiles=record["grafted_smiles"], precursor_smiles=precursor_from_frag(group["frag"]),
                   coverage="already registered" if record["inchikey"] in known else "new exact graft",
                   recipes=json.dumps(record["recipes"], sort_keys=True),
                   **describe(record["grafted_smiles"], 2, faraday),
                   usable_electrons_status="not established; declared two-electron hypothesis",
                   solution_properties_status="not evaluated in this enumeration",
                   route_status=config.ROUTE_STATUS)
        rows.append(row)
    atomic_write_json(output / "structures.json", records)
    with atomic_text_writer(output / "search_space.csv", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    counts = Counter((r["core"], r["coverage"]) for r in rows)
    summary = dict(exact_grafts=len(rows), recipes=sum(len(r["recipes"]) for r in records),
                   by_family={core:dict(total=sum(r["core"] == core for r in rows),
                                         new=counts[core, "new exact graft"],
                                         registered=counts[core, "already registered"]) for core in config.CORES},
                   axes=dict(cores=list(config.CORES), caps=list(config.CAPS), anchors=list(config.ANCHORS),
                             ring_substituents=list(config.RING_SUBSTITUENTS)),
                   config=str(ROOT / "config/family_search.py"),
                   config_sha256=hashlib.sha256((ROOT / "config/family_search.py").read_bytes()).hexdigest(),
                   quantum_jobs_submitted=0, registered_candidates=0,
                   selection_status="pool prepared; source/route review and family-balanced screening selection pending")
    atomic_write_json(output / "summary.json", summary)
    return summary


def select_comparisons(records, known=frozenset(), panel="initial"):
    """Select declared contrasts per core, with an unchanged ethyl control.

    Equivalent ring products retain the lowest template atom index, so ring
    comparisons use the same core site across panels. This is coverage selection, not a
    prediction of electrochemical performance.
    """
    selected = []
    config = load_config("family_search")
    if panel not in config.COMPARISON_PANELS:
        raise ValueError("Unknown comparison panel: " + panel)
    contrasts = config.COMPARISON_PANELS[panel]
    for core in config.CORES:
        for name, anchor, cap, ring in contrasts:
            matches = [(recipe, record) for record in records
                       for recipe in record["recipes"] if (recipe["core"], recipe["anchor"], recipe["cap"], recipe["ring"])
                       == (core, anchor, cap, ring)]
            if not matches:
                raise ValueError(f"No new candidate for {core}/{name}; do not silently substitute another contrast")
            recipe, record = min(matches, key=lambda pair: (
                -1 if pair[0]["template_atom"] is None else pair[0]["template_atom"], pair[1]["inchikey"]))
            if record["inchikey"] in known:
                raise ValueError(f"No new candidate for {core}/{name}; declared site already registered, do not substitute another positional isomer")
            controls = [r for r in records if any((p["core"], p["anchor"], p["cap"], p["ring"])
                        == (core, "direct_N", "ethyl", "H") for p in r["recipes"])]
            if len(controls) != 1:
                raise ValueError(f"Unique ethyl control missing for {core}")
            selected.append(dict(record=record, contrast=name, recipe=recipe, control=controls[0]))
    if len({s["record"]["inchikey"] for s in selected}) != len(selected):
        raise ValueError("Contrasts select the same exact graft")
    return selected


def prepare_comparisons(output, pool, panel="initial"):
    output, pool = Path(output).resolve(), Path(pool).resolve()
    records_path = pool / "structures.json"
    records_bytes = records_path.read_bytes()
    registered = {Chem.MolToInchiKey(validate_group(scaffold, group)): group["id"]
                  for scaffold, group in registry_definitions()}
    selected = select_comparisons(json.loads(records_bytes), set(registered), panel)
    proposals, rows = [], []
    for item in selected:
        record, control, recipe = item["record"], item["control"], item["recipe"]
        if control["inchikey"] not in registered:
            raise ValueError("Comparison control must already be registered")
        group = dict(record["group"], comparison_baseline=registered[control["inchikey"]],
                     search_axes=recipe, comparison_contrast=item["contrast"],
                     hypothesis=f"Measure {item['contrast']} on {record['core']} against the unchanged ethyl graft; no performance improvement assumed.")
        group["provenance"] = dict(group.get("provenance", {}), source_pool=str(records_path),
            source_pool_sha256=hashlib.sha256(records_bytes).hexdigest(),
            type="family-balanced computational comparison", exact_grafted_experiment=False,
            production_qm=False, parent_metrics_inherited=False)
        molecule = validate_group(load_config("family_search").SCAFFOLD, group)
        proposals.append(group)
        rows.append(dict(id=group["id"], core=record["core"], contrast=item["contrast"],
                         baseline=group["comparison_baseline"], inchikey=record["inchikey"],
                         total_atoms=Chem.AddHs(molecule).GetNumAtoms(),
                         selected_template_atom=recipe["template_atom"],
                         properties_status="UMA screening pending; solution properties uncomputed"))
    output.mkdir(parents=True, exist_ok=False)
    atomic_write_json(output / "proposals.json", proposals)
    with atomic_text_writer(output / "comparison_matrix.csv", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    report = dict(proposals=len(proposals), states=sum(len(g["states"]) for g in proposals),
                  proposals_sha256=hashlib.sha256((output / "proposals.json").read_bytes()).hexdigest(),
                  source_pool=str(records_path), source_pool_sha256=hashlib.sha256(records_bytes).hexdigest(),
                  panel=panel, selection_rule="declared contrasts per core; identical ethyl control; matched ring site across panels",
                  planned_first_stage="UMA gas preoptimization after new scope/resource qualification",
                  quantum_jobs_submitted=0)
    atomic_write_json(output / "preparation.json", report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--select-from", type=Path, help="prepare the declared comparison panel from a preserved pool")
    parser.add_argument("--panel", choices=tuple(load_config('family_search').COMPARISON_PANELS), default="initial",
                        help="declared contrast panel; donor_rings adds OMe/Me at the initial F/CN site")
    args = parser.parse_args(argv)
    if args.panel != 'initial' and not args.select_from:
        parser.error('--panel requires --select-from')
    print(json.dumps(prepare_comparisons(args.out, args.select_from, args.panel) if args.select_from else prepare(args.out), indent=2))


if __name__ == "__main__":
    main()
