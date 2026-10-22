"""Stream exact family recipes into deduplicated, parallel descriptor chunks.

This extends the design/evaluation loop without registering molecules or
launching quantum calculations. Its cheap outputs cannot bypass the route gate.
"""

import argparse
from functools import cache
from importlib import metadata
import json
import os
from pathlib import Path

from redox.core.common import ROOT, load_config
from redox.core.io import atomic_write_json
from redox.workflows import chunked


def family_records(config=None):
    from rdkit import Chem
    from redox.build.contracts import validate_group
    from redox.screening.family_search import ring_variants
    from redox.screening.discovery import precursor_from_frag

    config = load_config("family_search") if config is None else config
    for core, spec in config.CORES.items():
        for anchor, anchor_fragment in config.ANCHORS.items():
            for cap, cap_fragment in config.CAPS.items():
                fragment = spec["template"].format(
                    anchor=anchor_fragment, cap=cap_fragment
                )
                for ring, site, variant in ring_variants(
                    fragment, config.RING_SUBSTITUENTS
                ):
                    group = dict(
                        id="streamed_design",
                        name=core,
                        family=spec["family"],
                        frag=variant,
                        states=spec["states"],
                    )
                    molecule = validate_group(config.SCAFFOLD, group)
                    product = Chem.MolToSmiles(molecule, isomericSmiles=True)
                    precursor = precursor_from_frag(variant)
                    fragment_graph = Chem.MolToSmiles(
                        Chem.MolFromSmiles(variant), isomericSmiles=True
                    )
                    identity = dict(
                        product_smiles=product,
                        precursor_smiles=precursor,
                        scaffold_smiles=config.SCAFFOLD,
                        attachment_fragment=fragment_graph,
                        new_bond="Join the neighbors of the two deleted [*:1] atoms",
                    )
                    yield dict(
                        identity=identity,
                        payload=dict(
                            product_smiles=product,
                            precursor_smiles=precursor,
                            attachment_fragment=fragment_graph,
                            states=spec["states"],
                            route_status=config.ROUTE_STATUS,
                            family=spec["family"],
                        ),
                        provenance=dict(
                            core=core,
                            anchor=anchor,
                            cap=cap,
                            ring=ring,
                            template_atom=site,
                        ),
                    )


def control_records(count):
    for index in range(count):
        yield dict(
            identity=dict(synthetic_index=index),
            payload=dict(index=index),
            provenance=dict(role="scheduler/storage control; no molecule"),
        )


def stage_spec(handler):
    if handler not in ("descriptors", "scheduler_control"):
        raise ValueError("Unregistered stage handler")
    model_files = {}
    if handler == "descriptors":
        from rdkit.Chem import RDConfig

        model_files = {
            str(p): chunked.file_hash(p)
            for p in (
                Path(RDConfig.RDContribDir) / "SA_Score/fpscores.pkl.gz",
                Path(RDConfig.RDContribDir) / "SA_Score/sascorer.py",
            )
        }
    return dict(
        name=handler,
        quality_scope=(
            "bookkeeping_descriptors_only"
            if handler == "descriptors"
            else "synthetic_scheduler_control_no_chemistry"
        ),
        protocol=dict(
            schema=1,
            handler=handler,
            source_sha256=chunked.file_hash(__file__),
            model_sha256=model_files,
            software={
                name: metadata.version(name) for name in ("rdkit", "numpy", "scipy")
            },
            dependency_sha256={
                str(p): chunked.file_hash(p)
                for p in (
                    ROOT / "src/redox/screening/discovery.py",
                    ROOT / "src/redox/build/contracts.py",
                    ROOT / "src/redox/screening/family_search.py",
                    ROOT / "src/redox/build/structures.py",
                    ROOT / "src/redox/properties/capacity_and_proxies.py",
                    ROOT / "config/family_search.py",
                    ROOT / "config/electrolyte.py",
                    ROOT / "config/project.json",
                )
            },
        ),
    )


@cache
def descriptor_tools():
    from redox.properties.capacity_and_proxies import _sascorer

    return _sascorer(), load_config("electrolyte").FARADAY


def descriptor(task):
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from redox.properties.capacity_and_proxies import _sa_prep

    payload = task["payload"]
    states = payload["states"]
    electrons = max(s[3] for s in states) - min(s[3] for s in states)
    scorer, faraday = descriptor_tools()
    molecule = Chem.MolFromSmiles(payload["product_smiles"])
    if molecule is None:
        raise ValueError("Invalid exact product graph")
    mass = Descriptors.MolWt(molecule)
    prepared, _, note = _sa_prep(payload["product_smiles"])
    sa = scorer.calculateScore(prepared) if prepared is not None else None
    result = dict(
        MW=round(mass, 3),
        n_e_nominal=electrons,
        capacity_nominal_mAh_g=round(electrons * faraday / (3.6 * mass), 1),
        SA_grafted=round(sa, 2) if sa is not None else None,
        sa_note=note,
        charge_resting=Chem.GetFormalCharge(molecule),
        heavy_atoms=molecule.GetNumHeavyAtoms(),
    )
    return dict(
        result,
        solution_properties_status="not computed",
        route_status=payload["route_status"],
        nominal_electron_count_is_hypothesis=True,
        synthesis_eligibility="not established",
    )


def scheduler_control(task):
    index = task["payload"]["index"]
    return dict(
        index=index, square=index * index, role="synthetic control; no chemistry"
    )


def handler_for(plan):
    name = plan["stage"]["name"]
    if plan["stage"] != stage_spec(name):
        raise ValueError("Registered handler, protocol or dependency source changed")
    if name == "descriptors":
        # Cache the scorer once; the legacy helper appends to sys.path on each call.
        scorer, _ = descriptor_tools()
        from rdkit import Chem

        scorer.calculateScore(Chem.MolFromSmiles("CC"))
        return descriptor
    return scheduler_control


def verify(campaign):
    campaign = Path(campaign)
    plan = chunked.read_plan(campaign)
    counts = dict(completed=0, failed=0, done_chunks=0)
    for part in plan["chunks"]:
        path = campaign / "done" / (part["key"] + ".json")
        if path.exists():
            record = chunked.verify_done(path, plan=plan)
            counts["completed"] += record["completed"]
            counts["failed"] += record["failed"]
            counts["done_chunks"] += 1
    result = dict(
        status="execution_complete"
        if counts["done_chunks"] == len(plan["chunks"])
        else "incomplete_execution",
        **counts,
        expected_tasks=plan["unique_tasks"],
        expected_chunks=len(plan["chunks"]),
        stage=plan["stage"],
        scientific_accuracy_established=False,
    )
    atomic_write_json(campaign / "summary.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action",
        choices=(
            "prepare-family",
            "prepare-control",
            "worker",
            "verify",
            "recover-local",
        ),
    )
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--count", type=int, default=1000000)
    parser.add_argument("--owner", type=Path)
    parser.add_argument("--checkpoint-every", type=int, default=100)
    args = parser.parse_args()
    if args.action.startswith("prepare"):
        from redox.workflows.resources import (
            allocation_guard,
            cpu_allocation,
            host_capacity,
        )

        with allocation_guard():
            host_capacity(1)
            os.sched_setaffinity(0, cpu_allocation(1))
        if args.count < 0:
            parser.error("Negative control size")
        real = args.action == "prepare-family"
        result = chunked.prepare(
            family_records() if real else control_records(args.count),
            args.campaign,
            stage_spec("descriptors" if real else "scheduler_control"),
            chunk_size=args.chunk_size,
        )
        print(
            json.dumps(
                {
                    k: result[k]
                    for k in (
                        "input_records",
                        "unique_tasks",
                        "duplicate_recipes",
                        "prepared_seconds",
                    )
                }
            )
        )
    elif args.action == "worker":
        plan = chunked.read_plan(args.campaign)
        chunked.run_worker(
            args.campaign, handler_for(plan), checkpoint_every=args.checkpoint_every
        )
    elif args.action == "recover-local":
        if args.owner is None:
            parser.error("recover-local requires --owner")
        print(chunked.recover_local(args.campaign, args.owner))
    else:
        print(json.dumps(verify(args.campaign), indent=2))


if __name__ == "__main__":
    main()
