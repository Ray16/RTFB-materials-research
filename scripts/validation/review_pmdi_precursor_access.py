#!/usr/bin/env python
"""Bind the PMDI precursor-access source review to its exact graph and route.

The patent's identity assignment is author-reported. Text and rendered scheme
inspection support the source binding; this is not a synthesis replication or
an experimental validation of the proposed polymer attachment.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import time

from pypdf import PdfReader
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from redox.core.common import ROOT
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--pdf", type=Path, default=ROOT / "calcs/literature/2026-10-23_precursor_access/EP0010440B1.pdf")
    args = parser.parse_args()
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    if sha(args.pdf) != "8e793ae73ad13b5d06102e0d203627fc6ff522db01def6935a7a586ff9660ada":
        raise ValueError("Primary patent differs from the inspected source")
    text = re.sub(r"\s+", " ", PdfReader(args.pdf).pages[9].extract_text())
    example = text.split("Example 19", 1)[1].split("Example 20", 1)[0]
    if "N-Ethyl pyromellitic diimide" not in example or "nuclear magnetic resonance" not in example:
        raise ValueError("Expected source identity/assignment passage missing")
    evidence_path = ROOT / "results/grafting/route_evidence.json"
    evidence = json.loads(evidence_path.read_text())
    row = evidence["candidates"]["design_pmdi_ethyl"]
    before_path = args.out / "before_route_evidence.json"
    before = json.loads(before_path.read_text())["candidates"]["design_pmdi_ethyl"]
    if row["route_binding"] != before["route_binding"]:
        raise ValueError("The source review must not alter the modeled product or attachment")
    molecule = Chem.MolFromSmiles(row["route_binding"]["precursor_smiles"])
    nitrogens = [a for a in molecule.GetAtoms() if a.GetAtomicNum() == 7]
    if len(nitrogens) != 2 or sum(a.GetTotalNumHs() for a in nitrogens) != 1 or rdMolDescriptors.CalcMolFormula(molecule) != "C12H8N2O4":
        raise ValueError("Precursor is not the reviewed monoethyl/NH PMDI graph")
    for name in ("attachment_selectivity_reviewed", "conditions_defined", "balanced_route_reviewed", "redox_core_retention_reviewed"):
        if row[name] != before[name]:
            raise ValueError("A precursor source cannot silently change other route assessments")
    script = ROOT / "scripts/pipeline/screen_grafting_routes.py"
    spec = importlib.util.spec_from_file_location("precursor_route_check", script)
    gate = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gate)
    feed = json.loads((ROOT / "results/candidates/current.json").read_text())
    candidate = next(r for r in feed["original39"] + feed["additional_candidates"] if r["id"] == "design_pmdi_ethyl")
    assessment = json.loads((ROOT / "results/grafting/selection_assessment.json").read_text())
    decision = gate.screen(candidate, evidence["candidates"], assessment)
    if decision["route_status"] != "hold" or decision["proceed_to_full_redox"] or "exact_precursor_access_reviewed" in decision["missing_route_evidence"]:
        raise ValueError("Access review was lost or incorrectly promoted to full route admission")
    sources = [Path(__file__).resolve(), args.pdf.resolve(), evidence_path, before_path, script]
    sources += [args.pdf.parent / f"EP0010440B1_page{page}.png" for page in (3, 10)]
    report = dict(
        status="exact_precursor_source_identity_reviewed_route_still_hold", checked=time.time(),
        candidate_id=candidate["id"], precursor_smiles=Chem.MolToSmiles(molecule),
        precursor_inchikey=Chem.MolToInchiKey(molecule), formula="C12H8N2O4",
        source_id="EP0010440B1", source_location="Example 19, printed page 10; scaffold scheme, page 3",
        identity_basis="Reported mono-N-ethyl product and sequence interpreted against the scheme with ethyl/H nitrogen substitution",
        rendered_pages_inspected=[3, 10], independent_synthesis_replication=False,
        current_commercial_availability_verified=False, exact_polymer_grafting_demonstrated=False,
        prospective_route_decision=decision, source_sha256={str(p): sha(p) for p in sources},
    )
    with (args.out / "summary.json").open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    (args.out / "producer_source.py").write_bytes(Path(__file__).read_bytes())
    print(json.dumps({key: report[key] for key in ("status", "precursor_inchikey", "formula")}))


if __name__ == "__main__":
    main()
