#!/usr/bin/env python
"""Draw exact PMDI occupancy graphs used for the candidate feasibility searches."""
import argparse
import json
import os
from pathlib import Path
from xml.etree import ElementTree as ET

import cairosvg
from rdkit import Chem
from rdkit.Chem import rdCoordGen
from rdkit.Chem.Draw import rdMolDraw2D

from redox.core.common import ROOT
from redox.build.polymer_contacts import pi_core_groups
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=ROOT / "results/figures/polymer_models/pmdi_occupancy_construction")
    args = parser.parse_args()
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    plan_path = args.campaign / "plan.json"
    plan = json.loads(plan_path.read_text())
    validation_path = args.campaign / "input_validation.json"
    validation = json.loads(validation_path.read_text())
    if validation["status"] != "exact_candidate_occupancy_inputs_checked":
        raise ValueError("Prepared structures require exact graph/route input review")
    api.verify_sources(validation["source_sha256"])
    candidate = json.loads((args.campaign / "candidate_snapshot.json").read_text())
    paths = [Path(plan["common_empty_reference"]["case_path"])] + [Path(row["case"]) for row in plan["tasks"]]
    rows = {tuple(json.loads(path.read_text())["occupied_sites"]): (json.loads(path.read_text()), path) for path in paths}
    svg = ET.Element("svg", xmlns="http://www.w3.org/2000/svg", width="7.2in", height="6.5in", viewBox="0 0 1500 1350")
    ET.SubElement(svg, "rect", width="1500", height="1350", fill="white")
    shown = []
    for panel, (occupied, label) in enumerate((((), "a  Shared empty backbone"), ((2,), "b  Graft at site 2"), ((2, 3), "c  Grafts at sites 2 and 3"))):
        case, path = rows[occupied]
        mapped = Chem.MolFromSmiles(case["mapped_smiles"])
        backbone = {atom.GetAtomMapNum() - 1000: atom.GetIdx() for atom in mapped.GetAtoms() if 1000 < atom.GetAtomMapNum() <= 1010}
        if len(backbone) != 10 or any(not mapped.GetBondBetweenAtoms(backbone[i], backbone[i + 1]) for i in range(1, 10)):
            raise ValueError("Drawing lacks the exact open ten-carbon backbone")
        plain = Chem.Mol(mapped)
        for atom in plain.GetAtoms():
            atom.SetAtomMapNum(0)
        if Chem.MolToSmiles(plain) != Chem.MolToSmiles(Chem.MolFromSmiles(case["smiles"])):
            raise ValueError("Drawing differs from the prepared exact molecular graph")
        for site in range(1, 6):
            plain.GetAtomWithIdx(backbone[2 * site]).SetProp("atomNote", str(site))
        groups = pi_core_groups(plain, candidate["frag"])
        core = set(i for group in groups for i in group)
        for i in list(core):
            for bond in plain.GetAtomWithIdx(i).GetBonds():
                other = bond.GetOtherAtomIdx(i)
                if bond.GetBondType() == Chem.BondType.DOUBLE and plain.GetAtomWithIdx(other).GetAtomicNum() == 8:
                    core.add(other)
        if len(groups) != len(occupied):
            raise ValueError("Drawn motif count differs from occupancy")
        rdCoordGen.AddCoords(plain)
        if Chem.MolToSmiles(plain) != Chem.MolToSmiles(Chem.MolFromSmiles(case["smiles"])):
            raise ValueError("Depiction changed the declared stereochemical graph")
        drawer = rdMolDraw2D.MolDraw2DSVG(1500, 390)
        options = drawer.drawOptions()
        options.bondLineWidth = 1.8
        options.minFontSize = 24
        options.maxFontSize = 28
        options.annotationFontScale = 0.8
        options.padding = 0.08
        colors = {i: (0.52, 0.75, 0.88) for i in backbone.values()}
        colors.update({i: (0.82, 0.65, 0.82) for i in core})
        bonds = [plain.GetBondBetweenAtoms(backbone[i], backbone[i + 1]).GetIdx() for i in range(1, 10)]
        drawer.DrawMolecule(plain, highlightAtoms=sorted(colors), highlightBonds=bonds,
            highlightAtomColors=colors, highlightBondColors={i: (0.52, 0.75, 0.88) for i in bonds})
        drawer.FinishDrawing()
        child = ET.fromstring(drawer.GetDrawingText())
        child.set("x", "0")
        child.set("y", str(panel * 450 + 50))
        svg.append(child)
        text = ET.SubElement(svg, "text", x="15", y=str(panel * 450 + 34), style="font-family:sans-serif;font-size:30px;fill:black")
        text.text = label
        shown.append(dict(occupied_sites=list(occupied), smiles=case["smiles"], source=str(path), backbone_carbons=10, backbone_CC_bonds=9, chain_ends=2))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.with_suffix(".svg").write_text(ET.tostring(svg, encoding="unicode"))
    cairosvg.svg2pdf(url=str(args.out.with_suffix(".svg")), write_to=str(args.out.with_suffix(".pdf")))
    cairosvg.svg2png(url=str(args.out.with_suffix(".svg")), write_to=str(args.out.with_suffix(".png")), output_width=2880)
    metadata = dict(status="exact_prepared_candidate_graphs_drawn", candidate_id=plan["candidate_id"], enumerated_aliases=plan["enumerated_aliases"], cases=shown,
        source_sha256={str(path): api.digest(path) for path in [Path(__file__).resolve(), plan_path, validation_path, *paths]},
        output_sha256={str(args.out.with_suffix("." + suffix)): api.digest(args.out.with_suffix("." + suffix)) for suffix in ("svg", "pdf", "png")},
        caption="Exact five-repeat PMDI feasibility models. Blue marks the open backbone and purple the grafted PMDI cores; numbers identify repeat units. Sites 1 and 5 are styrene padding and unoccupied sites 2–4 retain benzyl chloride. A separate site-3 single-graft search is also prepared. These are graph drawings, not DFT minima, representative 3D conformations, measured loading or validated synthesis.")
    args.out.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
