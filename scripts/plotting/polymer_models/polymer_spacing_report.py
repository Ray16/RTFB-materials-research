#!/usr/bin/env python
"""Render config-defined occupancy cases and field-qualified penalty results."""

import argparse
import hashlib
import html
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdCoordGen
from rdkit.Chem.Draw import rdMolDraw2D

from redox.core.io import atomic_write_json


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pattern(states):
    symbols = {"styrene": "S", "benzyl_chloride": "Cl", "grafted": "A"}
    return " — ".join(symbols[state] for state in states)


def render(campaign):
    current_path = campaign / "current.json"
    current = json.loads(current_path.read_text())
    cases = current["cases"]
    drawings = campaign / "structures"
    drawings.mkdir(exist_ok=True)
    cards = []
    for case in cases:
        molecule = Chem.MolFromSmiles(case["mapped_smiles"])
        backbone = [
            atom.GetIdx()
            for atom in molecule.GetAtoms()
            if 1000 < atom.GetAtomMapNum() <= 1000 + 2 * case["repeat_units"]
        ]
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
        # CoordGen resolves crowded pendant rings in larger exact graphs.
        # These 2D coordinates never enter geometry or energy calculations.
        rdCoordGen.AddCoords(molecule)
        bonds = [
            bond.GetIdx()
            for bond in molecule.GetBonds()
            if bond.GetBeginAtomIdx() in backbone and bond.GetEndAtomIdx() in backbone
        ]
        assert len(backbone) == 2 * case["repeat_units"]
        assert len(bonds) == len(backbone) - 1
        drawer = rdMolDraw2D.MolDraw2DSVG(900, 520)
        drawer.drawOptions().bondLineWidth = 2
        drawer.DrawMolecule(
            molecule,
            highlightAtoms=backbone,
            highlightBonds=bonds,
            highlightAtomColors={index: (0.42, 0.30, 0.58) for index in backbone},
            highlightBondColors={index: (0.42, 0.30, 0.58) for index in bonds},
        )
        drawer.FinishDrawing()
        filename = case["key"] + ".svg"
        (drawings / filename).write_text(drawer.GetDrawingText())
        cards.append(
            "<article><h3>"
            + html.escape(case["key"])
            + "</h3><p><b>"
            + html.escape(pattern(case["site_states"]))
            + "</b></p>"
            + f'<img src="structures/{filename}" alt="Exact {case["key"]} graph">'
            + f"<p>{case['atoms']} atoms including H; open backbone with "
            + f"{case['backbone_carbon_count']} carbons and "
            + f"{case['backbone_CC_bonds']} bonds.</p>"
            + "<p>Geometry: "
            + html.escape(case["geometry_status"])
            + " · G: "
            + html.escape(case["G_status"])
            + "</p></article>"
        )

    def number(value):
        return "pending" if value is None else f"{value:+.2f}"

    rows = []
    conditional_rows = []
    for cycle in current["cycles"]:
        for model, values in cycle["conditional_grafting"].items():
            fields = (
                "graft_i_on_empty_kJ_mol",
                "graft_j_on_empty_kJ_mol",
                "graft_j_given_i_kJ_mol",
                "graft_i_given_j_kJ_mol",
            )
            cells = [html.escape(cycle["name"]), html.escape(model)]
            cells += [number(values[field] if values else None) for field in fields]
            conditional_rows.append(
                "<tr>" + "".join("<td>" + cell + "</td>" for cell in cells) + "</tr>"
            )
        rows.append(
            "<tr><td>"
            + html.escape(cycle["name"])
            + "</td><td>"
            + str(cycle["sites"])
            + "</td><td>"
            + str(cycle["intervening_repeat_units"])
            + "</td><td>"
            + number(cycle["electronic_neighbor_penalty_deltaDeltaE_kJ_mol"])
            + "</td><td>"
            + number(cycle["free_energy_neighbor_penalty_deltaDeltaG_kJ_mol"])
            + "</td></tr>"
        )
    status = current["status"]
    page = """<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Polymer spacing and conditional grafting</title><style>
body{font:17px/1.55 system-ui,sans-serif;max-width:1350px;margin:35px auto;padding:0 20px;color:#15304a}
a{color:#0072b2}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:20px}
article{border:1px solid #dbe2e8;border-radius:8px;padding:14px;overflow-wrap:anywhere}img{width:100%}
.note{padding:15px;background:#fff4df}table{border-collapse:collapse;width:100%;font-size:15px}td,th{padding:10px;border-bottom:1px solid #ddd;text-align:left}.table{overflow-x:auto}</style>
<h1>Polymer spacing and conditional grafting</h1>
<p><a href="../../README.md">Polymer-model guide</a> · <a href="current.json">Current fields and QC</a> · <a href="plan.json">Frozen calculation plan</a></p>
<p class="note"><b>Conditional modeling scenarios, not a measured polymer loading.</b>
S = plain styrene padding; Cl = a reactive benzyl-chloride site; A = a grafted motif.
Each symbol denotes one repeat unit. Purple bonds show the actual open carbon backbone.
The same tagged backbone is preserved across matched cases. A drawing is a connectivity diagram, not an optimized conformation.</p>"""
    page += (
        f"<p><b>{status['completed_DFT']}/{status['DFT_tasks']} qualified DFT structures; "
        f"{status['completed_native_G']}/{status['graph_cases']} qualified native G values.</b>"
        " Blank scientific quantities remain pending until the required records pass QC.</p>"
    )
    page += """<h2>Two dimensions</h2><p>Electronic site nonadditivity measures the sampled effect of occupying both sites, including backbone/conformational response. Conditional grafting ΔG adds the qualified precursor/base/byproduct terms under defined conditions.</p>
<p><b>Neighbor penalty:</b> ΔΔG = G(AA) + G(00) − G(A0) − G(0A) = ΔG(graft j | i occupied) − ΔG(graft j | empty).</p>
<p>A positive penalty makes the next graft less favorable. Its full grafting ΔG can still be negative. Electronic ΔΔE is never relabeled as ΔΔG.</p>
<div class="table"><table><thead><tr><th>Scenario</th><th>Occupied sites</th><th>Intervening units</th><th>ΔΔE (kJ/mol)</th><th>ΔΔG (kJ/mol)</th></tr></thead><tbody>"""
    page += (
        "".join(rows)
        + "</tbody></table></div><h2>Conditional grafting ΔG</h2>"
        + "<p>Sites i and j follow the scenario order above. All values are kJ/mol; pending values require qualified polymer and reagent thermochemistry.</p>"
        + '<div class="table"><table><thead><tr><th>Scenario</th><th>Byproduct model</th><th>i on empty</th><th>j on empty</th><th>j with i occupied</th><th>i with j occupied</th></tr></thead><tbody>'
        + "".join(conditional_rows)
        + '</tbody></table></div><h2>Exact modeled cases</h2><div class="cards">'
        + "".join(cards)
        + "</div>"
    )
    page += "<p>Native G uses the selected sampled conformer, not an equilibrated conformational ensemble. Model favorability does not establish attachment rates, yields, graft spacing, swelling, electronic coupling or whole-bead performance.</p></html>"
    (campaign / "index.html").write_text(page)
    atomic_write_json(
        campaign / "figure_manifest.json",
        dict(
            source_sha256=digest(current_path),
            producer_sha256=digest(__file__),
            outputs={
                str(path): digest(path)
                for path in [campaign / "index.html", *drawings.glob("*.svg")]
            },
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    args = parser.parse_args()
    render(args.campaign.resolve())
