#!/usr/bin/env python
"""Original chemical TOC artwork from the exact PMDI candidate identities.

State arrows and neighbor cartoons describe assessment dimensions. No curves,
numerical results, passed gates or experimentally established route are implied.
"""
import argparse
import io
import json
import os
from pathlib import Path
import time
from xml.etree import ElementTree as ET

import cairosvg
from PIL import Image
from pypdf import PdfReader
from rdkit import Chem
from rdkit.Chem import rdCoordGen
from rdkit.Chem.Draw import rdMolDraw2D

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", NS)


def element(parent, tag, **attributes):
    return ET.SubElement(parent, "{" + NS + "}" + tag,
                         {k.replace("_", "-"): str(v) for k, v in attributes.items()})


def label(parent, x, y, value, size=40, **kwargs):
    item = element(parent, "text", x=x, y=y, font_family="DejaVu Sans",
                   font_size=size, text_anchor="middle", fill="#222222", **kwargs)
    item.text = value
    return item


def molecule(parent, smiles, x, width, height=305):
    mol = Chem.MolFromSmiles(smiles)
    before = Chem.MolToSmiles(mol)
    rdCoordGen.AddCoords(mol)
    if Chem.MolToSmiles(mol) != before:
        raise ValueError("Chemical depiction changed the exact graph")
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opt = drawer.drawOptions()
    opt.minFontSize, opt.maxFontSize, opt.bondLineWidth = 38, 42, 2.4
    opt.padding = .06
    atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in (7, 8)]
    drawer.DrawMolecule(mol, highlightAtoms=atoms,
                        highlightAtomColors={i: (.86, .73, .86) for i in atoms})
    drawer.FinishDrawing()
    child = ET.fromstring(drawer.GetDrawingText())
    child.set("x", str(x))
    child.set("y", "82")
    parent.append(child)
    return dict(smiles=before, inchikey=Chem.MolToInchiKey(mol))


def create(source, stem):
    candidate = json.loads(source.read_text())
    if candidate["id"] != "design_pmdi_ethyl" or candidate["exact_grafted_inchikey"] != "YPRQKNJFGRGANA-UHFFFAOYSA-N":
        raise ValueError("TOC artwork requires the reviewed exact candidate")
    states = candidate["states"]
    states = json.loads(states) if isinstance(states, str) else states
    if [(row[1], row[2]) for row in states] != [(0, 1), (-1, 2), (-2, 1)]:
        raise ValueError("Illustrated charge ladder differs from the recorded candidate states")
    policy_path = ROOT / "config/candidate_feasibility.json"
    policy = json.loads(policy_path.read_text())
    if len(policy["required_axes"]) != 10 or "All required axes must pass" not in policy["rule"]:
        raise ValueError("All-axis policy differs from the illustrated decision")
    svg = ET.Element("{" + NS + "}svg", width="3.25in", height="1.75in", viewBox="0 0 1300 700")
    element(svg, "rect", x=0, y=0, width=1300, height=700, fill="white")
    defs = element(svg, "defs")
    marker = element(defs, "marker", id="arrow", viewBox="0 0 10 10", refX=9, refY=5,
                     markerWidth=8, markerHeight=8, orient="auto-start-reverse")
    element(marker, "path", d="M 0 0 L 10 5 L 0 10 z", fill="#333333")
    label(svg, 230, 52, "Exact precursor", size=42)
    label(svg, 960, 52, "Graft-site model", size=42)
    precursor = molecule(svg, candidate["precursor_smiles"], 10, 445)
    product = molecule(svg, candidate["exact_grafted_smiles"], 630, 660)
    element(svg, "line", x1=470, y1=225, x2=610, y2=225, stroke="#333333",
            stroke_width=4, marker_end="url(#arrow)")
    label(svg, 540, 151, "Identity", size=36)
    label(svg, 540, 192, "mapping", size=36)
    element(svg, "line", x1=25, y1=405, x2=1275, y2=405, stroke="#b0b0b0", stroke_width=2)
    for x, title in ((215, "Redox states"), (650, "Defined route"), (1085, "Neighbor effects")):
        label(svg, x, 454, title, size=40)
    for x, charge in ((95, "0"), (215, "−1"), (335, "−2")):
        element(svg, "circle", cx=x, cy=530, r=43, fill="#e4f3fb", stroke="#0072b2", stroke_width=4)
        label(svg, x, 545, charge, size=43)
    for a, b in ((143, 164), (263, 284)):
        element(svg, "line", x1=a, y1=530, x2=b, y2=530, stroke="#333333", stroke_width=3,
                marker_end="url(#arrow)")
    label(svg, 650, 525, "ΔG (route)", size=44)
    label(svg, 650, 578, "Base / MeCN", size=37)
    points = [(950 + 30 * i, 500 + (i % 2) * 24) for i in range(10)]
    element(svg, "polyline", points=" ".join(f"{x},{y}" for x, y in points),
            fill="none", stroke="#0072b2", stroke_width=7)
    for index in (3, 5):
        x, y = points[index]
        element(svg, "line", x1=x, y1=y, x2=x, y2=y + 35, stroke="#333333", stroke_width=3)
        element(svg, "circle", cx=x, cy=y + 50, r=22, fill="#cc79a7", stroke="#333333", stroke_width=2)
    label(svg, 1085, 636, "Spacing / coupling", size=36)
    label(svg, 650, 683, "All criteria must pass", size=46, font_weight="bold")
    stem.parent.mkdir(parents=True, exist_ok=True)
    if any(stem.with_suffix("." + ext).exists() for ext in ("svg", "pdf", "eps", "png", "tif", "json")):
        raise FileExistsError("Preserve the existing artwork; use a new output stem")
    data = ET.tostring(svg, encoding="utf-8", xml_declaration=True)
    stem.with_suffix(".svg").write_bytes(data)
    cairosvg.svg2pdf(bytestring=data, write_to=str(stem.with_suffix(".pdf")))
    cairosvg.svg2eps(bytestring=data, write_to=str(stem.with_suffix(".eps")))
    png = cairosvg.svg2png(bytestring=data, output_width=1950, output_height=1050)
    stem.with_suffix(".png").write_bytes(png)
    im = Image.open(io.BytesIO(png)).convert("RGB")
    im.save(stem.with_suffix(".tif"), dpi=(600, 600), compression="tiff_lzw")
    pdf = PdfReader(stem.with_suffix(".pdf"))
    box = pdf.pages[0].mediabox
    if abs(float(box.width) - 234) > 1e-5 or abs(float(box.height) - 126) > 1e-5:
        raise ValueError("TOC PDF dimensions differ from the target physical size")
    report = dict(status="chemical_TOC_artwork_prepared_requires_author_review", created=time.time(),
        width_in=3.25, height_in=1.75, raster_dpi=600, minimum_explicit_text_pt=6.48,
        precursor=precursor, capped_product=product, exact_candidate_id=candidate["id"],
        redox_state_scope="Recorded computed charge/multiplicity assignments; no atom-localized charge or resonance claim",
        neighbor_cartoon=dict(backbone_vertices=10, bonds=9, ends=2,
            scope="Assessment schematic: substituent chemistry, stereochemistry and spatial contacts are not represented"),
        scientific_claim="Three linked assessment dimensions; all ten required axes must pass. No quantitative results or successful synthesis implied.",
        source_sha256={str(p): api.digest(p) for p in (Path(__file__).resolve(), source, policy_path)},
        output_sha256={str(stem.with_suffix("." + ext)): api.digest(stem.with_suffix("." + ext))
                       for ext in ("svg", "pdf", "eps", "png", "tif")},
        publisher_specification="https://pubsapp.acs.org/paragonplus/submission/toc_abstract_graphics_guidelines.pdf",
        caption="Exact PMDI precursor and capped graft-site identities link molecular redox assessment to a defined grafting route and conditional neighbor analysis. The bottom-right open-chain cartoon represents spacing/coupling assessment; it is not a molecular geometry or loading measurement. All required criteria must pass before selection.")
    stem.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, default=ROOT / "results/validation/pmdi_occupancy_2026-10-23/candidate_snapshot.json")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    report = create(args.candidate.resolve(), args.out.resolve())
    print(json.dumps(dict(status=report["status"], size_in=[report["width_in"], report["height_in"]])))


if __name__ == "__main__":
    main()
