#!/usr/bin/env python
"""Compare a defined acid-capture term with exact-identity MeCN pKa anchors.

The hybrid grafting values are term sensitivities with a fixed DFT reference.
They are not experimentally calibrated full grafting reaction free energies.
"""

import argparse
from html.parser import HTMLParser
import json
import math
from pathlib import Path
import time

from rdkit import Chem, rdBase

from redox.core.io import atomic_write_json
from redox.workflows.chunked import file_hash


class BasicityTable(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows, self.row, self.cell, self.urls = [], None, None, None

    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self.row = []
        elif tag in ("td", "th") and self.row is not None:
            self.cell, self.urls = [], []
        elif tag == "a" and self.cell is not None:
            href = dict(attrs).get("href")
            if href:
                self.urls.append(href)

    def handle_data(self, value):
        if self.cell is not None:
            self.cell.append(value)

    def handle_endtag(self, tag):
        if tag in ("td", "th") and self.cell is not None:
            self.row.append(dict(text="".join(self.cell).strip(), urls=self.urls))
            self.cell, self.urls = None, None
        elif tag == "tr" and self.row is not None:
            if len(self.row) == 6:
                self.rows.append(self.row)
            self.row = None


def identity(smiles):
    with rdBase.BlockLogs():
        molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
        return None
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(molecule, isomericSmiles=True)


def exact_base(table, smiles):
    expected = identity(smiles)
    if expected is None:
        raise ValueError("Invalid target base graph")
    matches = []
    for row in table.rows:
        if identity(row[5]["text"]) != expected:
            continue
        try:
            value = float(row[3]["text"])
        except ValueError:
            continue
        matches.append(dict(CAS=row[0]["text"], name=row[1]["text"], pKaH_MeCN=value,
                            smiles=row[5]["text"], source_urls=row[3]["urls"], exact_graph=expected))
    if len(matches) > 1:
        raise ValueError("Multiple exact-base records require source review")
    return matches[0] if matches else None


def capture_deltaG(pKa_acid, pKa_conjugate_base, temperature_K=298.15):
    if not all(math.isfinite(x) for x in (pKa_acid, pKa_conjugate_base, temperature_K)) or temperature_K <= 0:
        raise ValueError("Finite pKa values and positive temperature required")
    return -.00831446261815324 * temperature_K * math.log(10) * (pKa_conjugate_base - pKa_acid)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basicity-html", type=Path, required=True)
    parser.add_argument("--acid-anchor", type=Path, required=True)
    parser.add_argument("--base-results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    table = BasicityTable()
    table.feed(args.basicity_html.read_text())
    acid = json.loads(args.acid_anchor.read_text())
    if acid["solvent"] != "acetonitrile" or identity(acid["smiles"]) != identity("Cl"):
        raise ValueError("The current source term requires the exact HCl MeCN anchor")
    source = json.loads(args.base_results.read_text())
    rows = []
    for key, base in source["bases"].items():
        match = exact_base(table, base["neutral"])
        value = source["corrections"][key]["free_ions_kJ_mol"]
        row = dict(base=key, model_smiles=base["neutral"], DFT_free_ion_capture_kJ_mol=value,
                   reference=match, reference_status="exact_base_in_author_collection" if match else
                   "no_exact_base_anchor_in_checked_collection", external_capture_kJ_mol=None,
                   DFT_minus_external_kJ_mol=None, hybrid_term_sensitivity=[])
        if match is not None:
            reference = capture_deltaG(acid["pKa_MeCN"], match["pKaH_MeCN"])
            row.update(external_capture_kJ_mol=reference, DFT_minus_external_kJ_mol=value-reference)
            for candidate in source["candidates"]:
                original = candidate["base_models"][key]["free_ions_kJ_mol"]
                baseline = candidate["reference_deltaG_kJ_mol"]
                if original is None or baseline is None:
                    continue
                hybrid = baseline + reference
                row["hybrid_term_sensitivity"].append(dict(id=candidate["id"], fixed_DFT_reference_kJ_mol=baseline,
                    all_DFT_base_assisted_kJ_mol=original, fixed_reference_with_external_capture_kJ_mol=hybrid,
                    displayed_sign_retained=bool((original < 0) == (hybrid < 0))))
        rows.append(row)
    result = dict(status="external_term_comparison_scope_limited", evaluated_at=time.time(),
                  temperature_K=298.15, reaction="B + HCl -> BH+ + Cl-", delta_n_solutes=0,
                  acid_anchor=acid, rows=rows, full_grafting_accuracy_validated=False,
                  calibration_fitted=False, contact_pair_comparison="No measured association anchor; not inferred from pKa",
                  limits=["Author-tabulated solvent pKa anchors; source temperature/activity details still require review",
                          "Exact DIPEA may not be replaced by diisopropylamine or triethylamine",
                          "HCl cancels in the complete all-DFT cycle; term error cannot be assigned to an individual species",
                          "Replacing one term leaves DFT-reference error and shared-error cancellation unresolved",
                          "Hybrid values are a sensitivity calculation, not measured or calibrated full grafting energies"],
                  source_sha256={str(p.resolve()): file_hash(p) for p in [args.basicity_html, args.acid_anchor,
                                  args.base_results, Path(__file__)]})
    args.out.mkdir(parents=True, exist_ok=False)
    atomic_write_json(args.out / "comparison.json", result)
    print(json.dumps([dict(base=r["base"], reference_status=r["reference_status"],
                          DFT_minus_external_kJ_mol=r["DFT_minus_external_kJ_mol"]) for r in rows]))


if __name__ == "__main__":
    main()
