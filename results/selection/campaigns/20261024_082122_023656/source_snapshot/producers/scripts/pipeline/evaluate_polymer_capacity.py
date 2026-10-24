#!/usr/bin/env python
"""Composition scenarios and necessary capacity bounds for an exact candidate.

No experimental loading, accessibility or positional probability is inferred.
One- and two-electron budgets are assumptions, not polymer redox measurements.
"""
import argparse
from itertools import combinations
import json
import os
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, rdMolDescriptors

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.properties.polymer_capacity import fragment_capacity_bound, FARADAY_C_MOL
from redox.build.polymer_cases import build_case
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    args = parser.parse_args()
    base = args.campaign.resolve()
    output = base / "capacity_scenarios.json"
    if output.exists():
        raise FileExistsError("Preserve the existing scenario publication; use a new output campaign for changes")
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    candidate_path = base / "candidate_snapshot.json"
    candidate = json.loads(candidate_path.read_text())
    plan = json.loads((base / "plan.json").read_text())
    fragment = candidate["frag"]
    n = 5
    grid = []
    for handles in range(n + 1):
        for grafts in range(handles + 1):
            row = dict(repeat_units=n, reactive_handles=handles, grafts=grafts)
            for electrons in (1, 2):
                row[str(electrons) + "e"] = fragment_capacity_bound(fragment, n, handles, grafts, electrons)
            grid.append(row)
    exact_graphs = []
    for grafts in range(4):
        for occupied in combinations([2, 3, 4], grafts):
            case = build_case(fragment, n, [2, 3, 4], list(occupied), ["CW"] * 4)
            molecule = Chem.MolFromSmiles(case["smiles"])
            mass = Descriptors.MolWt(molecule)
            bound = fragment_capacity_bound(fragment, n, 3, grafts, 2)
            if abs(bound["mass_g_mol"] - mass) > 1e-8:
                raise ValueError("Analytic mass differs from the exact occupied polymer graph")
            exact_graphs.append(dict(occupied_sites=list(occupied), smiles=case["smiles"],
                formula=rdMolDescriptors.CalcMolFormula(molecule), polymer_graph_key=Chem.MolToInchiKey(molecule),
                mass_g_mol=mass, one_e_upper_bound_mAh_g=bound["capacity_upper_bound_mAh_g"] / 2,
                two_e_upper_bound_mAh_g=bound["capacity_upper_bound_mAh_g"],
                geometry_search_prepared=list(occupied) in ([2], [3], [2, 3]),
                extra_quantum_jobs_submitted_by_this_analysis=0))
    capped_mass = Descriptors.MolWt(Chem.MolFromSmiles(candidate["exact_grafted_smiles"]))
    source = [Path(__file__).resolve(), ROOT / "src/redox/properties/polymer_capacity.py", ROOT / "src/redox/build/polymer_cases.py",
        candidate_path, base / "plan.json", base / "input_validation.json"]
    report = dict(status="conditional_stoichiometric_upper_bounds_checked_against_exact_graphs",
        candidate_id=candidate["id"], exact_capped_product_key=plan["exact_capped_product_key"],
        enumerated_aliases=plan["enumerated_aliases"], reference_capped_mass_g_mol=capped_mass,
        reference_capped_two_e_nominal_mAh_g=2 * FARADAY_C_MOL / (3.6 * capped_mass),
        composition_grid=grid, exact_graphs_at_handle_density_0_6=exact_graphs,
        selected_application_target_mAh_g=None, capacity_lower_bound_mAh_g=0,
        accessible_capacity_established=False, experimental_loading_inferred=False,
        acceptance_rule="A target above the bound rejects that declared composition/electron-budget scenario. A target below it does not admit the material.",
        mass_convention="Dry neutral H-capped linear fragment; residual chloride retained; unlisted crosslinker/electrolyte mass excluded",
        electron_budget_scope="One or two electrons per graft assumed; accessible polymer redox path not measured or computed here",
        rdkit_version=rdBase.rdkitVersion, source_sha256={str(p): api.digest(p) for p in source})
    atomic_write_json(output, report)
    print(json.dumps(dict(candidate=candidate["id"], conditional_graphs=len(exact_graphs),
        at_three_handles=[dict(grafts=row["grafts"], upper_mAh_g=row["2e"]["capacity_upper_bound_mAh_g"]) for row in grid if row["reactive_handles"] == 3])))


if __name__ == "__main__":
    main()
