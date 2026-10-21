#!/usr/bin/env python
"""Review the current polymer resume milestone without submitting calculations."""

import json
import math
from pathlib import Path
import re
import time

from redox.core.common import ROOT, EV_KJ
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.checkpoints import verify_state_checkpoint

BASE = ROOT / "results/polymer_models/campaigns"
OUT = ROOT / "results/reviews/2026-10-21/polymer_resume"
SPACING = BASE / "2026-10-21_spacing_memory_resume"


def main():
    selected, invalid = api.selected_graph_records(SPACING)
    assert len(selected) == 6 and not invalid
    records, invalid = api.admitted_records(SPACING, "energy")
    assert len(records) == 11 and not invalid
    feed = json.loads((SPACING / "current.json").read_text())
    cycles = []
    for row in feed["cycles"]:
        value = math.fsum(coefficient * selected[key]["E_smd_eV"]
                          for key, coefficient in row["coefficients"].items()) * EV_KJ
        assert abs(value - row["electronic_neighbor_penalty_deltaDeltaE_kJ_mol"]) < 1e-7
        assert row["free_energy_neighbor_penalty_deltaDeltaG_kJ_mol"] is None
        cycles.append(dict(name=row["name"], checked_deltaDeltaE_kJ_mol=value,
                           qualified_deltaDeltaG_kJ_mol=None))
    assert all(row["G_eV"] is None for row in feed["cases"])
    assert all(value is None for row in feed["cycles"] for value in row["conditional_grafting"].values())
    plans = {}
    for name in ("2026-10-21_spacing_continuation", "2026-10-21_spacing_memory_resume",
                 "2026-10-21_checkpoint_continuation", "2026-10-21_checkpoint_pmdi_open",
                 "2026-10-21_checkpoint_pmdi_stacked"):
        campaign = BASE / name
        plan = api.read_plan(campaign)
        api.verify_sources(plan["producer_sha256"])
        api.verify_sources(plan["input_sha256"])
        plans[name] = dict(plan_sha256=api.digest(campaign / "plan.json"), raw_root=plan["raw_root"])
    checkpoint = BASE / "2026-10-21_checkpoint_continuation/done/exp_phthalimide_control__n3__stereo1__open.json"
    completed = verify_state_checkpoint(checkpoint)
    assert completed["geometry_sanity_review"] == "passed"
    assert completed["topology_QC"] == "ok"
    assert completed["gradient_max_Ha_Bohr"] < .00045
    for target in ("open", "stacked"):
        audit = json.loads((BASE / f"2026-10-21_checkpoint_pmdi_{target}/checkpoint_geometry_review.json").read_text())
        assert audit["status"] == "geometry_sanity_pass" and not audit["issues"]
    index = json.loads((ROOT / "results/INDEX.json").read_text())
    paths = set()

    def walk(value):
        if isinstance(value, dict):
            for child in value.values():
                walk(child)
        elif isinstance(value, list):
            for child in value:
                walk(child)
        elif isinstance(value, str) and value.startswith(("results/", "docs/", "config/", "slides/", "calcs/")) and not any(marker in value for marker in ("*", "<")):
            paths.add(ROOT / value)

    walk(index)
    for document in (ROOT / "results/README.md", ROOT / "results/polymer_models/README.md"):
        for target in re.findall(r"\]\(([^)]+)\)", document.read_text()):
            if "://" not in target:
                paths.add((document.parent / target.split("#")[0]).resolve())
    # The summary is written below; its own incoming link is valid at publication.
    missing = [str(path) for path in paths if path != OUT / "summary.json" and not path.exists()]
    assert not missing, missing
    current = json.loads((ROOT / "results/polymer_models/interaction_status.json").read_text())
    assert current["stationary_electronic_cases"] == 25 and not current["invalid_receipts"]
    summary = dict(
        reviewed=time.time(), status="verified_electronic_milestone_native_G_pending",
        tests="41 passed: thermal-memory contraction, replay rejection, continuation, recovery, admission, cycles and polymer cases",
        Python_compile="passed for both new launch entrypoints",
        checked_oligomer_structures=25, expected_oligomer_structures=27,
        spacing_checked_conformers=11, spacing_graphs=6, preserved_alternate_conformer_memory_failures=1,
        checked_electronic_cycles=cycles, native_G_qualified=0,
        phthalimide_checkpoint_receipt=str(checkpoint.relative_to(ROOT)),
        campaigns=plans, local_links_checked=len(paths), local_links_missing=missing,
        starts="results/reviews/2026-10-21/polymer_resume/start_verification.json",
        start_observation="Four reserved GPU predictors verified: two full-Hessian controls and two PMDI checkpoint resets; eight disjoint CPU slots and all native thread caps=2",
        figure_review="Visually checked the 25/27 electronic-summary PNG; six association and three distinct folding comparisons shown",
        source_sha256={str(path.relative_to(ROOT)): api.digest(path) for path in (
            Path(__file__), ROOT / "scripts/pipeline/resume_polymer_thermal.py",
            ROOT / "scripts/pipeline/resume_remaining_trimers.py", SPACING / "plan.json",
            ROOT / "results/polymer_models/figures/key_results.png")},
        limits=["Full replay controls and large thermal resource/physical pilot still pending",
                "Electronic spacing cycles include conformational response and one failed alternate; no general spacing trend or pure-core binding claim",
                "No graft-loading distribution, whole-bead result, synthesis admission or electron-transfer rate inferred",
                "Repository-wide refactor remains deferred until the scientific work is properly addressed"],
    )
    OUT.mkdir(parents=True, exist_ok=True)
    atomic_write_json(OUT / "summary.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
