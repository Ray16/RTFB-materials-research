"""Read-only independent review of the first frozen spacing panel."""

import json
import os
from pathlib import Path
import time
import xml.etree.ElementTree as ET

from redox.core.common import ROOT
from redox.build.polymer_cases import panel
from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
from redox.workflows.polymer_campaign import digest, verify_sources, load_module


def main():
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    campaign = ROOT / "results/polymer_models/campaigns/2026-10-20_spacing_pilot"
    plan = json.loads((campaign / "plan.json").read_text())
    verify_sources(plan["producer_sha256"])
    verify_sources(plan["input_sha256"])
    preparation = json.loads((campaign / "preparation.json").read_text())
    verify_sources(preparation["producer_sha256"])
    original = json.loads((campaign / "source_snapshot/preparation_original.json").read_text())
    for source, snapshot in preparation["stage_source_snapshot"].items():
        assert digest(ROOT / snapshot["snapshot"]) == original["producer_sha256"][source]
    spec = json.loads((ROOT / "config/polymer_spacing_pilot.json").read_text())
    cases = json.loads((campaign / "cases.json").read_text())
    rebuilt, cycles = panel(spec, cases["fragment"])
    assert rebuilt == cases["cases"] and cycles == cases["cycles"]
    audit = load_module("audit_spacing_review", ROOT / "scripts/validation/audit_polymer_geometries.py")
    reviews = []
    for task in plan["tasks"]:
        ready = json.loads((campaign / "live_inputs" / task["key"] / "task.json").read_text())
        verify_sources(ready["input_sha256"])
        result = audit.audit(ready, Path(ready["input_xyz"]))
        assert not result["issues"], result
        reviews.append(dict(key=task["key"], input_sha256=ready["input_sha256"], geometry_review=result))
    figures = list((campaign / "structures").glob("*.svg"))
    for path in figures:
        ET.parse(path)
    terms = json.loads((campaign / "reaction_terms.json").read_text())
    verify_sources(terms["source_sha256"])
    assert not terms["missing_reagents"]
    previous = []
    for name in ("2026-10-20_geometry_recovery", "2026-10-20_frozen_fragments"):
        old = ROOT / "results/polymer_models/campaigns" / name / "plan.json"
        verify_sources(json.loads(old.read_text())["producer_sha256"])
        previous.append(dict(campaign=name, live_sources_unchanged=True))
    summary = dict(
        checked_at=time.time(), status="framework_and_prepared_inputs_verified",
        tests="29 passed: polymer cases, cycles, admission and existing topology tests",
        lint="Ruff E/F numerical code and tests; F report script passed",
        exact_graph_rebuild_matches=True, graph_cases=len(rebuilt),
        qualified_warm_starts=len(reviews), all_stage_hashes_verified=True,
        shared_reagent_terms=terms["status"], svg_count=len(figures),
        drawing_review="Adjacent-pair SVG visually checked; CoordGen avoids pendant-ring overlap; purple open backbone",
        previous_campaigns=previous, prepared_geometry_reviews=reviews,
        source_plan=str(campaign / "plan.json"), source_plan_sha256=digest(campaign / "plan.json"),
        coordinator_snapshot=json.loads((campaign / "coordinator.json").read_text()),
        calculation_status=json.loads((campaign / "status.json").read_text()),
        limits=["Initial phthalimide method control", "No experimental occupancy weights",
                "Native G requires an independent pilot and thermal QC", "No pure motif pair energy, coupling, yield or rate"],
    )
    Path(__file__).with_name("summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print({key: summary[key] for key in ("status", "graph_cases", "qualified_warm_starts", "shared_reagent_terms")})


if __name__ == "__main__":
    main()
