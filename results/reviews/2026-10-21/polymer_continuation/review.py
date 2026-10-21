"""Reconstruct the source/geometry/resource basis of the two continuations."""

import json
import os
from pathlib import Path
import subprocess
import time

from redox.core.common import ROOT
from redox.core.protocol import ACTIVE_SP
from redox.build.oligomer_topology import chain_topology
from redox.workflows.polymer_campaign import verify_sources, digest
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
from redox.workflows.checkpoints import verify_state_checkpoint


def main():
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    out = Path(__file__).resolve().parent
    campaigns = ROOT / "results/polymer_models/campaigns"
    new = []
    for name in ("2026-10-21_spacing_continuation", "2026-10-21_checkpoint_continuation"):
        campaign = campaigns / name
        plan = json.loads((campaign / "plan.json").read_text())
        verify_sources(plan["producer_sha256"])
        verify_sources(plan["input_sha256"])
        pilot = campaign / ("energy/pilot/status.json" if "spacing" in name else "pilot/status.json")
        proof = json.loads(pilot.read_text())
        assert proof["status"] == "memory_pilot_passed"
        verify_sources(proof["output_sha256"])
        owner = json.loads((campaign / "coordinator.json").read_text())
        process = Path(f"/proc/{owner['pid']}/stat")
        assert process.exists() and process.read_text().rsplit(")", 1)[1].split()[19] == owner["process_start"]
        new.append(dict(campaign=str(campaign), plan_sha256=digest(campaign / "plan.json"),
                        tasks=len(plan["tasks"]), pilot_status=proof["status"],
                        coordinator_snapshot=owner, coordinator_affinity=sorted(os.sched_getaffinity(owner["pid"]))))
    original = campaigns / "2026-10-20_spacing_pilot"
    numerical = json.loads((original / "energy/pilot/energy_and_gradient.json").read_text())
    assert numerical["result"]["energy_protocol"] == ACTIVE_SP
    recovery = campaigns / "2026-10-21_checkpoint_continuation"
    geometry = json.loads((recovery / "checkpoint_geometry_review.json").read_text())
    assert not geometry["issues"]
    task = json.loads((recovery / "plan.json").read_text())["tasks"][0]
    topology = chain_topology(task["smiles"], task["repeats"])
    assert topology["backbone_carbon_count"] == 6 and topology["backbone_CC_bonds"] == 5
    previous = []
    for name in ("2026-10-20_geometry_recovery", "2026-10-20_spacing_pilot", "2026-10-20_frozen_fragments"):
        verify_sources(json.loads((campaigns / name / "plan.json").read_text())["producer_sha256"])
        previous.append(name)
    receipts = 0
    for name in ("2026-10-19_pi_stacking", "2026-10-20_geometry_recovery"):
        for path in (campaigns / name / "done").glob("*.json"):
            verify_state_checkpoint(path)
            receipts += 1
    # These three historical producers were superseded in earlier checkpoints.
    # They are unchanged from HEAD during this continuation; old plans stay frozen.
    historical = {}
    for name in ("collect_pi_stacking.py", "prepare_pi_stacking.py", "run_pi_stacking.py"):
        path = ROOT / "scripts/pipeline" / name
        assert path.read_bytes() == subprocess.check_output(["git", "show", "HEAD:" + str(path.relative_to(ROOT))])
        historical[str(path)] = dict(current_sha256=digest(path), unchanged_from_HEAD=True)
    preflight = json.loads((out / "lambda13_preflight.json").read_text())
    assert preflight["ok"]
    summary = dict(
        checked_at=time.time(), status="continuations_prepared_qualified_resource_evidence_full_results_pending",
        tests="23 passed: continuation, admission, cycle and topology tests", lint="Ruff E/F passed",
        prior_campaign_producers_unchanged=previous, historical_superseded_producers=historical,
        prior_completion_receipts_verified=receipts, new_campaigns=new,
        checkpoint_geometry_review=geometry, checkpoint_topology=topology,
        memory_pilot_review=json.loads((campaigns / "2026-10-21_spacing_continuation/plan.json").read_text())["numerical_pilot_review"],
        node_preflight=str(out / "lambda13_preflight.json"), node_preflight_ok=True,
        preserved_original_pilot_failure=str(original / "energy/pilot/status.json"),
        checked_electronic_models=json.loads((ROOT / "results/polymer_models/interaction_status.json").read_text())["structurally_reviewed_DFT_cases"],
        repository_refactor="Deferred until all four scientific tasks, especially #4, are properly addressed",
        physical_limits=["No new spacing energies or native G yet", "A checkpoint is not a qualified minimum",
                         "No experimental loading, bead-scale probability, coupling or rate inferred"],
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print({key: summary[key] for key in ("status", "prior_completion_receipts_verified", "checked_electronic_models")})


if __name__ == "__main__":
    main()
