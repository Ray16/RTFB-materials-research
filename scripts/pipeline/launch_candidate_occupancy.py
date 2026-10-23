#!/usr/bin/env python
"""Canonical candidate search launch: exact route/geometry review first."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--host", default="lambda13")
    args = parser.parse_args()
    base = args.campaign.resolve()
    reviewer = ROOT / "scripts/validation/review_candidate_occupancy_inputs.py"
    preparer = ROOT / "scripts/pipeline/prepare_candidate_occupancy_searches.py"
    subprocess.run([sys.executable, str(reviewer), "--campaign", str(base)], check=True)
    receipt = base / "input_validation.json"
    checked = json.loads(receipt.read_text())
    if checked["status"] != "exact_candidate_occupancy_inputs_checked":
        raise ValueError("Candidate occupancy inputs are not qualified for a source search")
    api.verify_sources(checked["source_sha256"])
    atomic_write_json(base / "launch_qualification.json", dict(status="feasibility_search_only",
        source_sha256={str(path): api.digest(path) for path in (Path(__file__).resolve(), reviewer, preparer, receipt)},
        prospective_route_status=checked["prospective_route_status"], full_redox_admission=False, synthesis_admission=False))
    api.load_module("qualified_candidate_occupancy_launch", preparer).launch(base, args.host)


if __name__ == "__main__":
    main()
