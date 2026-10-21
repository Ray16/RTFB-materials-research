#!/usr/bin/env python
"""Apply the checked bounded optimizer reset to the two remaining trimer slots."""

import argparse
import os
from pathlib import Path

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api

SCRIPT = Path(__file__).resolve()
SOURCE = ROOT / "scripts/pipeline/continue_polymer_models.py"
BASE = ROOT / "results/polymer_models/campaigns"
KEYS = ("design_pmdi_ethyl__n3__stereo0__open", "design_pmdi_ethyl__n3__stereo0__stacked")


def configure(key):
    original = api.load_module("remaining_trimer_continuation", SOURCE)
    original.SCRIPT = SCRIPT
    original.KEY = key
    original.RECOVERY = BASE / ("2026-10-21_checkpoint_pmdi_" + key.split("__")[-1])
    sources = original.sources

    def source_identity():
        return dict(sources(), **{str(SOURCE): api.digest(SOURCE)})

    original.sources = source_identity
    return original


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "launch", "worker"))
    parser.add_argument("--key", choices=KEYS, default=os.environ.get("REDOX_TRIMER_KEY"))
    parser.add_argument("--role", default="queue_large", choices=("queue_large",))
    args = parser.parse_args()
    if args.key is None:
        parser.error("Supply --key or inherit the launcher's REDOX_TRIMER_KEY")
    os.environ["REDOX_TRIMER_KEY"] = args.key
    module = configure(args.key)
    if args.action == "prepare":
        module.prepare_recovery()
    else:
        runner = module.recovery_runner()
        if args.action == "launch":
            from redox.workflows import resources

            available = resources.free_slots
            # This campaign has exactly one pending optimizer continuation.
            resources.free_slots = lambda *a, **kw: available(*a, **kw)[:1]
            runner.launch(args.role)
        else:
            runner.worker(args.role)
