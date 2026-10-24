#!/usr/bin/env python
"""Exact-model response resource/derivative pilot with full weighted spin controls."""
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    source = ROOT / 'scripts/validation/qualify_pcm_first_response.py'
    checker_path = ROOT / 'scripts/validation/qualify_streamed_weighted.py'
    checker = api.load_module('weighted_response_full_spin_evidence', checker_path)
    module = api.load_module('preserved_exact_response_pilot', source)
    module.ADAPTER = checker.ADAPTER
    module.checked_controls = checker.checked_controls
    module.main()
    if sys.argv[1] == 'prepare':
        out = Path(sys.argv[sys.argv.index('--out') + 1]).resolve()
        plan = api.read_plan(out)
        for path in (Path(__file__).resolve(), source, checker_path):
            plan['producer_sha256'][str(path)] = api.digest(path)
        plan['execution_profile'] = 'FP64 tiled weighted first response with independently replayed streamed curvature'
        atomic_write_json(out / 'plan.json', plan)


if __name__ == '__main__':
    main()
