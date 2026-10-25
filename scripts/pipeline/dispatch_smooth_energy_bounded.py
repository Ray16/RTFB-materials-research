#!/usr/bin/env python
"""Use the existing guarded energy placement with the source-bound RAM branch."""
from pathlib import Path
from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


if __name__ == '__main__':
    module = api.load_module('smooth_energy_measured_RAM_placement', ROOT / 'scripts/pipeline/dispatch_smooth_energy_pilot.py')
    module.SCRIPT = Path(__file__).resolve()
    module.WORKER = ROOT / 'scripts/pipeline/smooth_polymer_energy_bounded.py'
    module.main()
