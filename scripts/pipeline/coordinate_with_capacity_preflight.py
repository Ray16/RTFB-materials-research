#!/usr/bin/env python
"""Resume declared smooth energies with CPU-only, bound-aware placement.

The original controller, scientific plans and worker bytes remain unchanged.
Only process-local preflight/placement functions are adapted, with their source
hashes checked throughout execution. Uncertain existing submissions stay held.
"""
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def main():
    upstream=ROOT/'scripts/pipeline/coordinate_smooth_electronic_batch.py'
    adapter_path=ROOT/'scripts/pipeline/dispatch_smooth_energy_capacity.py'
    proof={str(p):api.digest(p) for p in (Path(__file__).resolve(),upstream,adapter_path)}
    original_load,original_verify=api.load_module,api.verify_sources
    adapter=original_load('bound_aware_energy_capacity_adapter',adapter_path)
    def verify(mapping):original_verify(proof);original_verify(mapping)
    def load(name,path):
        verify({});module=original_load(name,path)
        if Path(path).resolve()==ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py':
            module.placement_module=adapter.capacity_placement
            module.preflight=adapter.capacity_preflight(module)
        return module
    api.load_module,api.verify_sources=load,verify
    module=load('unchanged_electronic_controller_bound_aware_placement',upstream)
    module.SCRIPT=Path(__file__).resolve()
    sys.argv=[str(Path(__file__).resolve()),*sys.argv[1:]];module.main()


if __name__=='__main__':main()
