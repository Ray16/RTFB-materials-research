#!/usr/bin/env python
"""Verify the production Fe seed policy against the independently repaired Fc+ energy."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from redox.core.common import DFT, ROOT
from redox.qm.dft import dft_smd


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--backend',choices=['gpu','cpu'],default='gpu');a=ap.parse_args()
    if a.out.exists():ap.error('previous validation exists; preserve it')
    xyz=DFT/'ferrocene/ox/opt.xyz';start=time.time()
    expected=json.loads((ROOT/'results/validation/physical_thermal_spin_2026-10-10/fc_gas_repair.json').read_text())
    result=dft_smd(xyz,1,2,do_opt=False,do_freq=False,do_smd=False,do_gas=True,backend=a.backend)
    error=result['e_gas_eV']-expected['new_energy_eV']
    record=dict(status='passed' if result['converged_gas'] and abs(error)<1e-6 else 'failed',
                difference_eV=error,elapsed_seconds=time.time()-start,
                source_sha256=hashlib.sha256((ROOT/'src/redox/qm/dft.py').read_bytes()).hexdigest(),
                geometry_sha256=hashlib.sha256(xyz.read_bytes()).hexdigest(),result=result,
                production_data_modified=False)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)
    if record['status']!='passed':raise SystemExit(1)


if __name__=='__main__':main()
