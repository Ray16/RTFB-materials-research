#!/usr/bin/env python
"""Real measured regression of automatic Fe occupation guesses in SMD."""
import hashlib,json,time
from pathlib import Path
from redox.core.common import ROOT
from redox.core.protocol import ACTIVE_SP,software_versions
from redox.qm.dft import dft_smd
BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10'


def main():
    reference=json.loads((BASE/'refined_fc_smd_occupation_checks/minus.json').read_text())
    geometry=Path(reference['geometry'])
    if hashlib.sha256(geometry.read_bytes()).hexdigest()!=reference['geometry_sha256']:raise RuntimeError('regression input changed')
    t=time.time()
    result=dft_smd(geometry,1,2,backend='gpu',do_opt=False,do_freq=False,do_gas=False,do_smd=True)
    error=result['e_smd_eV']-reference['best']['e_eV']
    report=dict(status='verified' if abs(error)<1e-5 and result['converged_smd'] and
        abs(result.get('smd_spin_contam',99))<.1 and result.get('smd_occupation_guesses')==['rm0','rm1'] else 'failed',
        error_from_independent_seeded_eV=error,independent_energy_eV=reference['best']['e_eV'],
        protocol=ACTIVE_SP,result=result,software=software_versions(),wall_s=time.time()-t,production_modified=False)
    (BASE/'automatic_fc_smd_policy_regression.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({k:v for k,v in report.items() if k!='result'},indent=2),flush=True)
    if report['status']!='verified':raise SystemExit(1)


if __name__=='__main__':main()
