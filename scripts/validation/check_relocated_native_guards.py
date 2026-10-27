#!/usr/bin/env python
"""Adverse-input checks for the explicit relocated native-record reader."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import shutil
import tempfile


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original-source-root',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    from redox.core.common import ROOT
    from redox.qm import native_thermal as native
    from redox.qm.relocated_native import relocated_native_reader
    source=ROOT/'calcs/dft/formaldehyde/rest'
    raw_path=next((source/'native_thermal').glob('*.json'))
    original=raw_path.read_bytes();raw=json.loads(original)
    original_reader=native._load_native_surface_record
    if original_reader(source,0,1) is not None:
        raise RuntimeError('Expected the unadapted reader to refuse the moved certificate paths')
    checks=[]
    with relocated_native_reader(args.original_source_root):
        if native._load_native_surface_record(source,0,1) is None:
            raise RuntimeError('Byte-identical relocated input failed admission')
        mutations=[('gradient',lambda r:r.update(gradient_max=1.)),
                   ('stationarity',lambda r:r.update(stationary_minimum=False)),
                   ('spin_contamination',lambda r:r['smd_spin_diagnostics'].update(spin_contam=.5)),
                   ('missing_derivative_checks',lambda r:r.update(derivative_checks=[])),
                   ('identity',lambda r:r['identity'].update(mult=3)),
                   ('certificate_hash',lambda r:r['implementation_controls'][0].update(sha256='changed')),
                   ('imaginary_frequency',lambda r:r['frequencies_cm'].__setitem__(0,-1000.)),
                   ('software',lambda r:r['software'].update(gpu4pyscf='changed')),
                   ('asymmetry',lambda r:r.update(maximum_asymmetry=.1)),
                   ('undeclared_original_prefix',lambda r:r['implementation_controls'][0].update(path='/unrelated/status.json'))]
        for name,mutate in mutations:
            with tempfile.TemporaryDirectory(dir=ROOT/'calcs',prefix='adverse_native_') as work:
                state=Path(work)/'rest';state.mkdir()
                for f in ('result.json','opt.xyz'):shutil.copyfile(source/f,state/f)
                shutil.copytree(source/'native_thermal',state/'native_thermal')
                record=copy.deepcopy(raw);mutate(record)
                (state/'native_thermal'/raw_path.name).write_text(json.dumps(record))
                rejected=native._load_native_surface_record(state,0,1) is None
                checks.append(dict(case=name,rejected=rejected))
                if not rejected:raise RuntimeError('Altered physical evidence was admitted: '+name)
    if native._load_native_surface_record is not original_reader or raw_path.read_bytes()!=original:
        raise RuntimeError('Reader or source record was not preserved')
    result=dict(status='relocated_native_reader_adverse_inputs_rejected',checks=checks,
                original_reader_restored=True,original_raw_record_unchanged=True,
                producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    with args.out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(status=result['status'],negative_checks=len(checks))))


if __name__=='__main__':main()
