#!/usr/bin/env python
"""Report stored physical-QC values for the twelve frozen spacing conformers."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

PRODUCER_BYTES=Path(__file__).read_bytes()


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(source,out):
    if out.exists():raise FileExistsError('Preserve previous QC audits')
    data=json.loads(source.read_text());rows=[];proof={str(source):sha(source)}
    for case in data['cases']:
        for conformer in case['conformers']:
            receipt=Path(conformer['thermal_receipt']);r=json.loads(receipt.read_text());status=Path(r['status_path'])
            if sha(status)!=r['status_sha256']:raise ValueError('Thermal receipt no longer matches its status')
            record=json.loads(status.read_text())
            if record['status']!='state_records_verified' or record['thermal_QC']!='ok':raise ValueError('Unqualified native conformer')
            for name,expected in record['output_sha256'].items():
                if sha(Path(name))!=expected:raise ValueError('Source output changed: '+name)
            paths=[Path(p) for p in record['output_sha256'] if Path(p).parent.name=='native_thermal' and p.endswith('.json')]
            partitions=[Path(p) for p in record['output_sha256'] if Path(p).parent.name=='projected_thermal' and p.endswith('.json')]
            if len(paths)!=1 or len(partitions)!=1:raise ValueError('One exact native record and projected partition required')
            native=json.loads(paths[0].read_text());part=json.loads(partitions[0].read_text())
            gradient=native['gradient_max'];asymmetry=native['maximum_asymmetry'];frequencies=part['frequencies_cm']
            error=max(c['relative_error'] for c in native['derivative_checks'])
            values=[gradient,asymmetry,error,*frequencies]
            if (not all(math.isfinite(x) for x in values) or not 0<=gradient<=.00045
                    or not 0<=asymmetry<=1e-6 or not 0<=error<=.01 or min(frequencies)<-5):
                raise ValueError('Stored native physical QC does not pass')
            rows.append(dict(key=conformer['key'],occupied_sites=case['case']['occupied_sites'],
                             maximum_gradient_Ha_Bohr=gradient,maximum_Hessian_asymmetry=asymmetry,
                             lowest_projected_internal_frequency_cm=min(frequencies),
                             largest_sampled_column_relative_error=error,native_record=str(paths[0]),
                             native_sha256=sha(paths[0]),partition=str(partitions[0]),partition_sha256=sha(partitions[0])))
            proof.update({str(receipt):sha(receipt),str(status):sha(status),str(paths[0]):sha(paths[0]),str(partitions[0]):sha(partitions[0])})
    if len(rows)!=12 or len({r['key'] for r in rows})!=12:raise ValueError('Exactly twelve distinct baseline conformers required')
    result=dict(status='twelve_spacing_native_QC_values_reconstructed_from_checked_records',checked=time.time(),
                rows=rows,maximum_gradient_range_Ha_Bohr=[min(r['maximum_gradient_Ha_Bohr'] for r in rows),max(r['maximum_gradient_Ha_Bohr'] for r in rows)],
                source_sha256=proof,producer_sha256=hashlib.sha256(PRODUCER_BYTES).hexdigest(),
                scope='Read-only extraction of already checked physical-QC values. No new optimization, Hessian or conformer convergence claim.')
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    r=audit(a.source.resolve(),a.out.resolve());print(json.dumps(dict(status=r['status'],gradient_range=r['maximum_gradient_range_Ha_Bohr'])))
