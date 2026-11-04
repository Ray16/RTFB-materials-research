"""Read checked electronic receipts without losing identity/context or raw values.

This verifies a previously produced physical review and its backing records; it
does not replace fresh SCF, derivative or experimental validation.
"""
import json
import math
from pathlib import Path

from redox.core.hashing import sha256_file
from redox.workflows.molecular_native_namespace import verify


def checked_energy(path,namespace,key,method_id):
    namespace=Path(namespace);report=json.loads(Path(path).read_text())
    if (report.get('status')!='corrected_sampling_energy_independently_checked' or report.get('key')!=key
            or report.get('namespace')!=str(namespace) or report.get('portable_method_id')!=method_id
            or not report.get('actual_checkpoint_frames_verified') or not report.get('occupied_wavefunctions_checked')
            or report.get('native_G_qualified') is not False or report.get('primary_candidate_admission') is not False
            or not report.get('source_sha256')
            or any(sha256_file(p)!=h for p,h in report['source_sha256'].items())):
        raise ValueError('Electronic receipt has changed identity, context, scope or backing')
    spec=verify(namespace);case=Path(spec['case']);plan_path=case/'plan.json';plan=json.loads(plan_path.read_text())
    raw_path=Path(plan['raw_root'])/'experimental_energy.json'
    if not {str(namespace),str(plan_path),str(raw_path)}<=set(report['source_sha256']):
        raise ValueError('Electronic receipt must bind its exact namespace, plan and raw energy')
    actual=json.loads(raw_path.read_text())
    if (report.get('task')!=plan['tasks'][0] or plan['tasks'][0]['key']!=key
            or report.get('experimental_model')!=plan['experimental_model']
            or report.get('electronic_protocol')!=plan['electronic_protocol']
            or actual.get('energy_protocol')!=plan['electronic_protocol']
            or actual.get('experimental_method_id')!=method_id
            or actual.get('converged_smd') is not True or actual.get('converged_gas') is not True
            or actual.get('optimized') is not False
            or (actual.get('charge'),actual.get('mult'))!=(plan['tasks'][0]['charge'],plan['tasks'][0]['mult'])):
        raise ValueError('Electronic receipt differs from its actual graph/state/protocol/convergence')
    for phase in ('smd','gas'):
        value=report.get('E_'+phase+'_eV');raw_value=actual.get('e_'+phase+'_eV')
        if (type(value) not in (int,float) or not math.isfinite(value) or value!=raw_value):
            raise ValueError('Checked electronic scalar does not reproduce its bound raw energy')
    return report
