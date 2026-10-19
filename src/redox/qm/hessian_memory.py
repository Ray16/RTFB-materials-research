"""Certified execution budgets for identical native Hessian equations.

The protocol's six-GiB value remains the historical default. A larger execution
budget is recorded separately and requires full-matrix replay controls; it does
not alter geometry, grids, solvent, frequencies, or admission thresholds.
"""
import hashlib
import json
import math
import os
from pathlib import Path
from redox.core.common import ROOT

CERTIFICATE = ROOT / 'config/native_hessian_memory_validation.json'
ALLOCATION_POLICY = 'six-GiB integral batches; qualified larger unrestricted response workspace'
CDS_CERTIFICATE = ROOT / 'config/native_cds_step_validation.json'
CDS_CERTIFICATE_1E6 = ROOT / 'config/native_cds_step_1e6_validation.json'
CERTIFICATE_16 = ROOT / 'config/native_hessian_memory_16_validation.json'
ALLOCATION_POLICY_16 = 'six-GiB integral batches; sixteen-GiB unrestricted response workspace'
PCM_CERTIFICATE = ROOT / 'config/native_pcm_curvature_memory_validation.json'
PCM_MODE = 'host_analytic_tensors_atom_block_transfer'


def pcm_certificate(protocol):
    data=json.loads(PCM_CERTIFICATE.read_text());seal=data.pop('content_sha256')
    if (hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()!=seal
            or data.get('status')!='execution_pcm_curvature_qualified'
            or data.get('protocol')!=protocol or data.get('mode')!=PCM_MODE
            or data.get('scientific_model_changed') is not False or data.get('surface_grid_changed') is not False
            or data.get('software')!={'pyscf':'2.14.0','gpu4pyscf':'1.8.1'}):
        raise RuntimeError('PCM curvature certificate protocol/integrity mismatch')
    controls=data.get('controls',[])
    if len(controls)!=2 or {c.get('mult') for c in controls}!={1,2}:
        raise RuntimeError('PCM replay requires singlet and radical controls')
    for control in controls:
        if control.get('native_gates_passed') is not True:raise RuntimeError('PCM native control gates failed')
        for field in ('relative_matrix_error','maximum_matrix_difference','maximum_thermal_difference_eV'):
            value=control.get(field)
            if not isinstance(value,(int,float)) or not math.isfinite(value) or not 0<=value<=1e-6:
                raise RuntimeError('PCM matrix/thermal replay failed')
        evidence=control.get('evidence_sha256')
        if not isinstance(evidence,dict) or not evidence:raise RuntimeError('PCM evidence missing')
        for name,expected in evidence.items():
            if digest(ROOT/name)!=expected:raise RuntimeError('PCM evidence changed')
    if digest(ROOT/data['isolated_proof_relative_path'])!=data['isolated_proof_sha256']:
        raise RuntimeError('PCM primitive/control proof changed')
    return seal,data


def record_pcm_ok(record,protocol):
    settings=record.get('execution_pcm_curvature')
    if settings is None:return True
    try:
        seal,certificate=pcm_certificate(protocol)
        if (settings.get('qualification_only') is not False or settings.get('production_admission') is not True
                or settings.get('scientific_model_changed') is not False or settings.get('surface_grid_changed') is not False
                or settings.get('mode')!=PCM_MODE or settings.get('qualification_sha256')!=seal
                or settings.get('source_geometry_sha256')!=record['identity']['geometry_sha256']):return False
        for name in ('adapter_sha256','cpu_pcm_library_sha256','gpu_pcm_library_sha256','isolated_proof_sha256'):
            if settings.get(name)!=certificate.get(name):return False
        evidence=settings.get('target_evidence_sha256')
        if not isinstance(evidence,dict) or not evidence:return False
        for name,expected in evidence.items():
            if digest(ROOT/name)!=expected:return False
        return True
    except (OSError,ValueError,KeyError,TypeError,AttributeError,RuntimeError):return False


def derivative_certificate(protocol,step=1e-5):
    if step not in (1e-5,1e-6):
        raise ValueError('Unsupported CDS derivative step')
    path=CDS_CERTIFICATE if step==1e-5 else CDS_CERTIFICATE_1E6
    data = json.loads(path.read_text())
    seal = data.pop('content_sha256')
    if (hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest()!=seal
            or data.get('status')!='execution_derivatives_qualified'
            or data.get('protocol')!=protocol or data.get('cds_step_Bohr')!=step
            or data.get('library_default_cds_step_Bohr')!=1e-4
            or data.get('scientific_model_changed') is not False
            or data.get('software')!={'pyscf':'2.14.0','gpu4pyscf':'1.8.1'}
            or data.get('worst_contact_validation')!=dict(steps_Bohr=[.0005,.001],relative_acceptance=.01,symmetry_acceptance=1e-6)):
        raise RuntimeError('CDS derivative certificate protocol/integrity mismatch')
    controls=data.get('controls',[])
    if len(controls)!=2 or {c.get('mult') for c in controls}!={1,2}:
        raise RuntimeError('CDS derivative qualification requires both spin controls')
    for control in controls:
        if control.get('native_gates_passed') is not True:
            raise RuntimeError('CDS replay did not pass native gates')
        for field in ('relative_matrix_error','maximum_matrix_difference','maximum_thermal_difference_eV'):
            value=control.get(field)
            if not isinstance(value,(int,float)) or not math.isfinite(value) or not 0<=value<=1e-6:
                raise RuntimeError('CDS matrix/thermal replay failed')
        evidence=control.get('evidence_sha256')
        if not isinstance(evidence,dict) or not evidence:
            raise RuntimeError('CDS control evidence missing')
        for name,expected in evidence.items():
            if digest(ROOT/name)!=expected:raise RuntimeError('CDS control evidence changed')
    if step==1e-6:
        evidence=data.get('qualification_evidence_sha256')
        if data.get('target_native_gates_passed') is not True or not isinstance(evidence,dict) or not evidence:
            raise RuntimeError('Finer CDS setting requires independently reviewed target evidence')
        for name,expected in evidence.items():
            if digest(ROOT/name)!=expected:raise RuntimeError('Finer CDS qualification evidence changed')
    return seal


def record_derivatives_ok(record,protocol,natm):
    settings=record.get('execution_derivatives')
    if settings is None:return True
    try:
        step=settings.get('cds_step_Bohr')
        if (settings.get('qualification_only') is not False or settings.get('scientific_model_changed') is not False
                or step not in (1e-5,1e-6)
                or settings.get('qualification_sha256')!=derivative_certificate(protocol,step)
                or settings.get('source_geometry_sha256')!=record['identity']['geometry_sha256']):return False
        indices=settings.get('worst_contact_indices')
        if (not isinstance(indices,list) or len(indices)!=2
                or any(type(i) is not int or not 0<=i<3*natm for i in indices)):return False
        checks=settings.get('worst_contact_checks',[])
        required={(i,s) for i in indices for s in (.0005,.001)}
        observed={(c.get('column'),c.get('step_Bohr')) for c in checks}
        if not required<=observed:return False
        for c in checks:
            error=c.get('relative_error')
            if not isinstance(error,(int,float)) or not math.isfinite(error) or not 0<=error<=.01:return False
        evidence=settings.get('target_evidence_sha256')
        if not isinstance(evidence,dict) or not evidence:return False
        for name,expected in evidence.items():
            if digest(ROOT/name)!=expected:return False
        return True
    except (OSError,ValueError,KeyError,TypeError,AttributeError,RuntimeError):return False


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def certificate(protocol,cap=24):
    path= CERTIFICATE_16 if cap==16 else CERTIFICATE
    policy=ALLOCATION_POLICY_16 if cap==16 else ALLOCATION_POLICY
    data = json.loads(path.read_text())
    seal = data.pop('content_sha256')
    if hashlib.sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest() != seal:
        raise RuntimeError('Hessian memory certificate integrity failed')
    if (data.get('status') != 'execution_budget_qualified' or data.get('protocol') != protocol
            or data.get('software') != {'pyscf': '2.14.0', 'gpu4pyscf': '1.8.1'}
            or data.get('allocation_policy') != policy
            or data.get('caps_GiB') != [cap] or data.get('scientific_model_changed') is not False):
        raise RuntimeError('Hessian memory certificate protocol mismatch')
    controls = data.get('controls', [])
    if len(controls) != 2 or {c.get('mult') for c in controls} != {1, 2}:
        raise RuntimeError('Hessian memory budget requires singlet and radical controls')
    for control in controls:
        for field, limit in [('relative_matrix_error', 1e-6), ('maximum_matrix_difference', 1e-6),
                             ('maximum_thermal_difference_eV', 1e-6)]:
            value = control.get(field)
            if not isinstance(value, (int, float)) or not 0 <= value <= limit:
                raise RuntimeError('Hessian memory replay failed')
        if not control.get('native_gates_passed'):
            raise RuntimeError('Hessian memory control did not pass native gates')
        for name, expected in control['evidence_sha256'].items():
            if digest(ROOT / name) != expected:
                raise RuntimeError('Hessian memory control evidence changed')
    return seal


def execution_budget(protocol, cap=None, qualification=False):
    value = int(cap if cap is not None else os.environ.get('REDOX_HESSIAN_MEMORY_GIB',
                                                         protocol['df_batch_cap_GiB']))
    if value not in (6,16,24):
        raise ValueError('Supported native Hessian budgets are 6 and qualified16/24GiB')
    seal = None if value == 6 or qualification else certificate(protocol,value)
    return dict(cap_GiB=value, protocol_default_cap_GiB=protocol['df_batch_cap_GiB'],
                qualification_sha256=seal, qualification_only=qualification,
                scientific_model_changed=False,allocation_policy=ALLOCATION_POLICY_16 if value==16 else ALLOCATION_POLICY)


def stage_budget_bytes(budget,module_name,function_name):
    """Large persistent UKS response arrays require space; integral batches stay capped."""
    expanded=module_name=='gpu4pyscf.df.hessian.uhf' and function_name=='_get_veff'
    return (budget['cap_GiB'] if expanded else budget['protocol_default_cap_GiB'])*1024**3


def record_budget_ok(record, protocol):
    budget = record.get('execution_memory')
    if budget is None:
        return True  # Historical records used the protocol's fixed six-GiB default.
    try:
        if budget.get('qualification_only') or budget.get('scientific_model_changed') is not False:
            return False
        if budget.get('protocol_default_cap_GiB') != protocol['df_batch_cap_GiB']:
            return False
        if budget.get('cap_GiB') == 6:
            return budget.get('qualification_sha256') is None
        value=budget.get('cap_GiB')
        return (value in (16,24)
                and budget.get('allocation_policy')==(ALLOCATION_POLICY_16 if value==16 else ALLOCATION_POLICY)
                and budget.get('qualification_sha256')==certificate(protocol,value))
    except (OSError, ValueError, KeyError, RuntimeError):
        return False
