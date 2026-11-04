"""Prospective T4 optimization of the exact, already prepared chloride repairs.

The fixed-gradient witness admits one optimization pilot. Its closed physical
and resource review is required before any additional repair is admitted. No
Hessian, native free energy, or general molecular GPU class is qualified here.
Keep top-level dependencies standard-library-only for the frozen worker helper.
"""
import argparse
import json
import os
from pathlib import Path
import shutil

PROFILE = 'primary_chloride_T4_geometry_pilot_v1'
METHOD = dict(xc='r2scan', basis='def2-svp', disp='d4', solvent='acetonitrile',
    solvent_model='primary SMD with legacy CDS', SCF_conv_tol=1e-12,
    SCF_conv_tol_grad=1e-8, SCF_max_cycle=200, grid_level=3, gradient_grid_response=True)
CONVERGENCE = dict(convergence_energy=1e-9, convergence_grms=5e-6,
    convergence_gmax=1e-5, convergence_drms=1e-4, convergence_dmax=2e-4)
SMILES = 'Cc1ccc(C[n+]2ccc(-c3cc[n+](C)cc3)cc2)cc1.[Cl-]'


def validate_plan(plan):
    expected = dict(role='primary_counterion_geometry_control', geometry_resource_profile=PROFILE,
        geometry_memory_pilot=True, atoms=42, nao=412, charge=1, mult=1,
        parent_formal_core_charge=2, counterion_charge=-1, smiles=SMILES,
        maximum_host_RSS_GiB=32, minimum_free_GiB=32, minimum_GPU_memory_MiB=15360,
        maximum_peak_GPU_MiB=12288, maximum_wall_seconds=14400,
        maximum_durable_output_bytes=2*2**30, Hessian_computed=False,
        native_G_qualified=False, material_admission=False)
    seed = plan.get('seed', {})
    if (any(plan.get(k) != v for k, v in expected.items()) or plan.get('method') != METHOD
            or plan.get('optimizer_convergence') != CONVERGENCE
            or plan.get('source_id') not in ('N0_below', 'N0_outer')
            or seed.get('source') != 'verified negative native mode' or seed.get('sign') not in (-1, 1)
            or seed.get('maximum_atom_displacement_A') != .1
            or not isinstance(plan.get('fixed_gradient_review'), str)):
        raise ValueError('Only the declared exact T4 chloride geometry repairs are admitted')
    return expected


def device_for_plan(plan, gpu):
    validate_plan(plan)
    from redox.workflows.resources import inspect_gpu
    device = inspect_gpu(gpu)
    if (device['memory_MiB'] < 15360 or float(device['compute_capability']) != 7.5
            or 'T4' not in device['name'].split()):
        raise RuntimeError('The independent geometry resource profile specifically tests a free T4')
    return device


def fixed_gradient_witness(path):
    from redox.core.hashing import sha256_file
    from redox.workflows.molecular_native_namespace import verify
    report = json.loads(Path(path).read_text())
    if (report.get('status') != 'T4_primary_chloride_fixed_gradient_independently_checked'
            or report.get('method') != METHOD or report.get('atoms') != 42
            or (report.get('charge'), report.get('mult')) != (1, 1)
            or report.get('actual_SCF_frames_checked') != 1
            or report.get('fixed_gradient_resource_control_passed') is not True
            or report.get('full_geometry_resource_class_qualified') is not False
            or report.get('Hessian_class_qualified') is not False
            or report.get('native_G_qualified') is not False
            or report.get('material_admission') is not False
            or any(sha256_file(p) != h for p, h in report['source_sha256'].items())):
        raise ValueError('A source-backed independently checked fixed-gradient T4 witness is required')
    verify(Path(report['namespace']))
    return report


def geometry_witness(path):
    """A report alone cannot release followups without its exact as-used plan."""
    from redox.core.hashing import sha256_file
    from redox.workflows.molecular_native_namespace import verify
    report = json.loads(Path(path).read_text())
    spec = verify(Path(report['namespace']))
    plan = json.loads((Path(spec['case'])/'plan.json').read_text())
    validate_plan(plan)
    if (report.get('status') != 'primary_counterion_geometry_independently_checked'
            or report.get('geometry_resource_profile') != PROFILE
            or report.get('geometry_resource_class_qualified') is not True
            or report.get('Hessian_class_qualified') is not False
            or report.get('native_G_qualified') is not False
            or report.get('material_admission') is not False
            or any(sha256_file(p) != h for p, h in report['source_sha256'].items())):
        raise ValueError('The first T4 optimization needs closed physical and resource review')
    fixed_gradient_witness(Path(plan['fixed_gradient_review']))
    return report


def prepare(namespace, fixed_review, output, pilot_review=None):
    from redox.core.common import ROOT
    from redox.core.hashing import sha256_file
    from redox.core.io import atomic_write_json
    from redox.workflows.molecular_native_namespace import verify
    source = verify(namespace)
    origin = Path(source['case'])
    old = json.loads((origin/'plan.json').read_text())
    fixed_gradient_witness(fixed_review)
    if pilot_review is None:
        if old['source_id'] != 'N0_below' or old.get('seed', {}).get('sign') != -1:
            raise ValueError('The first optimization pilot is the prepared N0_below minus repair')
    else:
        geometry_witness(pilot_review)
    if any((origin/n).exists() for n in ('status.json', 'residency_guard', 'host_guard')):
        raise ValueError('Only an unstarted prepared repair may get a new execution namespace')
    if output.exists() or any(output.is_relative_to(Path(source[k])) for k in ('prior_root', 'archive')):
        raise FileExistsError('Use a fresh durable namespace outside scientific source trees')
    plan = dict(old, geometry_resource_profile=PROFILE, geometry_memory_pilot=True,
        minimum_GPU_memory_MiB=15360, maximum_peak_GPU_MiB=12288, Hessian_computed=False,
        fixed_gradient_review=str(output/'sources/fixed_gradient_review.json'),
        resource_pilot_review=None if pilot_review is None else str(output/'sources/geometry_pilot_review.json'),
        execution_source_namespace=str(namespace), source_sha256=dict(old['source_sha256']))
    validate_plan(plan)
    case = output/'case'; sources = output/'sources'
    case.mkdir(parents=True); sources.mkdir()
    for name in ('dft', 'locks', 'repairs'): (case/name).mkdir()
    for path in Path(source['sources']).iterdir():
        if path.is_file(): shutil.copyfile(path, sources/path.name)
    for path, name in ((ROOT/'scripts/validation/run_counterion_geometry.py', 'worker.py'),
            (Path(__file__), 'resource_profile.py'), (fixed_review, 'fixed_gradient_review.json'),
            (namespace, 'source_namespace.json'), (origin/'plan.json', 'source_plan.json')):
        shutil.copyfile(path, sources/name)
    if pilot_review is not None: shutil.copyfile(pilot_review, sources/'geometry_pilot_review.json')
    shutil.copyfile(origin/'input.xyz', case/'input.xyz')
    plan['source_sha256'].update({str(p): sha256_file(p) for p in sources.iterdir() if p.is_file()})
    plan['source_sha256'][str(case/'input.xyz')] = sha256_file(case/'input.xyz')
    atomic_write_json(case/'plan.json', plan)
    atomic_write_json(output/'namespace.json', dict(source, case=str(case), sources=str(sources),
        plan_sha256=sha256_file(case/'plan.json')))
    verify(output/'namespace.json')
    return dict(status='T4_geometry_only_repair_prepared', namespace=str(output/'namespace.json'),
        pilot_required=pilot_review is None, quantum_jobs_submitted=0)


def launch(namespace, gpu, output):
    from redox.workflows.molecular_native_namespace import verify
    from redox.workflows.bounded_molecular import launch as execute
    spec = verify(namespace); plan = json.loads((Path(spec['case'])/'plan.json').read_text())
    validate_plan(plan); fixed_gradient_witness(Path(plan['fixed_gradient_review']))
    if plan['resource_pilot_review'] is not None: geometry_witness(Path(plan['resource_pilot_review']))
    elif plan['source_id'] != 'N0_below' or plan['seed']['sign'] != -1:
        raise ValueError('Additional repairs need the independent full geometry pilot')
    return execute(namespace, gpu, output, caller_source=Path(__file__))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'launch'))
    for name in ('spec', 'fixed-review', 'pilot-review', 'out'): parser.add_argument('--'+name, type=Path)
    parser.add_argument('--gpu', type=int); args = parser.parse_args()
    if args.action == 'prepare':
        from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
        with allocation_guard():
            host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
        result = prepare(args.spec.resolve(), args.fixed_review.resolve(), args.out.resolve(),
            None if args.pilot_review is None else args.pilot_review.resolve())
    else: result = launch(args.spec.resolve(), args.gpu, args.out.resolve())
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_sha256', 'argv')}))
