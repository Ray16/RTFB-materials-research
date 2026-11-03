"""Audit a frozen full-roster sampling snapshot against its checked receipts.

This reuses the maintained component readers. It does not rerun SCF/Hessian QC,
assign an energy to a pending slot, or form a partition from an incomplete set.
"""
import argparse
from collections import Counter
import json
import math
import os
from pathlib import Path

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.sampling_receipts import checked_energy
from redox.validation.sampling_native_receipts import ReceiptCache
from redox.workflows.molecular_native_namespace import verify

CONTEXT = ('task', 'experimental_method_id', 'experimental_model', 'electronic_protocol',
    'native_settings', 'thermal_protocol', 'standard_state', 'temperature_K', 'charge',
    'multiplicity', 'geometry_sha256')


def parent_agreement(native, electronic, geometry_sha256):
    """Compare separate checked executions without replacing either scalar.

    Native pilots can retain an earlier full-orbital electronic execution of
    the same exact input. The public energy inventory can contain a later
    independently checked repeat. Match their physics and report the numerical
    difference; native G always retains its own checked electronic parent.
    """
    difference=abs(native['electronic_E_smd_eV']-electronic['E_smd_eV'])
    if (native['task']!=electronic['task'] or native['experimental_model']!=electronic['experimental_model']
            or native['electronic_protocol']!=electronic['electronic_protocol']
            or native['geometry_sha256']!=geometry_sha256 or not math.isfinite(difference) or difference>1e-8):
        raise ValueError('Native and inventory electronic executions differ in input, context or energy')
    return dict(absolute_energy_difference_eV=difference,comparison_limit_eV=1e-8,
        exact_geometry_and_protocol_match=True,native_electronic_parent_replaced=False)


def roster(snapshot, status):
    rows = snapshot['states']
    keys = [r['key'] for r in rows]
    if (snapshot.get('status') != status or snapshot.get('declared_inputs') != 240
            or len(rows) != 240 or len(set(keys)) != 240
            or snapshot.get('independent_review_errors') or snapshot.get('occupied_ensemble_cycles_qualified')
            or snapshot.get('pending_values_are_null') is not True):
        raise ValueError('A complete distinct 240-slot snapshot with separate ensemble scope is required')
    return {r['key']: r for r in rows}


def source_rosters(snapshot, rows):
    """Bind even pending identities to the actual declared campaign manifest."""
    bindings=snapshot.get('source_sha256',{});manifests=[]
    if not bindings or any(sha256_file(p)!=h for p,h in bindings.items()):
        raise ValueError('Snapshot producer, receipt or declared roster backing changed')
    for name in bindings:
        path=Path(name)
        if path.suffix!='.json':continue
        data=json.loads(path.read_text())
        if not isinstance(data,dict) or data.get('target_states')!=240:continue
        declared=data.get('states',[])
        if (len(declared)!=240 or {r['key'] for r in declared}!=set(rows)
                or data.get('experimental_method_id')!=snapshot['experimental_method_id']):
            raise ValueError('Snapshot slots differ from their source-bound declared roster')
        manifests.append(str(path))
    if not manifests:raise ValueError('The source-bound complete campaign manifest is required')
    return manifests


def audit(electronic_path, native_path):
    before = {str(p): sha256_file(p) for p in (electronic_path, native_path)}
    energy = json.loads(electronic_path.read_text()); native = json.loads(native_path.read_text())
    energies = roster(energy, 'complete_roster_corrected_electronic_review_snapshot')
    natives = roster(native, 'complete_roster_corrected_native_review_snapshot')
    manifests=source_rosters(energy,energies)+source_rosters(native,natives)
    method = energy['experimental_method_id']
    if (energies.keys() != natives.keys() or native['experimental_method_id'] != method
            or energy.get('native_G_qualified') is not False):
        raise ValueError('Electronic/native rosters, model and qualification scopes differ')
    cache = ReceiptCache(); records = []; checked_electronic = {}; sources = dict(before)
    for key, row in energies.items():
        if row.get('G_eV') is not None or row.get('native_G_qualified') is not False:
            raise ValueError('An electronic-only snapshot cannot supply native G')
        if row['status'] != 'electronic_independently_checked':
            if any(row.get(k) is not None for k in ('E_smd_eV', 'E_gas_eV')):
                raise ValueError('Pending electronic slots must stay blank')
            continue
        path = Path(row['review'])
        if sha256_file(path) != row['review_sha256']: raise ValueError('Electronic receipt changed')
        report = checked_energy(path, Path(row['namespace']), key, method)
        if any(row[k] != report[k] for k in ('E_smd_eV', 'E_gas_eV')):
            raise ValueError('Snapshot energy differs from its exact checked electronic component')
        checked_electronic[key] = report
        sources[str(path)] = sha256_file(path)
    for key, row in natives.items():
        if row.get('native_G_qualified') is not True:
            if row['status'] == 'native_independently_checked' or any(row.get(k) is not None
                    for k in ('G_by_cutoff_eV', 'nuclear_by_cutoff_eV')):
                raise ValueError('Pending native slots must stay blank and unqualified')
            continue
        if row['status'] != 'native_independently_checked' or key not in checked_electronic:
            raise ValueError('Every checked native slot needs its exact checked electronic parent')
        path = Path(row['review'])
        if sha256_file(path) != row['review_sha256']: raise ValueError('Native receipt changed')
        report = cache.read(path, Path(row['namespace']), key, method)
        if (any(row[k] != report[k] for k in ('G_by_cutoff_eV', 'nuclear_by_cutoff_eV', 'electronic_E_smd_eV'))
                or row['context'] != {k: report[k] for k in CONTEXT}):
            raise ValueError('Snapshot G, components, graph or physical context differ from their receipts')
        energy_spec=verify(Path(checked_electronic[key]['namespace']))
        geometry=Path(energy_spec['case'])/'input.xyz'
        agreement=parent_agreement(report,checked_electronic[key],sha256_file(geometry))
        sources[str(path)] = sha256_file(path)
        records.append(dict(key=key, review=str(path), G_by_cutoff_eV=report['G_by_cutoff_eV'],
            minimum_frequency_cm=report['minimum_frequency_cm'], maximum_column_error=report['maximum_column_error'],
            actual_SCF_frames_checked=report['actual_SCF_frames_checked'],
            observed_peak_GPU_MiB=report['observed_peak_GPU_MiB'],electronic_inventory_comparison=agreement))
    if (energy['checked_electronic_inputs'] != len(checked_electronic)
            or native['checked_native_inputs'] != len(records)
            or any(sha256_file(p) != h for p, h in before.items())):
        raise ValueError('Frozen snapshot census differs or an input changed during review')
    return dict(status='frozen_sampling_component_receipts_rechecked', declared_inputs=240,
        checked_electronic_inputs=len(checked_electronic), checked_native_inputs=len(records),
        native_by_graph=dict(Counter(r['key'].rsplit('__conf', 1)[0] for r in records)),
        actual_SCF_frames_represented=sum(r['actual_SCF_frames_checked'] for r in records), records=records,
        experimental_method_id=method, context_source_and_component_sums_rechecked=True,
        original_physical_QC_repeated=False, occupied_ensemble_cycles_qualified=False,
        primary_candidate_admission=False, material_admission=False,
        source_manifests=manifests,source_sha256=sources, producer_sha256=sha256_file(__file__))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('electronic', 'native', 'out'): parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
    if args.out.exists(): raise FileExistsError('Preserve previous snapshot audits')
    result = audit(args.electronic.resolve(), args.native.resolve())
    atomic_write_json(args.out, result)
    print(json.dumps({k: result[k] for k in ('status', 'checked_electronic_inputs', 'checked_native_inputs', 'native_by_graph')}))
