"""Read native sampling values with their exact electronic and thermal context.

A receipt proves a previously completed independent reconstruction. Reading it
does not perform a new Hessian calculation or establish ensemble convergence.
"""
import copy
import json
import math
from pathlib import Path

from redox.core.hashing import sha256_file
from redox.validation.sampling_receipts import checked_energy
from redox.workflows.molecular_native_namespace import verify


def checked_native(path, namespace, key, method_id):
    namespace = Path(namespace)
    report = json.loads(Path(path).read_text())
    if (report.get('status') != 'corrected_sampling_native_independently_reconstructed'
            or report.get('namespace') != str(namespace) or report.get('key') != key
            or report.get('experimental_method_id') != method_id
            or report.get('native_G_qualified') is not True
            or report.get('primary_candidate_admission') is not False
            or report.get('conformer_ensemble_qualified') is not False
            or not report.get('source_sha256')
            or any(sha256_file(p) != h for p, h in report['source_sha256'].items())):
        raise ValueError('Native receipt identity, qualification or backing differs')
    spec = verify(namespace)
    case = Path(spec['case'])
    plan_path = case / 'plan.json'
    plan = json.loads(plan_path.read_text())
    partition_path = Path(plan['raw_root']) / 'experimental_partition.json'
    electronic_path = Path(plan['electronic_review'])
    required = {str(p) for p in (namespace, plan_path, partition_path, electronic_path)}
    if not required <= report['source_sha256'].keys():
        raise ValueError('Native receipt must bind its namespace, plan and both components')
    partition = json.loads(partition_path.read_text())
    preliminary = json.loads(electronic_path.read_text())
    electronic = checked_energy(electronic_path, Path(preliminary['namespace']), key, method_id)
    energy_spec = verify(Path(electronic['namespace']))
    for directory in (case, Path(energy_spec['case'])):
        if (directory / 'residency_guard/external_interruption.json').exists():
            raise ValueError('A residency interruption invalidates the native component or its electronic parent')
    if (plan.get('role') != 'corrected_sampling_native'
            or report.get('task') != plan['tasks'][0] or report['task'] != electronic['task']
            or report['task']['key'] != key
            or (report.get('charge'), report.get('multiplicity')) != (0, 1)
            or (report['task']['charge'], report['task']['mult']) != (0, 1)
            or report.get('experimental_model') != plan['experimental_model']
            or report['experimental_model'] != electronic['experimental_model']
            or report.get('electronic_protocol') != plan['electronic_protocol']
            or report['electronic_protocol'] != electronic['electronic_protocol']
            or report.get('native_settings') != plan['experimental_native_settings']
            or report.get('thermal_protocol') != plan['thermal_protocol']
            or report.get('standard_state') != plan['thermal_protocol']['pressure']
            or report.get('temperature_K') != plan['thermal_protocol']['temperature_K']
            or report.get('geometry_sha256') != plan['scored_geometry_sha256']
            or sha256_file(case / 'input.xyz') != report['geometry_sha256']
            or sha256_file(Path(energy_spec['case']) / 'input.xyz') != report['geometry_sha256']
            or report.get('thermal_source') != partition['nuclear']
            or partition.get('experimental_method_id') != method_id
            or partition.get('derivative_checks_passed') is not True
            or partition.get('primary_candidate_admission') is not False
            or report.get('electronic_E_smd_eV') != electronic['E_smd_eV']
            or report.get('spin_partition_eV') != partition.get('spin_partition_eV')
            or report['spin_partition_eV'] != 0.):
        raise ValueError('Native receipt loses exact graph, state, geometry or component context')
    for name in ('G_by_cutoff_eV', 'nuclear_by_cutoff_eV'):
        if set(report.get(name, {})) != {'25', '50', '100'}:
            raise ValueError('Every declared rotor setting is required')
    for name, low, high in (('minimum_frequency_cm', -5, float('inf')),
                            ('maximum_column_error', 0, .01),
                            ('observed_peak_GPU_MiB', 1, 31744),
                            ('peak_host_RSS_bytes', 1, 80 * 2**30),
                            ('actual_SCF_frames_checked', 1, float('inf'))):
        value = report.get(name)
        if type(value) not in (int, float) or not math.isfinite(value) or not low <= value <= high:
            raise ValueError('Native physical/resource review quantities are absent or outside their gates')
    for cutoff in ('25', '50', '100'):
        nuclear = report['nuclear_by_cutoff_eV'][cutoff]
        total = report['G_by_cutoff_eV'][cutoff]
        if (any(type(v) not in (int, float) or not math.isfinite(v) for v in (nuclear, total))
                or nuclear != partition['nuclear']['nuclear_by_cutoff_eV'][cutoff]
                or total != partition['G_by_cutoff_eV'][cutoff]
                or abs(total - electronic['E_smd_eV'] - nuclear - report['spin_partition_eV']) > 1e-8):
            raise ValueError('Native electronic/nuclear sum does not reproduce the raw partition')
    return report


class ReceiptCache:
    """Hash once, then invalidate on any backing-file identity or stat change.

    Native receipts contain large immutable Hessian/checkpoint files. A running
    collector checks their metadata on every read and rehashes after a change.
    Small definitions/values are always hashed, including same-timestamp writes.
    A fresh controller always performs the full hash validation again.
    """
    def __init__(self):
        self.entries = {}

    @staticmethod
    def signatures(paths):
        result = {}
        for path in paths:
            stat = Path(path).stat()
            result[str(path)] = (stat.st_dev, stat.st_ino, stat.st_size,
                                 stat.st_mtime_ns, stat.st_ctime_ns,
                                 sha256_file(path) if stat.st_size <= 2**20 else None)
        return result

    def read(self, path, namespace, key, method_id):
        identity = (str(path), str(namespace), key, method_id)
        old = self.entries.get(identity)
        if old is not None and not any(p.exists() for p in old[2]) and self.signatures(old[0]) == old[0]:
            return copy.deepcopy(old[1])
        preliminary = json.loads(Path(path).read_text())
        spec = json.loads(Path(namespace).read_text())
        plan_path = Path(spec['case']) / 'plan.json'
        plan = json.loads(plan_path.read_text())
        electronic_path = Path(plan['electronic_review'])
        electronic = json.loads(electronic_path.read_text())
        energy_namespace = Path(electronic['namespace'])
        energy_spec = json.loads(energy_namespace.read_text())
        energy_plan_path = Path(energy_spec['case']) / 'plan.json'
        energy_plan = json.loads(energy_plan_path.read_text())
        paths = {str(path), str(namespace), str(plan_path), str(electronic_path),
                 str(energy_namespace), str(energy_plan_path),
                 *preliminary['source_sha256'], *electronic['source_sha256'],
                 *plan['source_sha256'], *energy_plan['source_sha256']}
        before = self.signatures(paths)
        report = checked_native(path, namespace, key, method_id)
        after = self.signatures(paths)
        if before != after:
            raise ValueError('Receipt backing changed during hash validation')
        forbidden = [Path(spec['case']) / 'residency_guard/external_interruption.json',
                     Path(energy_spec['case']) / 'residency_guard/external_interruption.json']
        self.entries[identity] = (after, copy.deepcopy(report), forbidden)
        return report
