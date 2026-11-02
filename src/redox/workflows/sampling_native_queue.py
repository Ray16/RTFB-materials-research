"""Full-roster native collection and growing, source-bound preparation.

The 240 original slots remain visible, including those whose energy is pending.
Only independently checked energies become native inputs; the largest complete
native pilot is mandatory before preparing any followup.
"""
import json
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.sampling_native_receipts import ReceiptCache
from redox.validation.sampling_receipts import checked_energy
from redox.workflows.molecular_native_namespace import verify, load

SOURCE_BYTES = Path(__file__).read_bytes()


def roster(manifest_path):
    manifest = json.loads(Path(manifest_path).read_text())
    rows = manifest['states']
    if manifest.get('target_states') != 240 or len(rows) != 240 or len({r['key'] for r in rows}) != 240:
        raise ValueError('The complete 240 distinct energy slots are required')
    return manifest


def checked_stream(manifest, stream):
    rows = stream['states']
    expected = {r['key'] for r in manifest['states']}
    if (stream.get('status') != 'complete_roster_corrected_electronic_review_snapshot'
            or stream.get('declared_inputs') != 240 or len(rows) != 240
            or {r['key'] for r in rows} != expected
            or stream.get('experimental_method_id') != manifest['experimental_method_id']
            or stream.get('independent_review_errors')):
        raise ValueError('A complete matching, independently reviewed electronic stream is required')
    return [r for r in rows if r['status'] == 'electronic_independently_checked']


def collect(manifest_path, native_rows, reviews, cache=None):
    """Retain all slots; reconstruct only producers with closed execution."""
    manifest = roster(manifest_path)
    cache = cache or ReceiptCache()
    reviews = Path(reviews)
    reviews.mkdir(parents=True, exist_ok=True)
    native = {r['key']: r for r in native_rows}
    if len(native) != len(native_rows) or not native.keys() <= {r['key'] for r in manifest['states']}:
        raise ValueError('Native rows must be distinct declared sampling slots')
    reviewer = load(ROOT / 'scripts/validation/review_corrected_sampling_native.py', 'queue_native_review')
    owners = load(ROOT / 'scripts/pipeline/observe_archived_native_namespace.py', 'queue_native_owners')
    records, errors = [], []
    for row in manifest['states']:
        key = row['key']
        record = dict(key=key, status='native_not_prepared', native_G_qualified=False,
                      G_by_cutoff_eV=None, nuclear_by_cutoff_eV=None)
        if key not in native:
            records.append(record)
            continue
        prepared = native[key]
        namespace = Path(prepared['namespace'])
        if sha256_file(namespace) != prepared['namespace_sha256']:
            raise ValueError('Native namespace changed')
        case = Path(json.loads(namespace.read_text())['case'])
        receipt = reviews / (key + '.json')
        record['namespace'] = str(namespace)
        try:
            if not receipt.exists():
                status_path = case / 'status.json'
                if not status_path.exists():
                    record['status'] = 'native_prepared'
                    records.append(record)
                    continue
                stage = json.loads(status_path.read_text())
                record['status'] = stage['status']
                if stage['status'] != 'corrected_sampling_native_finished_independent_review_pending':
                    records.append(record)
                    continue
                guard = case / 'residency_guard/owner.json'
                if (owners.alive(stage['owner']) is not False
                        or not guard.exists() or not json.loads(guard.read_text()).get('guard_stopped')
                        or not list(case.parent.glob('placement*/completion.json'))):
                    record['status'] = 'native_execution_closure_pending'
                    records.append(record)
                    continue
                atomic_write_json(receipt, reviewer.review(namespace))
            report = cache.read(receipt, namespace, key, manifest['experimental_method_id'])
            record.update(status='native_independently_checked', native_G_qualified=True,
                          review=str(receipt), review_sha256=sha256_file(receipt),
                          G_by_cutoff_eV=report['G_by_cutoff_eV'],
                          nuclear_by_cutoff_eV=report['nuclear_by_cutoff_eV'],
                          electronic_E_smd_eV=report['electronic_E_smd_eV'],
                          context={k: report[k] for k in ('task', 'experimental_method_id', 'experimental_model',
                              'electronic_protocol', 'native_settings', 'thermal_protocol', 'standard_state',
                              'temperature_K', 'charge', 'multiplicity', 'geometry_sha256')})
        except Exception as error:
            record.update(status='native_independent_review_failed_preserved', error=repr(error))
            errors.append(dict(key=key, error=repr(error)))
        records.append(record)
    from redox.validation.execution_review import freeze_source
    from redox.validation import sampling_native_receipts
    producer = freeze_source(reviews.parent, Path(__file__), SOURCE_BYTES)
    receipt_reader = freeze_source(reviews.parent, Path(sampling_native_receipts.__file__))
    return dict(status='complete_roster_corrected_native_review_snapshot', checked=time.time(),
                declared_inputs=240, checked_native_inputs=sum(r['native_G_qualified'] for r in records),
                states=records, independent_review_errors=errors,
                experimental_method_id=manifest['experimental_method_id'],
                occupied_ensemble_cycles_qualified=False, pending_values_are_null=True,
                source_sha256={str(manifest_path): sha256_file(manifest_path),
                    str(producer): sha256_file(producer), str(receipt_reader): sha256_file(receipt_reader),
                    **{r['review']: r['review_sha256'] for r in records if r.get('native_G_qualified')}})


def prepare_ready(manifest_path, stream, rows, native_pilot, resource_review, archive, output):
    from redox.workflows.sampling_native import prepare
    if native_pilot is None:
        raise ValueError('The largest complete native pilot is required before followup preparation')
    manifest = roster(manifest_path)
    eligible = checked_stream(manifest, stream)
    expected = {r['key']: r for r in manifest['states']}
    existing = {r['key']: r for r in rows}
    if len(existing) != len(rows) or not existing.keys() <= expected.keys():
        raise ValueError('Distinct declared native preparation required')
    for row in eligible:
        if row['key'] in existing:
            continue
        path = Path(row['review'])
        if sha256_file(path) != row['review_sha256']:
            raise ValueError('Electronic stream receipt changed')
        report = checked_energy(path, Path(row['namespace']), row['key'], manifest['experimental_method_id'])
        # Recovery namespaces may change execution ownership, never the task,
        # scored geometry, physical model or scientific settings.
        original = verify(Path(expected[row['key']]['namespace']))
        actual = verify(Path(row['namespace']))
        before = json.loads((Path(original['case']) / 'plan.json').read_text())
        after = json.loads((Path(actual['case']) / 'plan.json').read_text())
        if (report['task'] != before['tasks'][0]
                or any(before[k] != after[k] for k in ('experimental_model', 'electronic_protocol'))
                or sha256_file(Path(original['case']) / 'input.xyz') != sha256_file(Path(actual['case']) / 'input.xyz')):
            raise ValueError('A recovery changed the declared scientific input')
        directory = Path(output) / 'states' / row['key']
        if directory.exists():
            raise FileExistsError('Unregistered native preparation is preserved for explicit reconciliation')
        prepared = prepare(path, resource_review, archive, directory, native_pilot)
        rows.append({k: prepared[k] for k in ('key', 'namespace', 'namespace_sha256')})
        # Persist each completed preparation before attempting the next input.
        atomic_write_json(Path(output) / 'prepared_current.json', dict(declared_sampling_inputs=240, states=rows))
    return rows


def extend_peers(path, rows):
    """Serialize append-only extensions shared with the electronic coordinator."""
    import fcntl
    path = Path(path)
    with path.with_name(path.name + '.append.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        peers = json.loads(path.read_text())
        existing = {r['namespace']: r['namespace_sha256'] for r in peers}
        if len(existing) != len(peers):
            raise ValueError('Duplicate storage peers')
        for row in rows:
            namespace, digest = row['namespace'], row['namespace_sha256']
            if sha256_file(namespace) != digest or existing.get(namespace, digest) != digest:
                raise ValueError('A storage-peer namespace changed')
            if namespace not in existing:
                peers.append(dict(namespace=namespace, namespace_sha256=digest))
                existing[namespace] = digest
        atomic_write_json(path, peers)
    return peers
