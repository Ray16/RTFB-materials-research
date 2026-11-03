"""Separate proven closed resource refusals from bounded prelaunch failures.

This is read-only accounting. It retires no claim and admits no scientific
result. Started attempts and uncertain ownership must remain preserved.
"""
import json
from pathlib import Path

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.workflows.submission_refusal import refusal_branch

RESOURCE_MESSAGES=('GPU occupied; no screening control started',
                   'Control RAM or output headroom unavailable')


def launcher_versions(output):
    """Previously frozen launcher versions retain the original traceback lines."""
    return sorted({p for name in ('sources','review_sources','source_snapshots')
                   for p in (Path(output)/name).glob('*/screening_control.py')})


def prelaunch_accounting(namespace,alive,source_versions=()):
    root=Path(namespace).parent
    if (root/'case/status.json').exists() or (root/'case/residency_guard/owner.json').exists():
        raise ValueError('Never reopen a started scientific or guard attempt through prelaunch accounting')
    source=ROOT/'src/redox/workflows/screening_control.py'
    versions=[source,*[Path(p) for p in source_versions if Path(p)!=source]]
    hashes={str(p):sha256_file(p) for p in versions}
    if any(p!=source and p.parent.name!=hashes[str(p)] for p in versions):
        raise ValueError('Historical launcher source must retain its content-addressed snapshot')
    failures=0;records=[]
    for path in sorted((root/'dispatch_claim').glob('try_*.json')):
        receipt=json.loads(path.read_text())
        if alive(receipt['owner']) is not False:
            raise ValueError('Every exact prelaunch owner must be confirmed closed')
        placement=Path(receipt['placement']);log=root/'dispatch_logs'/(path.stem+'.log')
        branch=None;matched_source=None
        if log.exists() and not (placement/'launch.json').exists():
            for version in versions:
                branch=refusal_branch(version,log.read_text(),RESOURCE_MESSAGES,executed_path=source)
                if branch is not None:matched_source=version;break
        if branch is None:failures+=1
        records.append(dict(receipt=str(path),receipt_sha256=sha256_file(path),
            classification='resource_contention_no_scientific_attempt' if branch else 'environment_or_unclassified_prelaunch_failure',
            log_sha256=sha256_file(log) if log.exists() else None,branch=branch,
            launcher_source=str(matched_source) if matched_source else None,
            launcher_source_sha256=hashes[str(matched_source)] if matched_source else None,
            owner=receipt['owner']))
    if (any(sha256_file(p)!=digest for p,digest in hashes.items()) or (root/'case/status.json').exists()
            or (root/'case/residency_guard/owner.json').exists()):
        raise ValueError('Launcher source or attempt state changed during accounting')
    return dict(submission_receipts=len(records),prelaunch_failures=failures,
        resource_refusals=len(records)-failures,records=records,
        launcher_sources_sha256=hashes,
        claims_modified=False,scientific_attempts_restarted=False)


def nonresource_prelaunch_failures(namespace,alive,source_versions=()):
    return prelaunch_accounting(namespace,alive,source_versions)['prelaunch_failures']
