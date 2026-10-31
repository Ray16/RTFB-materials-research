"""Append declared artifacts to an immutable manuscript evidence epoch.

This handles byte identity and ownership only. Callers must establish the
scientific qualification and interpretation of each added artifact separately.
"""
import hashlib
import json
from pathlib import Path
import shutil
import time

PRODUCER_BYTES=Path(__file__).read_bytes()


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def inside(root,relative):
    relative=Path(relative)
    if relative.is_absolute() or '..' in relative.parts:raise ValueError('Use a project-relative artifact path')
    path=(root/relative).resolve()
    if not path.is_relative_to(root):raise ValueError('Artifact escapes the project')
    return path


def append_epoch(root,name,extras,scope,source_commit,figure_extras=None):
    root=Path(root).resolve()
    if not name or Path(name).name!=name or name in ('.','..') or not scope or not source_commit:
        raise ValueError('An epoch name, declared scope and source revision are required')
    pointer_path=root/'manuscript/evidence/current.json';pointer_bytes=pointer_path.read_bytes();pointer=json.loads(pointer_bytes)
    parent=inside(root,pointer['manifest']);parent_hash=sha(parent);manifest=json.loads(parent.read_text());old_figures=inside(root,pointer['figures'])
    evidence=root/'manuscript/evidence'/name;figures=root/'results/figures/manuscript'/name
    if evidence.exists() or figures.exists():raise FileExistsError('Preserve previous evidence epochs')
    mapping={}
    for relative,expected in manifest['snapshot_sha256'].items():
        source=inside(root,relative)
        if sha(source)!=expected:raise ValueError('Parent evidence changed: '+relative)
        if source.parent==parent.parent:target=evidence/source.name
        elif source.parent==old_figures:target=figures/source.name
        else:raise ValueError('Parent artifact is outside its declared flat evidence/figure directories')
        mapping[target]=source
    for key,relative in extras.items():
        if not key or Path(key).name!=key or key in ('.','..','manifest.json','CAPTURE_FAILED.txt'):raise ValueError('Added evidence names must be flat filenames')
        target=evidence/key
        if target in mapping:raise ValueError('An addition would replace inherited evidence: '+key)
        mapping[target]=inside(root,relative)
    for key,relative in (figure_extras or {}).items():
        if not key or Path(key).name!=key or key in ('.','..'):
            raise ValueError('Added figure names must be flat filenames')
        target=figures/key
        if target in mapping:raise ValueError('An addition would replace an inherited figure: '+key)
        mapping[target]=inside(root,relative)
    identities={str(p.relative_to(root)):sha(p) for p in set(mapping.values())}
    evidence.mkdir();figures.mkdir()
    try:
        for target,source in mapping.items():shutil.copyfile(source,target)
        if (pointer_path.read_bytes()!=pointer_bytes or sha(parent)!=parent_hash
                or any(sha(root/path)!=expected for path,expected in identities.items())):
            raise ValueError('Evidence or its parent changed during capture')
        result=dict(kind='immutable_draft_evidence_snapshot',captured_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                    source_commit=source_commit,parent_manifest=pointer['manifest'],parent_manifest_sha256=parent_hash,
                    source_sha256=identities,snapshot_sha256={str(p.relative_to(root)):sha(p) for p in mapping},
                    scope=scope,extension_producer_sha256=hashlib.sha256(PRODUCER_BYTES).hexdigest(),
                    scientific_qualification_inferred_from_capture=False)
        (evidence/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')
        return dict(manifest=str((evidence/'manifest.json').relative_to(root)),figures=str(figures.relative_to(root)))
    except BaseException as error:
        (evidence/'CAPTURE_FAILED.txt').write_text(repr(error)+'\n')
        raise
