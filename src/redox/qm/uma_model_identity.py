"""Local UMA fingerprints for wheel and editable FAIRChem installations.

Registry discovery never imports FAIRChem or initializes a model/GPU. Loaded
module metadata and a local editable-project URL are alternative location hints.
"""
from importlib import metadata
import hashlib
import json
import os
from pathlib import Path
import sys
from urllib.parse import unquote,urlparse


def digest(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):value.update(block)
    return value.hexdigest()


def registry_file(distribution):
    candidates=[Path(distribution.locate_file('fairchem/core/calculate/pretrained_models.json'))]
    loaded=sys.modules.get('fairchem.core.calculate')
    if loaded is not None and getattr(loaded,'__file__',None):
        candidates.append(Path(loaded.__file__).parent/'pretrained_models.json')
    raw=distribution.read_text('direct_url.json')
    if raw:
        direct=json.loads(raw);url=urlparse(direct.get('url',''))
        if url.scheme=='file' and direct.get('dir_info',{}).get('editable'):
            project=Path(unquote(url.path))
            candidates.extend([project/'src/fairchem/core/calculate/pretrained_models.json',
                project/'fairchem/core/calculate/pretrained_models.json'])
    existing=sorted({p.resolve() for p in candidates if p.is_file()})
    if not existing:raise RuntimeError('No local FAIRChem registry found; no model was initialized')
    if len({digest(p) for p in existing})!=1:raise RuntimeError('Conflicting FAIRChem registry locations')
    return existing[0]


def model_identity(model):
    from huggingface_hub import hf_hub_download
    distribution=metadata.distribution('fairchem-core');registry=registry_file(distribution)
    entries=json.loads(registry.read_text());entry=entries.get(model)
    if entry is None and model=='uma-s-1p2p1' and 'uma-s-1p2' in entries:
        entry=dict(entries['uma-s-1p2'],filename='uma-s-1p2p1.pt')
    if entry is None:raise RuntimeError('Requested model is absent from the local registry')
    cache=os.environ.get('FAIRCHEM_CACHE_DIR',str(Path.home()/'.cache/fairchem'))
    def local(record):
        return Path(hf_hub_download(repo_id=entry['repo_id'],filename=record['filename'],
            subfolder=record.get('subfolder'),revision=entry.get('revision'),cache_dir=cache,local_files_only=True))
    checkpoint=local(entry)
    references={k:digest(local(entry[k])) for k in ('atom_refs','form_elem_refs') if entry.get(k)}
    return dict(model=model,checkpoint_sha256=digest(checkpoint),reference_sha256=references,
        loader_version=distribution.version,registry_entry=entry,task_name='omol',
        registry_path=str(registry),registry_sha256=digest(registry))
