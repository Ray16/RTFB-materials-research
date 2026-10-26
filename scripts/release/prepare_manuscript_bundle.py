#!/usr/bin/env python
"""Consolidate immutable manuscript reviewer packages without publishing them."""
import argparse
import hashlib
import json
from pathlib import Path,PurePosixPath
import shutil
import subprocess
import time

ROOT=Path(__file__).resolve().parents[2]
CONFIG=ROOT/'config/manuscript_reviewer_bundle.json'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def inside(root,relative):
    path=PurePosixPath(relative)
    if path.is_absolute() or '..' in path.parts or not path.parts:raise ValueError('Unsafe relative source path')
    result=(root/str(path)).resolve()
    if not result.is_relative_to(root.resolve()):raise ValueError('Source escapes project root')
    return result


def prepare(root,config,out):
    root=root.resolve();out=out.resolve()
    if out.exists():raise FileExistsError('Never overwrite a reviewer bundle')
    config_path=config.resolve();spec=json.loads(config_path.read_text())
    if (spec.get('schema')!=1 or spec.get('public_deposition_authorized_by_manifest')
            or spec.get('data_license') is not None or spec.get('author_metadata') is not None):
        raise ValueError('Prepared-only bundle requires explicit unresolved publication metadata')
    paths=spec['source_paths'];packages=spec['packages']
    if not paths or not packages or len({r['id'] for r in packages})!=len(packages):raise ValueError('Unique declared inputs required')
    for path in paths:inside(root,path)
    inputs=[]
    for entry in packages:
        record=inside(root,entry['record']);data=json.loads(record.read_text())
        archive=inside(root,data['compressed_archive'])
        if sha(archive)!=data['archive_sha256']:raise ValueError('Source archive hash differs: '+entry['id'])
        inputs.append((entry,record,data,archive))
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()
    out.mkdir(parents=True);(out/'archives').mkdir();(out/'records').mkdir()
    try:
        project_tar=out/'project-source.tar'
        subprocess.run(['git','archive','--format=tar','--output',str(project_tar),head,'--',*paths],cwd=root,check=True,timeout=180)
        package_rows=[]
        for entry,record,data,archive in inputs:
            target=out/'archives'/archive.name;shutil.copyfile(archive,target)
            if sha(target)!=data['archive_sha256'] or sha(archive)!=data['archive_sha256']:
                raise ValueError('Archive changed while copying: '+entry['id'])
            record_target=out/'records'/(entry['id']+'.json');shutil.copyfile(record,record_target)
            if sha(record_target)!=sha(record):raise ValueError('Release record changed while copying')
            package_rows.append(dict(id=entry['id'],archive='archives/'+archive.name,archive_sha256=sha(target),
                source_record='records/'+entry['id']+'.json',scope=data.get('scope'),overlap_with_other_packages=True))
        if subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()!=head:
            raise ValueError('Git HEAD changed during source capture')
        (out/'README.md').write_text('''# Manuscript reviewer bundle (local preparation)\n\nThe source tar contains the committed code, environment locks, examples, manuscript and its frozen evidence epoch. The four nested archives contain previously checked raw scientific records and overlap in molecular identity; their state counts are not additive. Extract each archive and use its own README for its reproducible CPU replay.\n\nThe source tar is a committed snapshot, not a fresh GPU-environment certification. Cached replay does not repeat SCF, stationarity, spin or original experimental validation. Separate third-party source articles are identified by DOI and remain subject to their own licenses. Data-license and author metadata, public deposition and DOI are pending. This bundle is not a publication-readiness or material-eligibility certificate.\n''')
        files=sorted(p for p in out.rglob('*') if p.is_file())
        manifest=dict(status='local_manuscript_reviewer_bundle_prepared_not_public',created=time.time(),
            source_commit=head,source_paths=paths,packages=package_rows,
            files_sha256={str(p.relative_to(out)):sha(p) for p in files},config_sha256=sha(config_path),
            producer_sha256=sha(__file__),data_license=None,author_metadata=None,public_DOI_created=False,
            public_archive_deposited=False,third_party_source_policy=spec['third_party_source_policy'],
            full_gpu_reproduction_inferred=False,experimental_accuracy_inferred=False,conformer_convergence_inferred=False)
        (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');return manifest
    except BaseException as error:
        (out/'PREPARATION_FAILED.txt').write_text(repr(error)+'\n');raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    import os
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    x=prepare(ROOT,CONFIG,a.out);print(json.dumps({'status':x['status'],'source_commit':x['source_commit'],'packages':len(x['packages'])}))
