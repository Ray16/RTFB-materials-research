#!/usr/bin/env python
"""Publish exact independently checked geometry reviews to declared consumers.

Source review bytes and physical provenance are preserved. Existing consumer
paths must already contain identical bytes; different files are never replaced.
"""
import argparse
import json
import os
from pathlib import Path
import tempfile
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def validated_copy(source,target,campaign):
    source,target,campaign=map(lambda p:Path(p).resolve(),(source,target,campaign))
    if not source.is_file() or not target.is_relative_to(ROOT/'results'):
        raise ValueError('Review source and a declared project-result destination required')
    data=source.read_bytes();review=json.loads(data)
    if (review.get('status')!='alternate_smooth_geometry_and_resource_pilot_independently_checked'
            or review.get('campaign')!=str(campaign) or review.get('primary_candidate_admission') is not False
            or review.get('active_electronic_energy_qualified') is not False
            or review.get('native_Hessian_or_G_qualified') is not False):
        raise ValueError('Independent geometry review has a different identity or qualification scope')
    api.verify_sources(review['source_sha256'])
    if source.read_bytes()!=data:raise ValueError('Review changed during source validation')
    target.parent.mkdir(parents=True,exist_ok=True)
    temporary=None
    try:
        with tempfile.NamedTemporaryFile(dir=target.parent,prefix='.geometry_review_',delete=False) as stream:
            temporary=Path(stream.name);stream.write(data);stream.flush();os.fsync(stream.fileno())
        try:os.link(temporary,target);kind='published'
        except FileExistsError:
            if target.read_bytes()!=data:raise ValueError('Preserve a different existing consumer review')
            kind='already_identical'
    finally:
        if temporary is not None:temporary.unlink(missing_ok=True)
    if target.read_bytes()!=data or source.read_bytes()!=data:raise ValueError('Review bytes changed while publishing')
    return dict(source=str(source),target=str(target),campaign=str(campaign),sha256=api.digest(target),status=kind)


def publish(batches,review_directory,out):
    if out.exists():raise FileExistsError('Preserve previous handoff audits')
    planned=[];sources={str(Path(__file__).resolve()):api.digest(Path(__file__).resolve())}
    for batch in batches:
        manifest=json.loads((batch/'batch.json').read_text());api.verify_sources(manifest['source_sha256'])
        if manifest['stage']!='electronic' or manifest['primary_candidate_admission']:
            raise ValueError('Only declared alternate electronic consumers are supported')
        sources[str(batch/'batch.json')]=api.digest(batch/'batch.json')
        for entry in manifest['entries']:
            if entry['existing_checked_reference']:continue
            campaign=Path(entry['geometry_campaign']);target=Path(entry['geometry_review']);source=review_directory/target.name
            if target.exists() or not source.exists():continue
            planned.append((source,target,campaign))
    # Validate every available addition before publishing the first one.
    for source,target,campaign in planned:
        review=json.loads(source.read_text())
        if review.get('campaign')!=str(campaign):raise ValueError('Review filename does not establish graph identity')
        api.verify_sources(review['source_sha256']);sources[str(source)]=api.digest(source)
    out.mkdir(parents=True,exist_ok=False)
    record=dict(status='checked_geometry_review_handoff_in_progress',started=time.time(),source_sha256=sources,
                entries=[],scientific_calculations_run=False,qualification_scope_changed=False)
    atomic_write_json(out/'summary.json',record)
    for source,target,campaign in planned:
        record['entries'].append(validated_copy(source,target,campaign));atomic_write_json(out/'summary.json',record)
    api.verify_sources(sources)
    record.update(status='checked_geometry_reviews_published_to_declared_consumers',finished=time.time(),published=len(record['entries']))
    atomic_write_json(out/'summary.json',record)
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--batch',type=Path,action='append',required=True)
    p.add_argument('--review-directory',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();r=publish([x.resolve() for x in a.batch],a.review_directory.resolve(),a.out.resolve())
    print(json.dumps(dict(status=r['status'],published=r['published'])))
