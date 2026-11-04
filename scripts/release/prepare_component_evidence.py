#!/usr/bin/env python3
"""Prepare a local review archive of the closed implementation-control evidence."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
from redox.core.common import ROOT,HARTREE_EV,EV_KJ
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def prepare(output, *, viologen_manifest, base_manifest, reviews):
    if output.exists() or output.is_relative_to(ROOT):raise FileExistsError('Fresh durable package outside checkout required')
    summary=json.loads((reviews/'summary.json').read_text())
    if summary['status']!='all_exact_component_actual_checkpoint_frames_verified':raise ValueError('Complete actual-frame review required')
    panels={'viologen':viologen_manifest.parent,'base':base_manifest.parent}
    ancestor=output.parent
    while not ancestor.exists():ancestor=ancestor.parent
    if shutil.disk_usage(ancestor).free<34*2**30:raise RuntimeError('Review archive requires the shared 32-GiB margin plus bounded copy headroom')
    from pyscf.lib import param
    output.mkdir(parents=True);(output/'code').mkdir();(output/'panels').mkdir()
    for source,name in [(ROOT/'scripts/release/replay_component_controls.py','replay.py'),(ROOT/'src/redox/validation/scf_frame.py','scf_frame.py'),(ROOT/'src/redox/validation/component_sensitivity.py','component_sensitivity.py')]:shutil.copyfile(source,output/'code'/name)
    shutil.copyfile(ROOT/'LICENSE',output/'CODE_LICENSE')
    shutil.copytree(ROOT/'tests/fixtures/licenses',output/'third_party_licenses')
    cases={};by_namespace={};panel_map={};expected={};protocol=None
    for name,folder in panels.items():
        manifest=json.loads((folder/'manifest.json').read_text());mapping={}
        for row in manifest['states']:
            review_path=reviews/name/(row['key']+'.json')
            matching=[h for p,h in summary['panels'][name]['source_sha256'].items() if Path(p).name==review_path.name]
            if len(matching)!=1 or sha256_file(review_path)!=matching[0]:
                raise ValueError('Independent component review changed')
            review=json.loads(review_path.read_text())
            if not review['actual_checkpoint_frames_verified'] or any(sha256_file(p)!=h for p,h in review['source_sha256'].items()):raise ValueError('Independent backing changed')
            ns=str(Path(row['namespace']).resolve())
            plan=json.loads((Path(ns).parent/'case/plan.json').read_text())
            if protocol is None:protocol=plan['electronic_protocol']
            elif plan['electronic_protocol']!=protocol:raise ValueError('Mixed electronic protocols in a review collection')
            if ns not in by_namespace:
                key='case_'+str(len(cases)).zfill(3);by_namespace[ns]=key;target=output/'cases'/key;target.mkdir(parents=True)
                origin=Path(ns).parent
                for sub in ('sources','case'):
                    if any(p.is_symlink() for p in (origin/sub).rglob('*')):raise ValueError('No source symlinks allowed')
                    shutil.copytree(origin/sub,target/sub)
                shutil.copyfile(row['namespace'],target/'namespace.json');shutil.copyfile(review_path,target/'independent_review.json')
                # Runtime records are descriptive provenance in a portable replay;
                # a historical PID cannot be checked for liveness after extraction.
                placements=[p for p in origin.glob('placement_*') if (p/'completion.json').exists() and json.loads((p/'completion.json').read_text()).get('returncode')==0]
                if len(placements)!=1:raise ValueError('One accepted closed placement required')
                (target/'execution_receipts').mkdir()
                for leaf in ('launch.json','completion.json'):shutil.copyfile(placements[0]/leaf,target/'execution_receipts'/leaf)
                cases[key]=dict(folder=str(target.relative_to(output)),original_namespace=ns,review_sha256=sha256_file(review_path))
            mapping[row['key']]=by_namespace[ns]
        package=folder/'inputs/package_manifest.json';dest=output/'panels'/(name+'_manifest.json');shutil.copyfile(package,dest)
        panel_map[name]=dict(package=str(dest.relative_to(output)),case_keys=mapping)
        report=(viologen_manifest if name=='viologen' else base_manifest).parent/'inputs/package_manifest.json'
        from redox.validation.component_sensitivity import reconstruct
        energies={}
        for row in manifest['states']:
            case=Path(row['namespace']).parent/'case';plan=json.loads((case/'plan.json').read_text())
            energies[row['key']]=dict(stored_Ha=plan['reference_energy_Ha'],installed_Ha=json.loads((case/'installed/result.json').read_text())['e_smd_Ha'],corrected_Ha=json.loads((case/'corrected/result.json').read_text())['e_smd_Ha'],charge=row['source_state']['charge'],mult=row['source_state']['mult'])
        expected[name]=reconstruct(json.loads(report.read_text()),energies)
    atomic_write_json(output/'expected_replay.json',expected)
    (output/'README.md').write_text('''# Closed electronic implementation controls — local review draft

This package contains 29 distinct completed execution namespaces for ten-state
viologen and twenty-one-input base-reaction panels; shared reference states
are mapped explicitly and the counts are not additive unique chemistry counts.
It includes actual SCF HDF5 checkpoints, geometry/method/state metadata, source
snapshots, installed/corrected scalar results and independent review receipts.
Failed/interrupted attempts remain at their original shared-storage locations
and are not promoted into this accepted-state archive.

Run `python code/replay.py PACKAGE` with NumPy, h5py and CPU PySCF installed.
The replay hashes every package file, checks the actual molecular checkpoint
frames and selected SCF energies, verifies atom/charge balance, and reconstructs
only electronic contributions to the declared potential/reaction shifts.
It never opens original checkout/calculation paths or reruns SCF. Producer spin
diagnostics are retained but not independently recalculated by this replay.
Nuclear terms, complete reaction G, material qualification and experimental
accuracy are not supplied. Historical process/resource receipts are preserved
as review provenance; extracted files cannot prove current process liveness.

Code retains MIT, with upstream files under the supplied third-party licenses.
Computed-data reuse terms, author metadata and public deposition remain pending.
This is a local review draft, not a DOI/public release or submission certificate.
''')
    metadata=dict(status='closed_component_evidence_local_review_package',cases=cases,panels=panel_map,
        hartree_eV=HARTREE_EV,eV_kJ_mol=EV_KJ,bohr_A=param.BOHR,
        electronic_protocol=protocol,
        files_sha256={str(p.relative_to(output)):sha256_file(p) for p in output.rglob('*') if p.is_file()},
        new_G_qualified=False,public_deposition_complete=False,data_license=None,author_metadata=None)
    atomic_write_json(output/'manifest.json',metadata);return metadata


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    for name in ('viologen-manifest','base-manifest','reviews'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    import os
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    r=prepare(a.out.resolve(),viologen_manifest=a.viologen_manifest.resolve(),base_manifest=a.base_manifest.resolve(),reviews=a.reviews.resolve());print(json.dumps(dict(cases=len(r['cases']),files=len(r['files_sha256']),public_deposition_complete=False)))
