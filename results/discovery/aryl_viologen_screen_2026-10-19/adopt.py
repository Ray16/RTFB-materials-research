"""Install verified new library/UMA trees without regenerating or replacing inputs."""
from contextlib import ExitStack
import fcntl
import json
from pathlib import Path
import shutil
import time

from redox.build.contracts import validate_group
from redox.build.manifest import locked_manifest, read_rows, write_rows
from redox.core.common import ROOT, load_config
from redox.core.io import atomic_write_json
from redox.qm import uma
from redox.validation.protonation import hydrogen_topology_status
from rdkit import Chem

base=Path(__file__).resolve().parent
receipt=base/'registration_prepared/adoption.json'
assert not receipt.exists(), 'Existing adoption requires verification rather than replacement'
review=json.loads((base/'uma_completion.json').read_text())
assert review['status']=='all_27_isolated_UMA_points_reviewed_terminal' and review['verified_states']==27
plan=json.loads((base/'uma_plan.json').read_text())
groups=json.loads((base/'registration_prepared/proposals.json').read_text())
registered={g['id']:g for g in load_config('custom_candidates').GROUPS}
assert all(registered[g['id']]==g for g in groups)
ids={g['id'] for g in groups}
rows=read_rows(base/'isolated_library/manifest.csv')
assert len(rows)==27 and {r['id'] for r in rows}==ids
for group in groups:
    declared=Chem.MolToSmiles(validate_group('Cc1ccc(C[*:1])cc1',group))
    group_rows=[r for r in rows if r['id']==group['id']]
    assert len(group_rows)==3 and all(r['smiles']==declared for r in group_rows)
    assert {(r['state'],int(r['charge']),int(r['mult_hint']),int(r['n_e'])) for r in group_rows}=={tuple(s) for s in group['states']}
for path,digest in plan['source_sha256'].items():assert uma.file_sha256(path)==digest
staging=base/'registration_prepared/install_staging';staging.mkdir(exist_ok=False)
installed=[];files={}
manifest=ROOT/'library/manifest.csv'
with ExitStack() as stack:
    for gid in sorted(ids):
        for suffix in ('','_native_thermal'):
            lock=stack.enter_context((ROOT/'logs/pipeline_locks'/f'{gid}{suffix}.lock').open('a+'))
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    with locked_manifest(manifest):
        original=manifest.read_bytes();old_rows=read_rows(manifest)
        assert not ids & {r['id'] for r in old_rows}
        (base/'registration_prepared/manifest_before.csv').write_bytes(original)
        for gid in sorted(ids):
            assert not (ROOT/'library'/gid).exists() and not (ROOT/'calcs/uma'/gid).exists()
            shutil.copytree(base/'isolated_library'/gid,staging/'library'/gid)
            shutil.copytree(base/'uma'/gid,staging/'uma'/gid)
        for task in plan['tasks']:
            row=task['row'];gid=row['id'];state=row['state']
            sd=staging/'uma'/gid/state;seeds=sorted((staging/'library'/gid/'conformers').glob('conf_*.xyz'))
            assert [uma.file_sha256(p) for p in seeds]==list(task['seed_sha256'].values())
            result=json.loads((sd/'result.json').read_text())
            assert not uma.cache_errors(result,row,sd/'relaxed.xyz',seeds,plan['model_identity'],.02)
            assert hydrogen_topology_status(row['smiles'],sd/'relaxed.xyz')=='ok'
        for kind,destination in [('library',ROOT/'library'),('uma',ROOT/'calcs/uma')]:
            for gid in sorted(ids):
                source=staging/kind/gid
                for p in source.rglob('*'):
                    if p.is_file():
                        relative=p.relative_to(source)
                        assert uma.file_sha256(p)==uma.file_sha256(base/('isolated_library' if kind=='library' else 'uma')/gid/relative)
                        files[str(destination/gid/relative)]=uma.file_sha256(p)
        try:
            for kind,destination in [('library',ROOT/'library'),('uma',ROOT/'calcs/uma')]:
                for gid in sorted(ids):
                    source=staging/kind/gid;target=destination/gid
                    source.replace(target);installed.append((target,source))
            assert manifest.read_bytes()==original
            write_rows(manifest,old_rows+rows)
        except BaseException:
            for target,source in reversed(installed):target.replace(source)
            raise
        assert read_rows(manifest)[:len(old_rows)]==old_rows
        assert all(uma.file_sha256(p)==h for p,h in files.items())
atomic_write_json(receipt,dict(status='nine_registered_graphs_27_warm_starts_installed_byte_identically',
    registered_ids=sorted(ids),library_rows_appended=27,previous_manifest_rows=len(old_rows),
    current_manifest_rows=len(read_rows(manifest)),canonical_sha256=files,
    input_regeneration=False,previous_rows_preserved=True,source_records_preserved=True,
    solution_properties_computed=0,checked=time.time()))
print('Installed nine new library graphs and 27 UMA warm starts without regeneration')
