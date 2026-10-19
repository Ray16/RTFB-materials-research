"""Review terminal isolated UMA records without promoting solution properties."""
import json
import time
from pathlib import Path

from redox.core.io import atomic_write_json
from redox.qm import uma
from redox.validation.protonation import hydrogen_topology_status

base=Path(__file__).resolve().parent
plan=json.loads((base/'uma_plan.json').read_text())
observations=[json.loads((base/f'{host}_observation.json').read_text()) for host in ('lambda5','lambda13')]
if not all(g['terminal'] for report in observations for g in report['gates']):
    raise SystemExit('Original gates remain live in last observation; refresh observations before terminal review')
assert sum(len(r['gates']) for r in observations)==27
for path,digest in plan['source_sha256'].items():assert uma.file_sha256(path)==digest
cases=[];failures=[]
for task in plan['tasks']:
    path=base/'tasks'/task['key']/'status.json'
    if not path.exists():raise SystemExit('Missing attempted point: '+task['key'])
    report=json.loads(path.read_text())
    if report['status']=='computing':raise SystemExit('Pending point: '+task['key'])
    assert report['task']==task
    if report['status']!='gas_UMA_state_and_hydrogen_topology_verified':
        failures.append(dict(key=task['key'],status=report['status'],error=report.get('error')))
        continue
    row=task['row'];directory=base/'uma'/row['id']/row['state']
    result=json.loads((directory/'result.json').read_text())
    assert not uma.cache_errors(result,row,directory/'relaxed.xyz',list(map(Path,task['seed_sha256'])),plan['model_identity'],.02)
    assert hydrogen_topology_status(row['smiles'],directory/'relaxed.xyz')=='ok'
    assert all(uma.file_sha256(p)==h for p,h in task['seed_sha256'].items())
    assert uma.file_sha256(directory/'result.json')==report['result_sha256']
    assert uma.file_sha256(directory/'relaxed.xyz')==report['geometry_sha256']
    peak=max(s['used_MiB'] for s in report['memory_samples'])
    assert peak==report['peak_device_used_MiB'] and peak*1024**2<.8*report['total_memory_bytes']
    assert report['peak_reserved_bytes']<.8*report['total_memory_bytes']
    cases.append(dict(key=task['key'],id=row['id'],state=row['state'],energy_eV=result['energy_eV'],
        chosen_mult=result['chosen_mult'],spin_gap_eV=result['spin_gap_eV'],
        DFT_spin_comparison_required=result['spin_gap_eV'] is not None and result['spin_gap_eV']<.5,
        fmax_eV_A=result['fmax'],peak_device_used_MiB=peak,
        result_sha256=report['result_sha256'],geometry_sha256=report['geometry_sha256']))
assert len(cases)+len(failures)==27
report=dict(status='all_27_isolated_UMA_points_reviewed_terminal' if not failures else 'terminal_UMA_screen_with_preserved_failures',
    cases=cases,failures=failures,verified_states=len(cases),failed_states=len(failures),
    complete_ladders=sum(sum(c['id']==gid for c in cases)==3 for gid in {t['row']['id'] for t in plan['tasks']}),
    source_sha256={str(base/'uma_plan.json'):uma.file_sha256(base/'uma_plan.json')},
    canonical_adoption=False,solution_properties_computed=0,source_measurements_inherited=False,
    scientific_accuracy_validated=False,checked=time.time())
atomic_write_json(base/'uma_completion.json',report)
print(json.dumps({k:v for k,v in report.items() if k not in ('cases','failures')}))
