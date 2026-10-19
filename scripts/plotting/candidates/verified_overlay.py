"""Read admitted expansion results for overlays on the existing figures."""
import hashlib
import json
import math
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'src'))
from redox.validation.provenance import (referenced_files, calculation_gid,
                                         scoped_spin_summary, input_changes_quarantined)
DEFAULT_FEED=ROOT/'results/candidates/current.json'
LEGACY_FEED=ROOT/'results/today_2026-10-11_verified_collection/latest.json'


def finite(x):return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)


def input_inventory_changed(summary):
    if 'spin_input_gids' not in summary or 'spin_input_paths' not in summary:return True
    current=set()
    for gid in summary['spin_input_gids']:
        current.update((ROOT/'calcs/uma'/gid).glob('*/result.json'))
        current.update((ROOT/'calcs/dft'/gid).glob('*/spin_validation/*.json'))
        for helper in (ROOT/'calcs/spincheck'/gid).glob('*_m*'):
            for pattern in ('result.json','opt.xyz','thermal_policy.json','sp/*.json','native_thermal/*.json',
                            'native_thermal/*/hessian_cartesian.npy','native_thermal/*/hessian_diagnostics.json'):
                current.update(helper.glob(pattern))
    current={p for p in current if p.is_file()}
    for p in list(current):
        if p.suffix=='.json':
            try:current.update(referenced_files(json.loads(p.read_text()),ROOT))
            except (OSError,ValueError):return True
    if {str(p) for p in current}!=set(summary['spin_input_paths']):return True
    source_files=set((ROOT/'src/redox').rglob('*.py'))|set((ROOT/'config').glob('*.py'))|set((ROOT/'config').glob('*.json'))
    expected={p for p in summary.get('source_sha256',{}) if Path(p).is_relative_to(ROOT/'src/redox') or Path(p).parent==ROOT/'config'}
    if {str(p) for p in source_files if p.is_file()}!=expected:return True
    for p in (ROOT/'config').glob('native*.json'):
        try:dependencies=referenced_files(json.loads(p.read_text()),ROOT)
        except (OSError,ValueError):return True
        if any(str(q) not in summary.get('source_sha256',{}) for q in dependencies):return True
    return False


def snapshot_inputs_current(summary):
    if summary.get('arithmetic_discrepancies')!=0 or not input_changes_quarantined(summary) or input_inventory_changed(summary):return False
    inputs={**summary.get('source_sha256',{}),**summary.get('admitted_input_sha256',{})}
    try:return bool(inputs) and all(Path(p).exists() and hashlib.sha256(Path(p).read_bytes()).hexdigest()==digest for p,digest in inputs.items())
    except OSError:return False


def load(feed=None):
    path=Path(feed) if feed else DEFAULT_FEED if DEFAULT_FEED.exists() else LEGACY_FEED
    if not path.exists():return None,[],[]
    payload=json.loads(path.read_text());summary=payload['summary']
    if summary.get('arithmetic_discrepancies')!=0 or not input_changes_quarantined(summary):
        return payload,[],[]
    rows=payload['original39']+payload.get('additional_candidates',payload.get('extra_NDI',[]))
    quarantined=set(summary.get('quarantined_input_gids',[]))
    display_gids={(r.get('registry_id') or r['id']) for r in rows
                  if r.get('solution_potentials_status')=='verified native calculation'
                  and (r.get('registry_id') or r['id']) not in quarantined}
    summary=scoped_spin_summary(summary,display_gids|{'ferrocene'},ROOT)
    changed=set()
    if input_inventory_changed(summary):changed.add('*')
    inputs={**summary.get('source_sha256',{}),**summary.get('admitted_input_sha256',{})}
    for name,digest in inputs.items():
        p=Path(name)
        try:valid=p.exists() and hashlib.sha256(p.read_bytes()).hexdigest()==digest
        except OSError:valid=False
        if not valid:
            changed.add(calculation_gid(p,ROOT) or '*')
    if '*' in changed:return payload,[],[]
    admitted=[r for r in rows if r.get('solution_potentials_status')=='verified native calculation'
              and (r.get('registry_id') or r['id']) not in changed and not ({'ferrocene','*'}&changed)]
    admitted=[r for r in admitted if (r.get('registry_id') or r['id']) not in quarantined]
    p=Path(payload['snapshot_path'])/'couples.csv'
    import csv
    couples=list(csv.DictReader(p.open())) if p.exists() else []
    return payload,rows,admitted


def short(gid):
    p=LEGACY_FEED.parent/'candidate_display_names.json'
    labels=json.loads(p.read_text()) if p.exists() else {}
    return labels.get(gid,{}).get('short_name',gid)


def family(name):
    if 'phenothiazine' in name:return 'phenothiazine (p-type)'
    if 'imide' in name:return 'imide (n-type)'
    if 'quinone' in name:return 'quinone (n-type)'
    if 'viologen' in name or 'pyridine-multi' in name:return 'pyridine-multi-e'
    return name
