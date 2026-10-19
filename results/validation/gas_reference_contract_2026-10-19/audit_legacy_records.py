"""Inventory existing reference evidence without changing or certifying legacy values."""
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
import time
from redox.core.common import ROOT
from redox.core.io import atomic_write_json

base=Path(__file__).resolve().parent
cohorts={
 'D3TaLES_matched_gas':list((ROOT/'results/d3tales_reorg_validation_d3level/calc').glob('*.json')),
 'literature_inner_acene':list((ROOT/'results/reorg_anchors/inner_acene').glob('*/calc/*.json'))+list((ROOT/'results/reorg_anchors/inner_acene/calc').glob('*.json')),
 'literature_other':list((ROOT/'results/reorg_anchors/electron').glob('*/calc/*.json')),
}
rows=[]
for cohort,files in cohorts.items():
 for path in sorted(set(files)):
  blob=path.read_bytes();digest=hashlib.sha256(blob).hexdigest()
  try:
   r=json.loads(blob)
   if not isinstance(r,dict):raise ValueError('record is not an object')
  except (ValueError,TypeError):
   rows.append(dict(cohort=cohort,path=str(path.relative_to(ROOT)),sha256=digest,classification='malformed'));continue
  values={k:r.get(k) for k in ('our_hole','our_electron','our_lambda') if r.get(k) is not None}
  finite=all(not isinstance(v,bool) and isinstance(v,(int,float)) and math.isfinite(v) for v in values.values())
  status=('nonfinite_calculated_value' if not finite else 'solver_verification_present_requires_source_check'
          if r.get('solver_verification') else 'legacy_solver_convergence_not_recorded')
  rows.append(dict(cohort=cohort,path=str(path.relative_to(ROOT)),sha256=digest,status=r.get('status'),
                   finite_calculated_values=finite,calculated_fields=list(values),classification=status,
                   recorded_converged=r.get('converged'),protocol=r.get('protocol')))
assert all(hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()==r['sha256'] for r in rows)
summary=dict(status='read_only_reference_evidence_inventory',observed=time.time(),records=len(rows),
 by_cohort=dict(Counter(r['cohort'] for r in rows)),by_classification=dict(Counter(r['classification'] for r in rows)),
 inputs_current_during_inventory=True,existing_records_changed=False,
 scientifically_incorrect_values_proven=False,legacy_values_retroactively_verified=False,
 scope='Only these exact saved gas-reference cohorts; not production candidate states or all on-disk calculations')
atomic_write_json(base/'legacy_inventory.json',rows);atomic_write_json(base/'legacy_inventory_summary.json',summary)
print(json.dumps(summary,indent=2))
