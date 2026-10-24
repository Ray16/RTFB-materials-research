"""Preserve the native pilot after proving a mandatory allocation exceeds VRAM."""
from pathlib import Path
import importlib.util,json,socket,time
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
C=ROOT/'results/polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_3_native_first_response'
worker=json.loads((C/'thermal/workers/lambda2_gpu0_pilot.json').read_text())
assert worker['host']==socket.gethostname().split('.')[0]
plan=json.loads((C/'plan.json').read_text());raw=Path(plan['raw_root'])
records=[]
for p in raw.rglob('build_execution/status.json'):
 for r in json.loads(p.read_text()).get('builds',[]):
  if r['atoms']==115:records.append((p,r))
assert records
n=max(r['surface_points'] for p,r in records);required=2*n*115*115*3*3*8
assert required>32768*2**20
source=ROOT/'scripts/fleet/watch_gpu_exclusivity.py';spec=importlib.util.spec_from_file_location('native_bound_owned_stop',source);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
proof=dict(status='native_device_allocation_lower_bound_exceeds_device',checked=time.time(),host=worker['host'],pid=worker['pid'],process_start=worker['process_start'],n_surface_points=n,atoms=115,two_required_dense_tensor_bytes=required,two_required_dense_tensor_GiB=required/2**30,physical_device_MiB=32768,source_function='gpu4pyscf.solvent.hessian.pcm.get_d2F_d2A',reason='Two explicit float64 device zeros([ngrids,natom,natom,3,3]) allocations exceed the device even before other storage',baseline_geometry_unchanged=True,scientific_QC_failure_inferred=False,G_qualified=False)
atomic_write_json(ROOT/'results/validation/pcm_first_response_workspace_v2_2026-10-23/native_allocation_bound.json',proof)
proof.update(m.stop_owned(worker['pid'],worker['process_start'],C))
atomic_write_json(ROOT/'results/validation/pcm_first_response_workspace_v2_2026-10-23/native_allocation_bound.json',proof)
print(json.dumps(proof))
