from pathlib import Path
import json,sys,os
from importlib import metadata
from redox.workflows.resources import verify_gate_owner
verify_gate_owner()
import ase
import torch
from fairchem.core import pretrained_mlip,FAIRChemCalculator
from fairchem.core.units.mlip_unit.api.inference import guess_inference_settings
assert guess_inference_settings('nitro').nitro is True
assert 'uma-s-1p2p1' in pretrained_mlip.available_models
report=dict(status='public_API_imports_registered_checkpoint_and_Nitro_preset_verified',python=sys.executable,
 base_prefix=sys.base_prefix,LD_LIBRARY_PATH=os.environ.get('LD_LIBRARY_PATH'),torch_version=torch.__version__,
 compiled_CUDA_version=torch.version.cuda,fairchem_version=metadata.version('fairchem-core'),ASE_version=ase.__version__,
 actual_Nitro_inference_executed=False,scope='Import and preset check only; no prediction or numerical GPU qualification')
p=Path('/nfs/lambda_stor_01/homes/rzhu/0_redox/results/validation/uma_nitro_2026-10-23/import_probe.json');p.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
