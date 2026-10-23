import ast
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def test_workspace_control_retains_original_models_and_threshold():
    module=load('scripts/validation/qualify_smd_workspace.py','workspace_qualifier')
    text=(ROOT/'scripts/validation/qualify_smd_memory_execution.py').read_text();node=next(n for n in ast.parse(text).body if isinstance(n,ast.FunctionDef) and n.name=='main')
    before='\n'.join(text.splitlines()[node.lineno-1:node.end_lineno]);after=module.adapted_body(before)
    assert 'value > 1e-6' in after and 'basis="def2-svp"' in after and '"r2scan"' in after
    assert 'budget_GiB=args.workspace_GiB' in after


def test_changed_control_template_is_rejected():
    module=load('scripts/validation/qualify_smd_workspace.py','workspace_qualifier')
    with pytest.raises(ValueError,match='template changed'):module.adapted_body('def main(): pass')


def test_SCF_checkpoint_identity_distinguishes_geometry_state_and_medium():
    module=load('scripts/pipeline/smd_workspace_execution.py','workspace_identity')
    molecule=SimpleNamespace(elements=['H','H'],atom_coords=lambda:np.array([[0,0,0],[0,0,1.5]]),basis='def2-tzvpd',charge=0,spin=0)
    mf=SimpleNamespace(mol=molecule,xc='wb97m-v',with_solvent=SimpleNamespace(solvent='acetonitrile'),init_guess='minao')
    a=module.molecular_identity(mf);molecule.charge=1;molecule.spin=1;b=module.molecular_identity(mf)
    assert a['charge']==0 and a['mult']==1 and b['charge']==1 and b['mult']==2
    mf.with_solvent=None;assert module.molecular_identity(mf)['solvent'] is None
