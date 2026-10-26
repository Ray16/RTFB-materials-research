import importlib.util
from pathlib import Path

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('native_resource_witness_test',ROOT/'scripts/validation/prepare_native_resource_witness.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_only_temporary_response_arrays_leave_successor_dependency_closure(tmp_path):
    raw=tmp_path/'raw'
    names=['execution/rks_make_h1_0.npy','execution/uks_make_h1_1.npy','execution/rks_make_h1_state_identity.json',
           'hessian.npy','gradient.npy','geometry.xyz','experimental_partition.json','SCF_evaluations/0000/evaluation.json',
           'execution/response_solver.json','execution/rks_partial_hess_elec_0.npy']
    original={str(raw/name):name for name in names}
    required,temporary=module.project_sources(original,raw)
    assert set(temporary)=={str(raw/'execution/rks_make_h1_0.npy'),str(raw/'execution/uks_make_h1_1.npy')}
    assert required[str(raw/'hessian.npy')]=='hessian.npy'
    assert set(required)|set(temporary)==set(original)
    assert all(value==original[name] for name,value in required.items())


def test_similarly_named_arrays_elsewhere_or_unknown_execution_files_are_retained(tmp_path):
    raw=tmp_path/'raw';names=[raw/'rks_make_h1_0.npy',tmp_path/'other/execution/rks_make_h1_0.npy',
        raw/'execution/rks_make_h1_unrecognized.npy',raw/'execution/hessian.npy']
    original={str(path):'hash' for path in names}
    required,temporary=module.project_sources(original,raw)
    assert required==original and not temporary
