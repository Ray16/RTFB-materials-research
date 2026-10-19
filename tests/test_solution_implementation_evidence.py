import hashlib
import json
import pytest
from redox.validation.solution_relaxation import verified_implementation_controls


@pytest.fixture
def controls(tmp_path):
    module=tmp_path/'method.py';module.write_text('qualified equations')
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    cpu=tmp_path/'cpu.json';gpu=tmp_path/'gpu.json'
    cpu.write_text(json.dumps(dict(status='passed',module_sha256=digest(module),residual_eV=0.,
        wrapper_Ha=-10.,hcore_reference_Ha=-10.,density_max_difference=0.)))
    gpu.write_text(json.dumps(dict(status='passed',backend='gpu',residual_eV=0.,
        wrapper_Ha=-10.,independent_hcore_Ha=-10.)))
    return module,cpu,gpu,{str(p):digest(p) for p in [cpu,gpu]}


def test_bound_cpu_and_backend_controls_do_not_validate_an_individual_case(controls):
    module,_,_,evidence=controls
    result=verified_implementation_controls(evidence,module,backend='gpu')
    assert result['implementation_reference_verified']
    assert not result['individual_calculation_independently_replayed']


@pytest.mark.parametrize('change',['missing_gpu','changed_method','bad_residual','failed_control'])
def test_control_scope_and_source_changes_cannot_be_ignored(controls,change):
    module,cpu,gpu,evidence=controls
    if change=='missing_gpu':evidence.pop(str(gpu))
    elif change=='changed_method':module.write_text('different equations')
    else:
        record=json.loads(gpu.read_text())
        if change=='bad_residual':record['residual_eV']=.001
        else:record['status']='failed'
        gpu.write_text(json.dumps(record));evidence[str(gpu)]=hashlib.sha256(gpu.read_bytes()).hexdigest()
    with pytest.raises(ValueError):verified_implementation_controls(evidence,module,backend='gpu')
