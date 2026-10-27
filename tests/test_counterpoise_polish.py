import importlib.util
import json
import os
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('counterpoise_polish',ROOT/'scripts/pipeline/counterpoise_cpu_polish.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


def fixture(tmp_path):
    prior=tmp_path/'prior';(prior/'AB').mkdir(parents=True)
    (prior/'AB/SCF.h5').write_bytes(b'checkpoint')
    for name in ('grid_coordinates_AU.npy','grid_weights.npy'):(prior/name).write_bytes(name.encode())
    raw=dict(status='CPU_remaining_term_converged_independent_review_pending',role='AB',converged=True,
             E_Ha=-10.,owner=dict(host=os.uname().nodename,pid=999999999,process_start='none'),source_sha256={},
             checkpoint_sha256=m.api.digest(prior/'AB/SCF.h5'))
    (prior/'status.json').write_text(json.dumps(raw));(prior/'AB/energy.json').write_text(json.dumps(raw))
    (prior/'launch.json').write_text(json.dumps(dict(source_sha256={})))
    proof={str(p):m.api.digest(p) for p in (prior/'launch.json',prior/'status.json',prior/'AB/energy.json',prior/'AB/SCF.h5',prior/'grid_coordinates_AU.npy',prior/'grid_weights.npy')}
    review=tmp_path/'review.json';review.write_text(json.dumps(dict(status='CPU_remaining_term_checkpoint_and_sampled_resources_reviewed_gradient_pending',role='AB',source_sha256=proof)))
    gradient=tmp_path/'gradient.json';gradient.write_text(json.dumps(dict(status='remaining_stationarity_failed_preserved',role='AB',passed=False,orbital_gradient_threshold_au=1e-7,orbital_gradient_norm_au=1.13e-7,energy_error_Ha=0.,checkpoint_energy_error_Ha=0.,reference_energy_Ha=-10.,source_sha256=proof)))
    return prior,review,gradient


def test_only_failed_independent_orbital_gradient_can_seed_refinement(tmp_path):
    prior,review,gradient=fixture(tmp_path)
    raw,_,previous=m.prior_checkpoint(prior,review,gradient)
    assert raw['role']=='AB' and previous['checkpoint_sha256']==m.api.digest(prior/'AB/SCF.h5')
    assert m.TARGETS['SCF_conv_tol_grad_au'] < m.TARGETS['independent_acceptance_threshold_au']


@pytest.mark.parametrize('field,value',[('passed',True),('orbital_gradient_norm_au',1e-8),('energy_error_Ha',1e-3),('orbital_gradient_threshold_au',1e-6),('reference_energy_Ha',-11.)])
def test_success_changed_threshold_or_energy_failure_is_not_relabelled(tmp_path,field,value):
    prior,review,gradient=fixture(tmp_path);d=json.loads(gradient.read_text());d[field]=value;gradient.write_text(json.dumps(d))
    with pytest.raises((ValueError,RuntimeError)):m.prior_checkpoint(prior,review,gradient)


def test_checkpoint_bytes_must_match_exact_independent_reviews(tmp_path):
    prior,review,gradient=fixture(tmp_path);(prior/'AB/SCF.h5').write_bytes(b'changed')
    with pytest.raises((ValueError,RuntimeError)):m.prior_checkpoint(prior,review,gradient)
