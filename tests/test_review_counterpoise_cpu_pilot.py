import importlib.util
from pathlib import Path

import numpy as np
import pytest
from pyscf.scf.chkfile import load_scf


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('cpu_pilot_review',ROOT/'scripts/validation/review_counterpoise_cpu_pilot.py')
reviewer=importlib.util.module_from_spec(spec);spec.loader.exec_module(reviewer)


def test_reviewed_water_checkpoint_has_consistent_closed_shell_orbitals():
    checkpoint=ROOT/'calcs/validation/2026-10-26_counterpoise/CPU_fixed_auxiliary_factory/A_own/SCF.h5'
    mol,state=load_scf(str(checkpoint))
    error=reviewer.verify_orbitals(state['mo_coeff'],state['mo_occ'],state['mo_energy'],
                                   mol.intor_symmetric('int1e_ovlp'),mol.nelectron)
    assert error < 1e-7


def test_rejects_wrong_occupation_and_nonorthogonal_orbitals():
    coeff=np.eye(3);occ=np.array([2.,0.,0.]);energy=np.array([-1.,1.,2.]);overlap=np.eye(3)
    assert reviewer.verify_orbitals(coeff,occ,energy,overlap,2)==0
    with pytest.raises(ValueError,match='closed-shell'):
        reviewer.verify_orbitals(coeff,np.array([1.,1.,0.]),energy,overlap,2)
    coeff[0,0]=2.
    with pytest.raises(ValueError,match='orthonormal'):
        reviewer.verify_orbitals(coeff,occ,energy,overlap,2)


def test_refuses_live_owner():
    import os
    from pathlib import Path
    start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]
    with pytest.raises(ValueError,match='still live'):
        reviewer.verify_owner_stopped(dict(host=os.uname().nodename,pid=os.getpid(),process_start=start))


def test_unfinished_pilot_is_never_promoted(tmp_path,monkeypatch):
    import json
    monkeypatch.setattr(reviewer,'ROOT',tmp_path)
    raw=tmp_path/'calcs'/'pilot';raw.mkdir(parents=True)
    (raw/'launch.json').write_text('{}')
    (raw/'status.json').write_text(json.dumps({'status':'CPU_pilot_SCF_running_unqualified'}))
    with pytest.raises(ValueError,match='Completed but unqualified'):
        reviewer.review(raw,tmp_path/'review.json')
    assert not (tmp_path/'review.json').exists()
