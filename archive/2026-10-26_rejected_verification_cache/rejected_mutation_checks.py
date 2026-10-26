import hashlib
import importlib.util
import os
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('verification_epoch_test',ROOT/'scripts/fleet/verification_epoch.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def checker(calls):
    def verify(mapping):
        for path,expected in mapping.items():
            calls.append(str(path))
            if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=expected:raise RuntimeError('Changed bytes')
    return verify


def test_large_object_is_hashed_once_per_pass_and_small_code_always(tmp_path):
    binary=tmp_path/'hessian.npy';binary.write_bytes(b'original');code=tmp_path/'code.py';code.write_bytes(b'original')
    expected=hashlib.sha256(b'original').hexdigest();calls=[];epoch=module.VerificationEpoch(checker(calls),minimum_bytes=1)
    epoch.begin();epoch.verify({binary:expected,code:expected});epoch.verify({binary:expected,code:expected})
    assert calls.count(str(binary))==1 and calls.count(str(code))==2
    epoch.begin();epoch.verify({binary:expected});assert calls.count(str(binary))==2


def test_in_place_change_with_restored_mtime_cannot_reuse_verification(tmp_path):
    path=tmp_path/'hessian.npy';path.write_bytes(b'original');expected=hashlib.sha256(path.read_bytes()).hexdigest()
    epoch=module.VerificationEpoch(checker([]),minimum_bytes=1);epoch.begin();epoch.verify({path:expected})
    old=path.stat();path.write_bytes(b'changed!');os.utime(path,ns=(old.st_atime_ns,old.st_mtime_ns))
    with pytest.raises(RuntimeError):epoch.verify({path:expected})


def test_changed_expectation_and_atomic_replacement_are_rechecked(tmp_path):
    path=tmp_path/'hessian.npy';path.write_bytes(b'original');expected=hashlib.sha256(path.read_bytes()).hexdigest()
    epoch=module.VerificationEpoch(checker([]),minimum_bytes=1);epoch.begin();epoch.verify({path:expected})
    with pytest.raises(RuntimeError):epoch.verify({path:'different expected hash'})
    old=path.stat();other=tmp_path/'replacement.npy';other.write_bytes(b'changed!');os.utime(other,ns=(old.st_atime_ns,old.st_mtime_ns));other.replace(path)
    with pytest.raises(RuntimeError):epoch.verify({path:expected})


def test_an_object_changed_during_hashing_is_never_cached(tmp_path):
    path=tmp_path/'hessian.npy';path.write_bytes(b'original');expected=hashlib.sha256(path.read_bytes()).hexdigest()
    def mutate(mapping):checker([])(mapping);path.write_bytes(b'changed!')
    epoch=module.VerificationEpoch(mutate,minimum_bytes=1);epoch.begin()
    with pytest.raises(RuntimeError,match='while hashing'):epoch.verify({path:expected})
    assert not epoch.cache
