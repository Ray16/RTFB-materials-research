import hashlib,json
import numpy as np
import pytest
from redox.validation.solution_relaxation import scf_cache_matches,verified_density


@pytest.mark.parametrize('schema',['legacy','paired'])
def test_both_cache_schemas_bind_geometry_species_and_frozen_reference(schema):
    state={'charge':1,'mult':2};method={'solvent_model':'native_SMD_frozen_slow'}
    if schema=='legacy':key=dict(identity=state,geometry={'geom_sha1':'short'},method=method,reference={'charge':0},refdm_sha256='density',source={'src/redox/qm/nonequilibrium_smd.py':'code'})
    else:key=dict(state=state,geometry_sha256='long',identity={'method':method,'source_sha256':{'/repo/src/redox/qm/nonequilibrium_smd.py':'code'}},reference={'charge':0},reference_density_sha256='density')
    record={'converged':True,'key':key}
    kwargs=dict(state=state,geometry_sha256='long',geometry_sha1='short',method=method,reference_density_sha256='density',numerical_sources={'src/redox/qm/nonequilibrium_smd.py':'code'})
    assert scf_cache_matches(record,**kwargs)
    assert not scf_cache_matches(record,**dict(kwargs,reference_density_sha256='other'))
    assert not scf_cache_matches(record,**dict(kwargs,numerical_sources={'src/redox/qm/nonequilibrium_smd.py':'changed'}))


@pytest.fixture
def density_case(tmp_path):
    p=tmp_path/'scf.json';dm=np.diag([2.,0.]);np.savez_compressed(p.with_suffix('.npz'),dm=dm)
    sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
    record=dict(e_Ha=-1.,e_eV=-27.211386245988,converged=True,electron_trace=2.)
    p.write_text(json.dumps(record));evidence={str(f):sha(f) for f in [p,p.with_suffix('.npz')]}
    return p,record,evidence,dm


def test_legacy_density_needs_review_of_both_json_and_array(density_case):
    p,r,evidence,dm=density_case;kwargs=dict(overlap=np.eye(2),electrons=2,mult=1)
    with pytest.raises(ValueError,match='Legacy density'):verified_density(r,p,evidence={},**kwargs)
    assert np.array_equal(verified_density(r,p,evidence=evidence,**kwargs),dm)


@pytest.mark.parametrize('change',['bytes','electron_count','units','spin','convergence'])
def test_density_cache_corruption_or_quantum_identity_change_is_refused(density_case,change):
    p,r,evidence,dm=density_case;kwargs=dict(overlap=np.eye(2),electrons=2,mult=1)
    if change=='bytes':np.savez_compressed(p.with_suffix('.npz'),dm=dm+.01)
    elif change=='electron_count':kwargs['electrons']=1
    elif change in ('units','convergence'):
        r['e_eV']=-1. if change=='units' else r['e_eV']
        if change=='convergence':r['converged']=False
        p.write_text(json.dumps(r));evidence[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    else:kwargs['mult']=3
    with pytest.raises(ValueError):verified_density(r,p,evidence=evidence,**kwargs)
