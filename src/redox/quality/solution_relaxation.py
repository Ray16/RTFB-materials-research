"""Read-only evidence checks for native-SMD relaxation implementation controls."""
import hashlib
import json
import math
from pathlib import Path


def verified_implementation_controls(evidence,source_module,*,backend):
    """Require bound CPU assembly proof and the requested backend's SCF control."""
    if backend not in ('cpu','gpu') or not isinstance(evidence,dict) or not evidence:
        raise ValueError('Explicit backend and implementation-control hashes required')
    module_hash=hashlib.sha256(Path(source_module).read_bytes()).hexdigest()
    found=set();bound={}
    for name,expected in evidence.items():
        path=Path(name)
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError('Implementation control changed')
        record=json.loads(path.read_text())
        if record.get('status')!='passed':raise ValueError('Implementation control did not pass')
        residual=record.get('residual_eV')
        if type(residual) not in (int,float) or not math.isfinite(residual) or abs(residual)>1e-5:
            raise ValueError('Implementation residual is unqualified')
        reference=record.get('hcore_reference_Ha',record.get('independent_hcore_Ha'))
        wrapper=record.get('wrapper_Ha')
        if any(type(x) not in (int,float) or not math.isfinite(x) for x in [reference,wrapper]):
            raise ValueError('Independent assembly energies missing')
        if abs((wrapper-reference)*27.211386245988-residual)>1e-8:
            raise ValueError('Implementation residual arithmetic differs')
        if record.get('module_sha256')==module_hash:
            delta=record.get('density_max_difference')
            if type(delta) not in (int,float) or not math.isfinite(delta) or not 0<=delta<=1e-6:
                raise ValueError('Independent density replay missing')
            found.add('cpu')
        if record.get('backend')=='gpu':found.add('gpu')
        bound[str(path)]=expected
    if 'cpu' not in found or backend not in found:
        raise ValueError('Current module and backend control are not both bound')
    return dict(implementation_reference_verified=True,backend=backend,
                module_sha256=module_hash,evidence_sha256=bound,
                individual_calculation_independently_replayed=False)


def scf_cache_matches(record,*,state,geometry_sha256,geometry_sha1,method,reference_density_sha256,
                      numerical_sources=None):
    """Accept both documented cache schemas, with exact species/geometry/response identity."""
    key=record.get('key') or {}
    if 'state' in key:
        species=key.get('state') or {};actual_method=(key.get('identity') or {}).get('method')
        same_geometry=key.get('geometry_sha256')==geometry_sha256
        reference_hash=key.get('reference_density_sha256')
    else:
        species=key.get('identity') or {};actual_method=key.get('method')
        same_geometry=(key.get('geometry') or {}).get('geom_sha1')==geometry_sha1
        reference_hash=key.get('refdm_sha256')
    if numerical_sources:
        saved=(key.get('identity') or {}).get('source_sha256',{}) if 'state' in key else key.get('source',{})
        normalized={}
        for path,digest in saved.items():
            relative=path
            for marker in ('/src/redox/','/scripts/'):
                if marker in path:relative=marker[1:]+path.split(marker,1)[1];break
            normalized[relative]=digest
        if any(normalized.get(path)!=digest for path,digest in numerical_sources.items()):return False
    return (record.get('converged') is True and actual_method==method and same_geometry
            and species.get('charge')==state['charge'] and species.get('mult')==state['mult']
            and (key.get('reference') is None)==(reference_density_sha256 is None)
            and reference_hash==reference_density_sha256)


def verified_density(record,path,*,evidence,overlap,electrons,mult):
    """Read density only when its bytes and quantum identity are independently bound."""
    import numpy as np
    path=Path(path);artifact=path.with_suffix('.npz')
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    if record.get('converged') is not True or json.loads(path.read_text())!=record:
        raise ValueError('Cached SCF record is unconverged or differs from its bound JSON')
    own_binding=record.get('density_artifact_sha256')
    if own_binding is not None:
        if sha(artifact)!=own_binding:raise ValueError('Cached density artifact changed')
    elif (not evidence or evidence.get(str(path))!=sha(path)
          or evidence.get(str(artifact))!=sha(artifact)):
        raise ValueError('Legacy density requires a source-bound independent review')
    for field in ('e_Ha','e_eV'):
        value=record.get(field)
        if type(value) not in (int,float) or not math.isfinite(value):raise ValueError('Invalid cached energy')
    if abs(record['e_eV']-record['e_Ha']*27.211386245988)>1e-8:raise ValueError('Cached energy units differ')
    with np.load(artifact,allow_pickle=False) as saved:dm=saved['dm'].copy()
    n=len(overlap);shape=(n,n) if mult==1 else (2,n,n)
    if dm.shape!=shape or not np.isfinite(dm).all() or np.iscomplexobj(dm):
        raise ValueError('Cached density shape/type differs from target species')
    if abs(dm-np.swapaxes(dm,-1,-2)).max()>1e-8:raise ValueError('Cached density is not symmetric')
    if record.get('density_matrix_sha256') not in (None,hashlib.sha256(dm.tobytes()).hexdigest()):
        raise ValueError('Cached density matrix changed')
    total=dm if dm.ndim==2 else dm.sum(axis=0)
    trace=float(np.einsum('ij,ji',total,overlap))
    if abs(trace-electrons)>1e-6:raise ValueError('Cached density electron count differs')
    if record.get('electron_trace') is not None and abs(record['electron_trace']-trace)>1e-6:
        raise ValueError('Cached electron trace differs')
    if mult>1:
        s2=record.get('s2');spin=(mult-1)/2
        if type(s2) not in (int,float) or not math.isfinite(s2) or abs(s2-spin*(spin+1))>.1:
            raise ValueError('Cached density spin diagnostics are unqualified')
    return dm
