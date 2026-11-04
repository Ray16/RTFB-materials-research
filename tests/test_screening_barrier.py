from pathlib import Path

import pytest
import numpy as np

from redox.qm.screening_barrier import standalone_source,reference_mask

SOURCE=(Path(__file__).parent/'fixtures/screen_index_v1_8_1.cu').read_text()
CONSTANTS=dict(ATM_SLOTS=6,PTR_COORD=1,BAS_SLOTS=8,ATOM_OF=0,ANG_OF=1,NPRIM_OF=2,PTR_EXP=5,PTR_COEFF=6)


def test_upstream_change_adds_only_the_reader_writer_barrier_to_numerical_body():
    original=standalone_source(SOURCE,CONSTANTS,corrected=False)
    changed=standalone_source(SOURCE,CONSTANTS,corrected=True)
    addition='        __syncthreads(); // upstream PR925: protect previous shared-coordinate readers\n'
    assert changed.replace(addition,'')==original
    assert changed.count('__syncthreads();')==2 and original.count('__syncthreads();')==1


@pytest.mark.parametrize('kind',['already_patched','wrong_loop','missing_ABI','invalid_ABI'])
def test_undefined_source_or_ABI_is_rejected(kind):
    source=SOURCE;constants=dict(CONSTANTS)
    if kind=='already_patched':source=source.replace('        __syncthreads();','        __syncthreads();\n        __syncthreads();')
    if kind=='wrong_loop':source=source.replace('grid0 += NG_PER_BLOCK','grid0 += 1')
    if kind=='missing_ABI':constants.pop('ATM_SLOTS')
    if kind=='invalid_ABI':constants['BAS_SLOTS']=False
    with pytest.raises(ValueError):standalone_source(source,constants)


def test_independent_bound_recognizes_near_s_shell_and_rejects_far_block():
    atm=np.zeros((1,6),dtype=np.int32);bas=np.zeros((1,8),dtype=np.int32)
    atm[0,1]=0;bas[0,2]=1;bas[0,5]=3;bas[0,6]=4
    env=np.array([0.,0.,0.,1.,1.]);coords=np.array([[.1,0.,0.],[.2,0.,0.],[100.,0.,0.],[101.,0.,0.]])
    result=reference_mask(coords,atm,bas,env,2,np.log(1e-12),CONSTANTS)
    assert np.array_equal(result,np.array([[1],[0]],dtype=np.int8))
