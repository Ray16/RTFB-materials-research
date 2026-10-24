# Copyright 2021-2024 The PySCF Developers. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

def get_exc_full_response(ni, mol, grids, xc_code, dms, relativity=0, hermi=1,
                          max_memory=2000, verbose=None):
    '''Full response including the response of the grids'''
    log = logger.new_logger(mol, verbose)
    t0 = log.init_timer()
    ni = numint.NumInt() # Don't mess up with the old numint object
    xctype = ni._xc_type(xc_code)

    grids = grids.copy()
    grids.build(sort_grids_of_each_atom = True)
    ngrids = grids.coords.shape[0]

    ni.gdftopt = None
    ni.build(mol, grids.coords)
    opt = ni.gdftopt

    natm = mol.natm
    mol = None
    _sorted_mol = opt._sorted_mol
    nao = _sorted_mol.nao
    dms = cupy.asarray(dms)
    assert dms.ndim == 2
    dms = opt.sort_orbitals(dms, axis=[0,1])

    de_grid_response_rho = cupy.zeros((natm, 3))
    dvmat_orbital_response = cupy.zeros((3, nao, nao))
    dm_mask_buf = cupy.empty(nao*nao)

    if xctype == 'LDA':
        ao_deriv = 0
        ncomp = 1
    elif xctype == 'GGA':
        ao_deriv = 1
        ncomp = 4
    elif xctype == 'MGGA':
        ao_deriv = 1
        ncomp = 5
    else:
        raise NotImplementedError(f"Unrecognized xctype = {xctype}")

    rho = cupy.empty([ncomp, ngrids])
    g1 = 0
    for ao, idx, weight, _ in ni.block_loop(_sorted_mol, grids, deriv = ao_deriv, strict_grid_order = True):
        g0, g1 = g1, g1 + weight.size
        dms_masked = take_last2d(dms, idx, out=dm_mask_buf)
        rho[:, g0:g1] = numint.eval_rho(_sorted_mol, ao, dms_masked, xctype = xctype, hermi = 1)
    assert g1 == ngrids
    _workspace_trace('density_loop_complete')
    cupy.cuda.get_current_stream().synchronize()
    del ao, idx, weight, dms_masked
    _workspace_trace('density_batch_released')

    exc, vxc = ni.eval_xc_eff(xc_code, rho, 1, xctype=xctype, spin=0)[:2]
    wv = grids.weights * vxc
    nonzero_weight_mask = cupy.abs(grids.weights) > 1e-14

    from gpu4pyscf.hessian.rks import get_dweight_dA

    de_grid_response_weight = cupy.zeros((natm, 3))
    dweightdA_right = rho[0] * exc
    del rho

    available_gpu_memory = get_avail_mem()
    available_gpu_memory = int(available_gpu_memory * 0.1) # Don't use too much gpu memory
    ao_nbytes_per_grid = ((2*3) * natm) * 8
    ngrids_per_batch = int(available_gpu_memory / ao_nbytes_per_grid)
    if ngrids_per_batch < 16:
        raise MemoryError(f"Out of GPU memory for XC energy first derivative, available gpu memory = {get_avail_mem()}"
                          f" bytes, nao = {nao}, natm = {natm}, ngrids (nonzero rho) = {ngrids}")
    ngrids_per_batch = (ngrids_per_batch + 16 - 1) // 16 * 16
    ### Don't split the batch too small for get_dweight_dA()
    # ngrids_per_batch = min(ngrids_per_batch, MIN_BLK_SIZE)

    for g0 in range(0, ngrids, ngrids_per_batch):
        g1 = min(g0 + ngrids_per_batch, ngrids)
        dweight_dA = get_dweight_dA(_sorted_mol, grids, (g0,g1))
        de_grid_response_weight += cupy.einsum("Adg->Ad", dweight_dA * dweightdA_right[g0:g1])
    del dweight_dA
    del dweightdA_right
    del exc

    g0 = 0
    for ao, idx, weight, _ in ni.block_loop(_sorted_mol, grids, nao, ao_deriv + 1, strict_grid_order = True):
        g1 = g0 + weight.shape[0]

        ao = ao[:, :, nonzero_weight_mask[g0:g1]]

        if ao.size == 0:
            g0 = g1
            continue

        dms_masked = take_last2d(dms, idx, out=dm_mask_buf)

        i_atom = int(grids.atm_idx[g0])
        assert cupy.max(cupy.abs(grids.atm_idx[g0:g1] - i_atom)) == 0 # Guaranteed by grids.build(sort_grids_of_each_atom = True)

        if xctype == 'LDA':
            split_wv = cupy.ascontiguousarray(wv[:, g0:g1][:, nonzero_weight_mask[g0:g1]])

            aow = numint._scale_ao(ao[0], split_wv[0])
            vtmp = _d1_dot_(ao[1:4], aow.T)

            dvmat_orbital_response[:, idx[:,None], idx] += vtmp
            de_grid_response_rho[i_atom] += cupy.einsum('xij,ji->x', vtmp, dms_masked) * 2

        elif xctype == 'GGA':
            split_wv = cupy.ascontiguousarray(wv[:, g0:g1][:, nonzero_weight_mask[g0:g1]])
            split_wv[0] *= .5

            vtmp = _gga_grad_sum_(ao, split_wv[:4])

            dvmat_orbital_response[:, idx[:,None], idx] += vtmp
            de_grid_response_rho[i_atom] += cupy.einsum('xij,ji->x', vtmp, dms_masked) * 2

        elif xctype == 'NLC':
            raise ValueError("You see a bug, please report to the developer team.")

        elif xctype == 'MGGA':
            split_wv = cupy.ascontiguousarray(wv[:, g0:g1][:, nonzero_weight_mask[g0:g1]])
            split_wv[0] *= .5
            split_wv[4] *= .5 # for the factor 1/2 in tau

            vtmp = _gga_grad_sum_(ao, split_wv[:4])
            _tau_grad_dot_(ao, split_wv[4], accumulate=True, out=vtmp)

            dvmat_orbital_response[:, idx[:,None], idx] += vtmp
            de_grid_response_rho[i_atom] += cupy.einsum('xij,ji->x', vtmp, dms_masked) * 2

        else:
            raise NotImplementedError(f"Unrecognized xctype = {xctype}")

        g0 = g1
    assert g1 == ngrids

    excsum = de_grid_response_weight + de_grid_response_rho
    excsum = excsum.get()
    # - sign because nabla_X = -nabla_x
    excsum -= rhf_grad.contract_h1e_dm(opt._sorted_mol, dvmat_orbital_response, dms, hermi=1)

    log.timer_debug1('rks grad vxc full response', *t0)
    return excsum, 0
