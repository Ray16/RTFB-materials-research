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

def analytical_hess_solver(pcmobj, dm, verbose=None):
    if (not pcmobj._intermediates or
        not any(isinstance(x, cupy.ndarray) for x in pcmobj._intermediates.values())):
        pcmobj.build()
    dm_cache = pcmobj._intermediates.get('dm', None)
    if dm_cache is not None and cupy.linalg.norm(dm_cache - dm) < 1e-10:
        pass
    else:
        pcmobj._get_vind(dm)
    mol = pcmobj.mol
    log = logger.new_logger(mol, verbose)
    t1 = log.init_timer()

    natom = mol.natm
    atmlst = range(natom) # Attention: This cannot be split

    gridslice    = pcmobj.surface['gslice_by_atom']
    v_grids      = pcmobj._intermediates['v_grids']
    q            = pcmobj._intermediates['q']
    f_epsilon    = pcmobj._intermediates['f_epsilon']
    if not pcmobj.if_method_in_CPCM_category:
        A = pcmobj._intermediates['A']
        D = pcmobj._intermediates['D']
        S = pcmobj._intermediates['S']

    ngrids = q.shape[0]

    vK_1 = pcmobj.left_solve_K(v_grids, K_transpose = True)

    if pcmobj.method.upper() in ['C-PCM', 'CPCM', 'COSMO']:
        dF, _ = get_dF_dA(pcmobj.surface, with_dA = False, surface_discretization_method = pcmobj.surface_discretization_method)
        dSii = get_dSii(pcmobj.surface, dF)

        available_gpu_memory = get_avail_mem()
        available_gpu_memory = int(available_gpu_memory * 0.8) # Don't use too much gpu memory
        dS_required_memory = 3 * ngrids * ngrids * cupy.dtype(cupy.float64).itemsize

        # dR = 0, dK = dS
        # d(S-1 R) = - S-1 dS S-1 R
        # d2(S-1 R) = (S-1 dS S-1 dS S-1 R) + (S-1 dS S-1 dS S-1 R) - (S-1 d2S S-1 R)
        if dS_required_memory <= available_gpu_memory:
            _, dS = get_dD_dS(pcmobj.surface, with_D=False, with_S=True)
            dSdx_dot_q = get_dS_dot_q(dS, dSii, q, atmlst, gridslice)
            S_1_dSdx_dot_q = einsum_ij_Adj_Adi_inverseK(pcmobj, dSdx_dot_q)
            del dSdx_dot_q
            VS_1_dot_dSdx = get_dST_dot_q(dS, dSii, vK_1, atmlst, gridslice)
            del dS
        else:
            dSdx_dot_q = get_dS_dot_q_lowmem(dSii, pcmobj.surface, q, atmlst, gridslice)
            S_1_dSdx_dot_q = einsum_ij_Adj_Adi_inverseK(pcmobj, dSdx_dot_q)
            del dSdx_dot_q
            VS_1_dot_dSdx = get_dST_dot_q_lowmem(dSii, pcmobj.surface, vK_1, atmlst, gridslice)

        del dSii
        d2e_from_d2KR = contract('Adi,BDi->ABdD', VS_1_dot_dSdx, S_1_dSdx_dot_q) * 2

        available_gpu_memory = get_avail_mem()
        available_gpu_memory = int(available_gpu_memory * 0.8) # Don't use too much gpu memory
        d2S_required_memory = (9 * ngrids * ngrids + 9 * ngrids * natom**2) * cupy.dtype(cupy.float64).itemsize

        if d2S_required_memory <= available_gpu_memory:
            _, d2S = get_d2D_d2S(pcmobj.surface, with_D=False, with_S=True)
            d2F, _ = get_d2F_d2A(pcmobj.surface, pcmobj.surface_discretization_method)
            d2Sii = get_d2Sii(pcmobj.surface, dF, d2F)
            del d2F
            v_d2S_q = get_v_dot_d2S_dot_q(d2S, d2Sii, vK_1, q, natom, gridslice)
            del d2S, d2Sii
        else:
            v_d2S_q = contract_d2S_offdiagonal(pcmobj.surface, vK_1, q, gridslice)
            v_d2S_q += contract_d2S_diagonal(pcmobj.surface, dF, vK_1 * q, pcmobj.surface_discretization_method)
        del dF
        d2e_from_d2KR -= v_d2S_q

        dK_1Rv = -S_1_dSdx_dot_q
        R = -f_epsilon
        dvK_1R = -einsum_Adi_ij_Adj_inverseK(VS_1_dot_dSdx, pcmobj) * R

    elif pcmobj.method.upper() in ['IEF-PCM', 'IEFPCM', 'SMD']:
        dD, dS = get_dD_dS(pcmobj.surface, with_D=True, with_S=True)
        dF, dA = get_dF_dA(pcmobj.surface, surface_discretization_method = pcmobj.surface_discretization_method)
        dSii = get_dSii(pcmobj.surface, dF)

        # dR = f_eps/(2*pi) * (dD*A + D*dA)
        # dK = dS - f_eps/(2*pi) * (dD*A*S + D*dA*S + D*A*dS)

        # d2R = f_eps/(2*pi) * (d2D*A + dD*dA + dD*dA + D*d2A)
        # d2K = d2S - f_eps/(2*pi) * (d2D*A*S + D*d2A*S + D*A*d2S + dD*dA*S + dD*dA*S + dD*A*dS + dD*A*dS + D*dA*dS + D*dA*dS)
        # The terms showing up twice on equation above (dD*dA + dD*dA for example) refer to dD/dx * dA/dy + dD/dy * dA/dx,
        # since D is not symmetric, they are not the same.

        # d(K-1 R) = - K-1 dK K-1 R + K-1 dR
        # d2(K-1 R) = (K-1 dK K-1 dK K-1 R) + (K-1 dK K-1 dK K-1 R) - (K-1 d2K K-1 R) - (K-1 dK K-1 dR)
        #             - (K-1 dK K-1 dR) + (K-1 d2R)
        f_eps_over_2pi = f_epsilon/(2.0*PI)

        dSdx_dot_q = get_dS_dot_q(dS, dSii, q, atmlst, gridslice)
        DA = D*A
        dKdx_dot_q = dSdx_dot_q - f_eps_over_2pi * contract('ij,Adj->Adi', DA, dSdx_dot_q)
        dAdx_dot_Sq = get_dA_dot_q(dA, S @ q, atmlst)
        dKdx_dot_q -= f_eps_over_2pi * contract('ij,Adj->Adi', D, dAdx_dot_Sq)
        AS = (A * S.T).T # It's just diag(A) @ S
        ASq = AS @ q
        dDdx_dot_ASq = get_dD_dot_q(dD, ASq, atmlst, gridslice, ngrids)
        dKdx_dot_q -= f_eps_over_2pi * dDdx_dot_ASq
        dDdx_dot_ASq = None

        K_1_dot_dKdx_dot_q = einsum_ij_Adj_Adi_inverseK(pcmobj, dKdx_dot_q)
        dKdx_dot_q = None

        vK_1_dot_dSdx = get_dST_dot_q(dS, dSii, vK_1, atmlst, gridslice)
        vK_1_dot_dKdx = vK_1_dot_dSdx
        vK_1_dot_dSdx = None
        vK_1_dot_dDdx = get_dDT_dot_q(dD, vK_1, atmlst, gridslice, ngrids)
        vK_1_dot_dKdx -= f_eps_over_2pi * contract('ij,Adj->Adi', AS.T, vK_1_dot_dDdx)
        AS = None
        vK_1D = D.T @ vK_1
        vK_1D_dot_dAdx = get_dA_dot_q(dA, vK_1D, atmlst)
        vK_1_dot_dKdx -= f_eps_over_2pi * contract('ij,Adj->Adi', S.T, vK_1D_dot_dAdx)
        vK_1DA = DA.T @ vK_1
        DA = None
        vK_1DA_dot_dSdx = get_dST_dot_q(dS, dSii, vK_1DA, atmlst, gridslice)
        dS = None
        dSii = None
        vK_1_dot_dKdx -= f_eps_over_2pi * vK_1DA_dot_dSdx
        vK_1DA_dot_dSdx = None

        d2e_from_d2KR  = contract('Adi,BDi->ABdD', vK_1_dot_dKdx, K_1_dot_dKdx_dot_q)
        d2e_from_d2KR += contract('Adi,BDi->BADd', vK_1_dot_dKdx, K_1_dot_dKdx_dot_q)

        d2F, d2A = get_d2F_d2A(pcmobj.surface, pcmobj.surface_discretization_method)
        vK_1_d2K_q  = get_v_dot_d2A_dot_q(d2A, vK_1D, S @ q)
        vK_1_d2R_V  = get_v_dot_d2A_dot_q(d2A, vK_1D, v_grids)
        d2A = None
        d2Sii = get_d2Sii(pcmobj.surface, dF, d2F)
        dF = None
        d2F = None
        d2D, d2S = get_d2D_d2S(pcmobj.surface, with_D=True, with_S=True)
        vK_1_d2K_q += get_v_dot_d2D_dot_q(d2D, vK_1, ASq, natom, gridslice)
        vK_1_d2R_V += get_v_dot_d2D_dot_q(d2D, vK_1, A * v_grids, natom, gridslice)
        d2D = None
        vK_1_d2K_q += get_v_dot_d2S_dot_q(d2S, d2Sii, vK_1DA, q, natom, gridslice)
        vK_1_d2K_q += contract('Adi,BDi->ABdD', vK_1_dot_dDdx, dAdx_dot_Sq)
        vK_1_d2K_q += contract('Adi,BDi->ABdD', vK_1_dot_dDdx * A, dSdx_dot_q)
        vK_1_d2K_q += contract('Adi,BDi->ABdD', vK_1D_dot_dAdx, dSdx_dot_q)
        vK_1_d2K_q += contract('Adi,BDi->BADd', vK_1_dot_dDdx, dAdx_dot_Sq)
        vK_1_d2K_q += contract('Adi,BDi->BADd', vK_1_dot_dDdx * A, dSdx_dot_q)
        vK_1_d2K_q += contract('Adi,BDi->BADd', vK_1D_dot_dAdx, dSdx_dot_q)
        vK_1_d2K_q *= -f_eps_over_2pi
        vK_1_d2K_q += get_v_dot_d2S_dot_q(d2S, d2Sii, vK_1, q, natom, gridslice)
        d2S = None
        d2Sii = None

        d2e_from_d2KR -= vK_1_d2K_q

        dAdx_dot_V = get_dA_dot_q(dA, v_grids, atmlst)
        dDdx_dot_AV = get_dD_dot_q(dD, A * v_grids, atmlst, gridslice, ngrids)
        dRdx_dot_V = f_eps_over_2pi * (dDdx_dot_AV + contract('ij,Adj->Adi', D, dAdx_dot_V))
        dDdx_dot_AV = None

        K_1_dot_dRdx_dot_V = einsum_ij_Adj_Adi_inverseK(pcmobj, dRdx_dot_V)
        dRdx_dot_V = None

        d2e_from_d2KR -= contract('Adi,BDi->ABdD', vK_1_dot_dKdx, K_1_dot_dRdx_dot_V)
        d2e_from_d2KR -= contract('Adi,BDi->BADd', vK_1_dot_dKdx, K_1_dot_dRdx_dot_V)

        vK_1_d2R_V += contract('Adi,BDi->ABdD', vK_1_dot_dDdx, dAdx_dot_V)
        vK_1_d2R_V += contract('Adi,BDi->BADd', vK_1_dot_dDdx, dAdx_dot_V)
        vK_1_d2R_V *= f_eps_over_2pi

        d2e_from_d2KR += vK_1_d2R_V

        dK_1Rv = -K_1_dot_dKdx_dot_q + K_1_dot_dRdx_dot_V

        VK_1D_dot_dAdx = get_dA_dot_q(dA, (D.T @ vK_1).T, atmlst)
        VK_1_dot_dDdx = get_dDT_dot_q(dD, vK_1, atmlst, gridslice, ngrids)
        VK_1_dot_dRdx = f_eps_over_2pi * (VK_1D_dot_dAdx + VK_1_dot_dDdx * A)

        dvK_1R = -_apply_R_batch(pcmobj,
            einsum_Adi_ij_Adj_inverseK(vK_1_dot_dKdx, pcmobj), True) + VK_1_dot_dRdx

    elif pcmobj.method.upper() in ['SS(V)PE']:
        dD, dS = get_dD_dS(pcmobj.surface, with_D=True, with_S=True)
        dF, dA = get_dF_dA(pcmobj.surface, surface_discretization_method = pcmobj.surface_discretization_method)
        dSii = get_dSii(pcmobj.surface, dF)

        # dR = f_eps/(2*pi) * (dD*A + D*dA)
        # dK = dS - f_eps/(4*pi) * (dD*A*S + D*dA*S + D*A*dS + dST*AT*DT + ST*dAT*DT + ST*AT*dDT)

        # d2R = f_eps/(2*pi) * (d2D*A + dD*dA + dD*dA + D*d2A)
        # d2K = d2S - f_eps/(4*pi) * (d2D*A*S + D*d2A*S + D*A*d2S + dD*dA*S + dD*dA*S + dD*A*dS + dD*A*dS + D*dA*dS + D*dA*dS
        #                           + d2ST*AT*DT + ST*d2AT*DT + ST*AT*d2DT + dST*dAT*DT + dST*dAT*DT + dST*AT*dDT + dST*AT*dDT + ST*dAT*dDT + ST*dAT*dDT)
        f_eps_over_2pi = f_epsilon/(2.0*PI)
        f_eps_over_4pi = f_epsilon/(4.0*PI)

        dSdx_dot_q = get_dS_dot_q(dS, dSii, q, atmlst, gridslice)
        DA = D*A
        dKdx_dot_q = dSdx_dot_q - f_eps_over_4pi * contract('ij,Adj->Adi', DA, dSdx_dot_q)
        dAdx_dot_Sq = get_dA_dot_q(dA, S @ q, atmlst)
        dKdx_dot_q -= f_eps_over_4pi * contract('ij,Adj->Adi', D, dAdx_dot_Sq)
        AS = (A * S.T).T # It's just diag(A) @ S
        ASq = AS @ q
        dDdx_dot_ASq = get_dD_dot_q(dD, ASq, atmlst, gridslice, ngrids)
        dKdx_dot_q -= f_eps_over_4pi * dDdx_dot_ASq
        dDdx_dot_ASq = None
        dDdxT_dot_q = get_dDT_dot_q(dD, q, atmlst, gridslice, ngrids)
        dKdx_dot_q -= f_eps_over_4pi * contract('ij,Adj->Adi', AS.T, dDdxT_dot_q)
        dAdxT_dot_DT_q = get_dA_dot_q(dA, D.T @ q, atmlst)
        dKdx_dot_q -= f_eps_over_4pi * contract('ij,Adj->Adi', S.T, dAdxT_dot_DT_q)
        AT_DT_q = DA.T @ q
        dSdxT_dot_AT_DT_q = get_dS_dot_q(dS, dSii, AT_DT_q, atmlst, gridslice)
        dKdx_dot_q -= f_eps_over_4pi * dSdxT_dot_AT_DT_q
        dSdxT_dot_AT_DT_q = None

        K_1_dot_dKdx_dot_q = einsum_ij_Adj_Adi_inverseK(pcmobj, dKdx_dot_q)
        dKdx_dot_q = None

        vK_1_dot_dSdx = get_dST_dot_q(dS, dSii, vK_1, atmlst, gridslice)
        vK_1_dot_dKdx = vK_1_dot_dSdx
        vK_1_dot_dSdx = None
        vK_1_dot_dDdx = get_dDT_dot_q(dD, vK_1, atmlst, gridslice, ngrids)
        vK_1_dot_dKdx -= f_eps_over_4pi * contract('ij,Adj->Adi', AS.T, vK_1_dot_dDdx)
        vK_1D_dot_dAdx = get_dA_dot_q(dA, D.T @ vK_1, atmlst)
        vK_1_dot_dKdx -= f_eps_over_4pi * contract('ij,Adj->Adi', S.T, vK_1D_dot_dAdx)
        vK_1DA = DA.T @ vK_1
        vK_1DA_dot_dSdx = get_dST_dot_q(dS, dSii, vK_1DA, atmlst, gridslice)
        vK_1_dot_dKdx -= f_eps_over_4pi * vK_1DA_dot_dSdx
        vK_1DA_dot_dSdx = None
        vK_1_dot_dSdxT = get_dS_dot_q(dS, dSii, vK_1, atmlst, gridslice)
        dS = None
        dSii = None
        vK_1_dot_dKdx -= f_eps_over_4pi * contract('ij,Adj->Adi', DA, vK_1_dot_dSdxT)
        DA = None
        vK_1_ST_dot_dAdxT = get_dA_dot_q(dA, (S @ vK_1).T, atmlst)
        vK_1_dot_dKdx -= f_eps_over_4pi * contract('ij,Adj->Adi', D, vK_1_ST_dot_dAdxT)
        vK_1_ST_AT = AS @ vK_1
        AS = None
        vK_1_ST_AT_dot_dDdxT = get_dD_dot_q(dD, vK_1_ST_AT, atmlst, gridslice, ngrids)
        vK_1_dot_dKdx -= f_eps_over_4pi * vK_1_ST_AT_dot_dDdxT
        vK_1_ST_AT_dot_dDdxT = None

        d2e_from_d2KR  = contract('Adi,BDi->ABdD', vK_1_dot_dKdx, K_1_dot_dKdx_dot_q)
        d2e_from_d2KR += contract('Adi,BDi->BADd', vK_1_dot_dKdx, K_1_dot_dKdx_dot_q)

        d2F, d2A = get_d2F_d2A(pcmobj.surface, pcmobj.surface_discretization_method)
        vK_1_d2K_q  = get_v_dot_d2A_dot_q(d2A, (D.T @ vK_1).T, S @ q)
        vK_1_d2K_q += get_v_dot_d2A_dot_q(d2A, (S @ vK_1).T, D.T @ q)
        vK_1_d2R_V  = get_v_dot_d2A_dot_q(d2A, (D.T @ vK_1).T, v_grids)
        d2A = None
        d2Sii = get_d2Sii(pcmobj.surface, dF, d2F)
        dF = None
        d2F = None
        d2D, d2S = get_d2D_d2S(pcmobj.surface, with_D=True, with_S=True)
        vK_1_d2K_q += get_v_dot_d2D_dot_q(d2D, vK_1, ASq, natom, gridslice)
        vK_1_d2K_q += get_v_dot_d2DT_dot_q(d2D, vK_1_ST_AT, q, natom, gridslice)
        vK_1_d2R_V += get_v_dot_d2D_dot_q(d2D, vK_1, A * v_grids, natom, gridslice)
        d2D = None
        vK_1_d2K_q += get_v_dot_d2S_dot_q(d2S, d2Sii, vK_1DA, q, natom, gridslice)
        vK_1_d2K_q += get_v_dot_d2ST_dot_q(d2S, d2Sii, vK_1, AT_DT_q, natom, gridslice)
        vK_1_d2K_q += contract('Adi,BDi->ABdD', vK_1_dot_dDdx, dAdx_dot_Sq)
        vK_1_d2K_q += contract('Adi,BDi->ABdD', vK_1_dot_dDdx * A, dSdx_dot_q)
        vK_1_d2K_q += contract('Adi,BDi->ABdD', vK_1D_dot_dAdx, dSdx_dot_q)
        vK_1_d2K_q += contract('Adi,BDi->ABdD', vK_1_dot_dSdxT, dAdxT_dot_DT_q)
        vK_1_d2K_q += contract('Adi,BDi->ABdD', vK_1_dot_dSdxT * A, dDdxT_dot_q)
        vK_1_d2K_q += contract('Adi,BDi->ABdD', vK_1_ST_dot_dAdxT, dDdxT_dot_q)
        vK_1_d2K_q += contract('Adi,BDi->BADd', vK_1_dot_dDdx, dAdx_dot_Sq)
        vK_1_d2K_q += contract('Adi,BDi->BADd', vK_1_dot_dDdx * A, dSdx_dot_q)
        vK_1_d2K_q += contract('Adi,BDi->BADd', vK_1D_dot_dAdx, dSdx_dot_q)
        vK_1_d2K_q += contract('Adi,BDi->BADd', vK_1_dot_dSdxT, dAdxT_dot_DT_q)
        vK_1_d2K_q += contract('Adi,BDi->BADd', vK_1_dot_dSdxT * A, dDdxT_dot_q)
        vK_1_d2K_q += contract('Adi,BDi->BADd', vK_1_ST_dot_dAdxT, dDdxT_dot_q)
        vK_1_d2K_q *= -f_eps_over_4pi
        vK_1_d2K_q += get_v_dot_d2S_dot_q(d2S, d2Sii, vK_1, q, natom, gridslice)
        d2S = None
        d2Sii = None

        d2e_from_d2KR -= vK_1_d2K_q

        dAdx_dot_V = get_dA_dot_q(dA, v_grids, atmlst)
        dDdx_dot_AV = get_dD_dot_q(dD, A * v_grids, atmlst, gridslice, ngrids)
        dRdx_dot_V = f_eps_over_2pi * (dDdx_dot_AV + contract('ij,Adj->Adi', D, dAdx_dot_V))
        dDdx_dot_AV = None

        K_1_dot_dRdx_dot_V = einsum_ij_Adj_Adi_inverseK(pcmobj, dRdx_dot_V)

        d2e_from_d2KR -= contract('Adi,BDi->ABdD', vK_1_dot_dKdx, K_1_dot_dRdx_dot_V)
        d2e_from_d2KR -= contract('Adi,BDi->BADd', vK_1_dot_dKdx, K_1_dot_dRdx_dot_V)

        vK_1_d2R_V += contract('Adi,BDi->ABdD', vK_1_dot_dDdx, dAdx_dot_V)
        vK_1_d2R_V += contract('Adi,BDi->BADd', vK_1_dot_dDdx, dAdx_dot_V)
        vK_1_d2R_V *= f_eps_over_2pi

        d2e_from_d2KR += vK_1_d2R_V

        dK_1Rv = -K_1_dot_dKdx_dot_q + K_1_dot_dRdx_dot_V

        VK_1D_dot_dAdx = get_dA_dot_q(dA, (D.T @ vK_1).T, atmlst)
        VK_1_dot_dDdx = get_dDT_dot_q(dD, vK_1, atmlst, gridslice, ngrids)
        VK_1_dot_dRdx = f_eps_over_2pi * (VK_1D_dot_dAdx + VK_1_dot_dDdx * A)

        DA = D*A
        R = -f_epsilon * (cupy.eye(DA.shape[0]) - 1.0/(2.0*PI)*DA)
        dvK_1R = -einsum_Adi_ij_Adj_inverseK(vK_1_dot_dKdx, pcmobj) @ R + VK_1_dot_dRdx

    else:
        raise RuntimeError(f"Unknown implicit solvent model: {pcmobj.method}")

    d2e = d2e_from_d2KR

    intopt_derivative = int3c1e.VHFOpt(mol)
    intopt_derivative.build(cutoff = 1e-14, aosym = False)

    dVdx = get_dvgrids(pcmobj, dm, range(mol.natm), intopt_derivative)
    d2e -= contract('Adi,BDi->BADd', dvK_1R, dVdx)
    d2e -= contract('Adi,BDi->ABdD', dVdx, dK_1Rv)

    d2e *= 0.5
    d2e = d2e.get()
    t1 = log.timer_debug1('solvent hessian d(V * dK-1R/dx * V)/dx contribution', *t1)
    return d2e
