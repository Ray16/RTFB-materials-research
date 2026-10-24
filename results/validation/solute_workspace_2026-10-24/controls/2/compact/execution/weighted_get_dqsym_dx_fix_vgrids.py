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

def get_dqsym_dx_fix_vgrids(pcmobj, atmlst):
    assert pcmobj._intermediates is not None

    gridslice    = pcmobj.surface['gslice_by_atom']
    v_grids      = pcmobj._intermediates['v_grids']
    q            = pcmobj._intermediates['q']
    q_sym        = pcmobj._intermediates['q_sym']
    f_epsilon    = pcmobj._intermediates['f_epsilon']
    if not pcmobj.if_method_in_CPCM_category:
        A = pcmobj._intermediates['A']
        D = pcmobj._intermediates['D']
        S = pcmobj._intermediates['S']

    ngrids = q_sym.shape[0]

    if pcmobj.method.upper() in ['C-PCM', 'CPCM', 'COSMO']:
        dF, _ = get_dF_dA(pcmobj.surface, with_dA = False, surface_discretization_method = pcmobj.surface_discretization_method)
        dSii = get_dSii(pcmobj.surface, dF)
        del dF

        available_gpu_memory = get_avail_mem()
        available_gpu_memory = int(available_gpu_memory * 0.9) # Don't use too much gpu memory
        dS_required_memory = 3 * ngrids * ngrids * cupy.dtype(cupy.float64).itemsize

        # dR = 0, dK = dS
        if dS_required_memory <= available_gpu_memory:
            _, dS = get_dD_dS(pcmobj.surface, with_D=False, with_S=True)
            dSdx_dot_q = get_dS_dot_q(dS, dSii, q_sym, atmlst, gridslice)
            del dS
        else:
            dSdx_dot_q = get_dS_dot_q_lowmem(dSii, pcmobj.surface, q_sym, atmlst, gridslice, stream = None)

        dqdx_fix_Vq = einsum_ij_Adj_Adi_inverseK(pcmobj, dSdx_dot_q)

    elif pcmobj.method.upper() in ['IEF-PCM', 'IEFPCM', 'SMD']:
        dF, dA = get_dF_dA(pcmobj.surface, surface_discretization_method = pcmobj.surface_discretization_method)
        dSii = get_dSii(pcmobj.surface, dF)
        dF = None

        dD, dS = get_dD_dS(pcmobj.surface, with_D=True, with_S=True)

        # dR = f_eps/(2*pi) * (dD*A + D*dA)
        # dK = dS - f_eps/(2*pi) * (dD*A*S + D*dA*S + D*A*dS)
        f_eps_over_2pi = f_epsilon/(2.0*PI)

        dSdx_dot_q = get_dS_dot_q(dS, dSii, q, atmlst, gridslice)

        DA = _scaled_matrix(D, right=A)
        dKdx_dot_q = dSdx_dot_q - f_eps_over_2pi * contract('ij,Adj->Adi', DA, dSdx_dot_q)

        dAdx_dot_Sq = get_dA_dot_q(dA, S @ q, atmlst)
        dKdx_dot_q -= f_eps_over_2pi * contract('ij,Adj->Adi', D, dAdx_dot_Sq)

        AS = _scaled_matrix(S, left=A)
        dDdx_dot_ASq = get_dD_dot_q(dD, AS @ q, atmlst, gridslice, ngrids)
        dKdx_dot_q -= f_eps_over_2pi * dDdx_dot_ASq

        dqdx_fix_Vq = -einsum_ij_Adj_Adi_inverseK(pcmobj, dKdx_dot_q)

        dAdx_dot_V = get_dA_dot_q(dA, v_grids, atmlst)

        dDdx_dot_AV = get_dD_dot_q(dD, A * v_grids, atmlst, gridslice, ngrids)

        dRdx_dot_V = f_eps_over_2pi * (dDdx_dot_AV + contract('ij,Adj->Adi', D, dAdx_dot_V))
        dqdx_fix_Vq += einsum_ij_Adj_Adi_inverseK(pcmobj, dRdx_dot_V)

        invKT_V = pcmobj.left_solve_K(v_grids, K_transpose = True)
        dDdxT_dot_invKT_V = get_dDT_dot_q(dD, invKT_V, atmlst, gridslice, ngrids)

        DT_invKT_V = D.T @ invKT_V
        dAdxT_dot_DT_invKT_V = get_dA_dot_q(dA, DT_invKT_V, atmlst)
        dqdx_fix_Vq += f_eps_over_2pi * (contract('i,Adi->Adi', A, dDdxT_dot_invKT_V) + dAdxT_dot_DT_invKT_V)

        dSdxT_dot_invKT_V = get_dST_dot_q(dS, dSii, invKT_V, atmlst, gridslice)
        dKdxT_dot_invKT_V = dSdxT_dot_invKT_V

        dKdxT_dot_invKT_V -= f_eps_over_2pi * contract('ij,Adj->Adi', AS.T, dDdxT_dot_invKT_V)
        dKdxT_dot_invKT_V -= f_eps_over_2pi * contract('ij,Adj->Adi', S.T, dAdxT_dot_DT_invKT_V)

        dSdxT_dot_AT_DT_invKT_V = get_dST_dot_q(dS, dSii, DA.T @ invKT_V, atmlst, gridslice)
        dKdxT_dot_invKT_V -= f_eps_over_2pi * dSdxT_dot_AT_DT_invKT_V
        invKT_dKdxT_dot_invKT_V = einsum_ij_Adj_Adi_inverseK(pcmobj, dKdxT_dot_invKT_V, K_transpose = True)

        del DA, AS
        dqdx_fix_Vq -= _apply_R_batch(pcmobj, invKT_dKdxT_dot_invKT_V, True)

        dqdx_fix_Vq *= -0.5

    elif pcmobj.method.upper() in ['SS(V)PE']:
        dF, dA = get_dF_dA(pcmobj.surface, surface_discretization_method = pcmobj.surface_discretization_method)
        dSii = get_dSii(pcmobj.surface, dF)
        dF = None

        dD, dS = get_dD_dS(pcmobj.surface, with_D=True, with_S=True)

        f_eps_over_4pi = f_epsilon/(4.0*PI)
        DA = D*A

        def dK_dot_q(q):
            dSdx_dot_q = get_dS_dot_q(dS, dSii, q, atmlst, gridslice)

            dKdx_dot_q = dSdx_dot_q - f_eps_over_4pi * contract('ij,Adj->Adi', DA, dSdx_dot_q)

            dAdx_dot_Sq = get_dA_dot_q(dA, S @ q, atmlst)
            dKdx_dot_q -= f_eps_over_4pi * contract('ij,Adj->Adi', D, dAdx_dot_Sq)

            AS = (A * S.T).T # It's just diag(A) @ S
            dDdx_dot_ASq = get_dD_dot_q(dD, AS @ q, atmlst, gridslice, ngrids)
            dKdx_dot_q -= f_eps_over_4pi * dDdx_dot_ASq

            dDdxT_dot_q = get_dDT_dot_q(dD, q, atmlst, gridslice, ngrids)
            dKdx_dot_q -= f_eps_over_4pi * contract('ij,Adj->Adi', AS.T, dDdxT_dot_q)

            dAdxT_dot_DT_q = get_dA_dot_q(dA, D.T @ q, atmlst)
            dKdx_dot_q -= f_eps_over_4pi * contract('ij,Adj->Adi', S.T, dAdxT_dot_DT_q)

            dSdxT_dot_AT_DT_q = get_dST_dot_q(dS, dSii, DA.T @ q, atmlst, gridslice)
            dKdx_dot_q -= f_eps_over_4pi * dSdxT_dot_AT_DT_q

            return dKdx_dot_q

        f_eps_over_2pi = f_epsilon/(2.0*PI)

        dKdx_dot_q = dK_dot_q(q)
        dqdx_fix_Vq = -einsum_ij_Adj_Adi_inverseK(pcmobj, dKdx_dot_q)

        dAdx_dot_V = get_dA_dot_q(dA, v_grids, atmlst)

        dDdx_dot_AV = get_dD_dot_q(dD, A * v_grids, atmlst, gridslice, ngrids)

        dRdx_dot_V = f_eps_over_2pi * (dDdx_dot_AV + contract('ij,Adj->Adi', D, dAdx_dot_V))
        dqdx_fix_Vq += einsum_ij_Adj_Adi_inverseK(pcmobj, dRdx_dot_V)

        invKT_V = pcmobj.left_solve_K(v_grids, K_transpose = True)
        dDdxT_dot_invKT_V = get_dDT_dot_q(dD, invKT_V, atmlst, gridslice, ngrids)

        DT_invKT_V = D.T @ invKT_V
        dAdxT_dot_DT_invKT_V = get_dA_dot_q(dA, DT_invKT_V, atmlst)
        dqdx_fix_Vq += f_eps_over_2pi * (contract('i,Adi->Adi', A, dDdxT_dot_invKT_V) + dAdxT_dot_DT_invKT_V)

        dKdx_dot_invKT_V = dK_dot_q(invKT_V)
        invKT_dKdx_dot_invKT_V = einsum_ij_Adj_Adi_inverseK(pcmobj, dKdx_dot_invKT_V, K_transpose = True)

        R = -f_epsilon * (cupy.eye(DA.shape[0]) - 1.0/(2.0*PI)*DA)
        dqdx_fix_Vq += -contract('ij,Adj->Adi', R.T, invKT_dKdx_dot_invKT_V)

        dqdx_fix_Vq *= -0.5

    else:
        raise RuntimeError(f"Unknown implicit solvent model: {pcmobj.method}")

    return dqdx_fix_Vq
