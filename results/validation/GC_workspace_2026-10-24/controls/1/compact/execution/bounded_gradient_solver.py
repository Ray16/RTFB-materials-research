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

def grad_solver(pcmobj, dm, v_grids = None, v_grids_l = None, q = None):
    '''
    dE = 0.5*v* d(K^-1 R) *v + q*dv
    v^T* d(K^-1 R)v = v^T*K^-1(dR - dK K^-1R)v = v^T K^-1(dR v - dK q)
    '''
    mol = pcmobj.mol
    log = logger.new_logger(mol, mol.verbose)
    t1 = log.init_timer()
    if (not pcmobj._intermediates or
        not any(isinstance(x, cupy.ndarray) for x in pcmobj._intermediates.values())):
        pcmobj.build()
    dm_cache = pcmobj._intermediates.get('dm', None)
    if dm_cache is not None and cupy.linalg.norm(dm_cache - dm) < 1e-10:
        pass
    else:
        pcmobj._get_vind(dm)

    gridslice    = pcmobj.surface['gslice_by_atom']
    if v_grids is None:
        v_grids = pcmobj._intermediates['v_grids']
    if v_grids_l is None:
        v_grids_l = pcmobj._intermediates['v_grids']
    if q is None:
        q = pcmobj._intermediates['q']
    f_epsilon = pcmobj._intermediates['f_epsilon']
    lowmem_mode = True
    if (not pcmobj.if_method_in_CPCM_category) and (not lowmem_mode):
        A = pcmobj._intermediates['A']
        D = pcmobj._intermediates['D']
        S = pcmobj._intermediates['S']

    if pcmobj.frozen_dm0_for_finite_difference_without_response is not None:
        # Note: The v_grids computed above actually use frozen_dm0 as input, so it's actually v_grids_right
        v_grids_l = pcmobj._get_vgrids(dm, with_nuc = True)[0]

        # TODO: In the case where left and right dm are not the same,
        #       we indeed need to compute the derivative of 0.5 * (K^-1 R + R^T (K^-1)^T).
        #       This is not the same as the derivative of K^-1 R, if IEFPCM or SSVPE is used.
        #       However the difference is too small, and in the use case of computing polarizability derivative,
        #       we are not able to observe the difference.
        #       If there are other use cases where the error is more significant,
        #       we probably need to fix this.
        #       This is the only term affected by this problem in energy and gradient calculation.
        #       In hessian calculation, similar problems occur in energy 2nd derivative and Fock derivative terms.

    vK_1 = pcmobj.left_solve_K(v_grids_l, K_transpose = True)

    def contract_bra(a, B, c):
        ''' i,xij,j->jx '''
        tmp = a.dot(B)
        return (tmp*c).T

    def contract_ket(a, B, c):
        ''' i,xij,j->ix '''
        tmp = B.dot(c)
        return (a*tmp).T

    natm = natm_without_ghost(pcmobj.mol)
    de = cupy.zeros([natm, 3])
    if pcmobj.method.upper() in ['C-PCM', 'CPCM', 'COSMO']:
        # dR = 0, dK = dS
        de_dS  = 0.5 * vK_1.reshape(-1, 1) * left_multiply_dS_offdiagonal(pcmobj.surface, q, stream = None)
        de_dS -= 0.5 * q.reshape(-1, 1) * left_multiply_dS_offdiagonal(pcmobj.surface, vK_1, transpose = True, stream = None)
        de -= cupy.asarray([cupy.sum(de_dS[p0:p1], axis=0) for p0,p1 in gridslice])

        de -= 0.5 * contract_dSii(pcmobj.surface, vK_1 * q, pcmobj.surface_discretization_method)

    elif pcmobj.method.upper() in ['IEF-PCM', 'IEFPCM', 'SMD']:
        # dR = f_eps/(2*pi) * (dD*A + D*dA),
        # dK = dS - f_eps/(2*pi) * (dD*A*S + D*dA*S + D*A*dS)
        fac = f_epsilon/(2.0*PI)

        if not lowmem_mode:
            dF, dA = get_dF_dA(pcmobj.surface, surface_discretization_method = pcmobj.surface_discretization_method)
            dSii = get_dSii(pcmobj.surface, dF)
            dF = None

            dD, dS = get_dD_dS(pcmobj.surface, with_D=True, with_S=True)

            Av = A*v_grids
            de_dR  = 0.5*fac * contract_ket(vK_1, dD, Av)
            de_dR -= 0.5*fac * contract_bra(vK_1, dD, Av)
            de_dR  = cupy.asarray([cupy.sum(de_dR[p0:p1], axis=0) for p0,p1 in gridslice])

            vK_1_D = vK_1.dot(D)
            vK_1_Dv = vK_1_D * v_grids
            de_dR += 0.5*fac * contract('j,xjn->nx', vK_1_Dv, dA)

            de_dS0  = 0.5*contract_ket(vK_1, dS, q)
            de_dS0 -= 0.5*contract_bra(vK_1, dS, q)
            de_dS0  = cupy.asarray([cupy.sum(de_dS0[p0:p1], axis=0) for p0,p1 in gridslice])

            vK_1_q = vK_1 * q
            de_dS0 += 0.5*contract('i,xin->nx', vK_1_q, dSii)

            vK_1_DA = vK_1_D*A
            de_dS1  = 0.5*contract_ket(vK_1_DA, dS, q)
            de_dS1 -= 0.5*contract_bra(vK_1_DA, dS, q)
            de_dS1  = cupy.asarray([cupy.sum(de_dS1[p0:p1], axis=0) for p0,p1 in gridslice])

            vK_1_DAq = vK_1_DA*q
            de_dS1 += 0.5*contract('j,xjn->nx', vK_1_DAq, dSii)

            Sq = cupy.dot(S,q)
            ASq = A*Sq
            de_dD  = 0.5*contract_ket(vK_1, dD, ASq)
            de_dD -= 0.5*contract_bra(vK_1, dD, ASq)
            de_dD  = cupy.asarray([cupy.sum(de_dD[p0:p1], axis=0) for p0,p1 in gridslice])

            de_dA = 0.5*contract('j,xjn->nx', vK_1_D*Sq, dA)   # 0.5*cupy.einsum('j,xjn,j->nx', vK_1_D, dA, Sq)

            de_dK = de_dS0 - fac * (de_dD + de_dA + de_dS1)
            de += de_dR - de_dK
        else:
            A = pcmobj._intermediates['A']
            Av = A*v_grids
            de_dR  = 0.5*fac * vK_1.reshape(-1, 1) * left_multiply_dD(pcmobj.surface, Av, stream = None)
            de_dR -= 0.5*fac * Av.reshape(-1, 1) * left_multiply_dD(pcmobj.surface, vK_1, transpose = True, stream = None)
            de_dR  = cupy.asarray([cupy.sum(de_dR[p0:p1], axis=0) for p0,p1 in gridslice])

            de_dS0  = 0.5 * vK_1.reshape(-1, 1) * left_multiply_dS_offdiagonal(pcmobj.surface, q, stream = None)
            de_dS0 -= 0.5 * q.reshape(-1, 1) * left_multiply_dS_offdiagonal(pcmobj.surface, vK_1, transpose = True, stream = None)
            de_dS0  = cupy.asarray([cupy.sum(de_dS0[p0:p1], axis=0) for p0,p1 in gridslice])

            vK_1_D = left_multiply_D(pcmobj.surface, vK_1, transpose = True)
            vK_1_DA = vK_1_D * A
            de_dS1  = 0.5 * vK_1_DA.reshape(-1, 1) * left_multiply_dS_offdiagonal(pcmobj.surface, q, stream = None)
            de_dS1 -= 0.5 * q.reshape(-1, 1) * left_multiply_dS_offdiagonal(pcmobj.surface, vK_1_DA, transpose = True, stream = None)
            de_dS1  = cupy.asarray([cupy.sum(de_dS1[p0:p1], axis=0) for p0,p1 in gridslice])

            Sq = left_multiply_S(pcmobj.surface, q)
            ASq = A * Sq
            de_dD  = 0.5 * vK_1.reshape(-1, 1) * left_multiply_dD(pcmobj.surface, ASq, stream = None)
            de_dD -= 0.5 * ASq.reshape(-1, 1) * left_multiply_dD(pcmobj.surface, vK_1, transpose = True, stream = None)
            de_dD  = cupy.asarray([cupy.sum(de_dD[p0:p1], axis=0) for p0,p1 in gridslice])

            de_dK = de_dS0 - fac * (de_dD + de_dS1)
            de_dK += 0.5 * contract_dSii(pcmobj.surface, vK_1 * q - fac * vK_1_DA * q, pcmobj.surface_discretization_method)

            de += de_dR - de_dK
            de += 0.5*fac * contract_dA(pcmobj.surface, vK_1_D * v_grids + vK_1_D * Sq,
                                        pcmobj.surface_discretization_method) # First term from de_dR, second from de_dK

    elif pcmobj.method.upper() in [ 'SS(V)PE' ]:
        # dR = f_eps/(2*pi) * (dD*A + D*dA),
        # dK = dS - f_eps/(2*pi) * (dD*A*S + D*dA*S + D*A*dS)
        fac = f_epsilon/(2.0*PI)

        if not lowmem_mode:
            dF, dA = get_dF_dA(pcmobj.surface, surface_discretization_method = pcmobj.surface_discretization_method)
            dSii = get_dSii(pcmobj.surface, dF)
            dF = None

            dD, dS = get_dD_dS(pcmobj.surface, with_D=True, with_S=True)

            Av = A*v_grids
            de_dR  = 0.5*fac * contract_ket(vK_1, dD, Av)
            de_dR -= 0.5*fac * contract_bra(vK_1, dD, Av)
            de_dR  = cupy.asarray([cupy.sum(de_dR[p0:p1], axis=0) for p0,p1 in gridslice])

            vK_1_D = vK_1.dot(D)
            vK_1_Dv = vK_1_D * v_grids
            de_dR += 0.5*fac * contract('j,xjn->nx', vK_1_Dv, dA)

            de_dS0  = 0.5*contract_ket(vK_1, dS, q)
            de_dS0 -= 0.5*contract_bra(vK_1, dS, q)
            de_dS0  = cupy.asarray([cupy.sum(de_dS0[p0:p1], axis=0) for p0,p1 in gridslice])

            vK_1_q = vK_1 * q
            de_dS0 += 0.5*contract('i,xin->nx', vK_1_q, dSii)

            vK_1_DA = vK_1_D*A
            de_dS1  = 0.5*contract_ket(vK_1_DA, dS, q)
            de_dS1 -= 0.5*contract_bra(vK_1_DA, dS, q)
            de_dS1  = cupy.asarray([cupy.sum(de_dS1[p0:p1], axis=0) for p0,p1 in gridslice])
            vK_1_DAq = vK_1_DA*q
            de_dS1 += 0.5*contract('j,xjn->nx', vK_1_DAq, dSii)

            DT_q = cupy.dot(D.T, q)
            ADT_q = A * DT_q
            de_dS1_T  = 0.5*contract_ket(vK_1, dS, ADT_q)
            de_dS1_T -= 0.5*contract_bra(vK_1, dS, ADT_q)
            de_dS1_T  = cupy.asarray([cupy.sum(de_dS1_T[p0:p1], axis=0) for p0,p1 in gridslice])
            vK_1_ADT_q = vK_1 * ADT_q
            de_dS1_T += 0.5*contract('j,xjn->nx', vK_1_ADT_q, dSii)

            Sq = cupy.dot(S,q)
            ASq = A*Sq
            de_dD  = 0.5*contract_ket(vK_1, dD, ASq)
            de_dD -= 0.5*contract_bra(vK_1, dD, ASq)
            de_dD  = cupy.asarray([cupy.sum(de_dD[p0:p1], axis=0) for p0,p1 in gridslice])

            vK_1_S = cupy.dot(vK_1, S)
            vK_1_SA = vK_1_S * A
            de_dD_T  = 0.5*contract_ket(vK_1_SA, -dD.transpose(0,2,1), q)
            de_dD_T -= 0.5*contract_bra(vK_1_SA, -dD.transpose(0,2,1), q)
            de_dD_T  = cupy.asarray([cupy.sum(de_dD_T[p0:p1], axis=0) for p0,p1 in gridslice])

            de_dA = 0.5*contract('j,xjn->nx', vK_1_D*Sq, dA)   # 0.5*cupy.einsum('j,xjn,j->nx', vK_1_D, dA, Sq)

            de_dA_T = 0.5*contract('j,xjn->nx', vK_1_S*DT_q, dA)

            de_dK = de_dS0 - 0.5 * fac * (de_dD + de_dA + de_dS1 + de_dD_T + de_dA_T + de_dS1_T)
            de += de_dR - de_dK
        else:
            A = pcmobj._intermediates['A']
            Av = A*v_grids
            de_dR  = 0.5*fac * vK_1.reshape(-1, 1) * left_multiply_dD(pcmobj.surface, Av, stream = None)
            de_dR -= 0.5*fac * Av.reshape(-1, 1) * left_multiply_dD(pcmobj.surface, vK_1, transpose = True, stream = None)
            de_dR  = cupy.asarray([cupy.sum(de_dR[p0:p1], axis=0) for p0,p1 in gridslice])

            de_dS0  = 0.5 * vK_1.reshape(-1, 1) * left_multiply_dS_offdiagonal(pcmobj.surface, q, stream = None)
            de_dS0 -= 0.5 * q.reshape(-1, 1) * left_multiply_dS_offdiagonal(pcmobj.surface, vK_1, transpose = True, stream = None)
            de_dS0  = cupy.asarray([cupy.sum(de_dS0[p0:p1], axis=0) for p0,p1 in gridslice])

            vK_1_D = left_multiply_D(pcmobj.surface, vK_1, transpose = True)
            vK_1_DA = vK_1_D * A
            de_dS1  = 0.5 * vK_1_DA.reshape(-1, 1) * left_multiply_dS_offdiagonal(pcmobj.surface, q, stream = None)
            de_dS1 -= 0.5 * q.reshape(-1, 1) * left_multiply_dS_offdiagonal(pcmobj.surface, vK_1_DA, transpose = True, stream = None)
            de_dS1  = cupy.asarray([cupy.sum(de_dS1[p0:p1], axis=0) for p0,p1 in gridslice])

            DT_q = left_multiply_D(pcmobj.surface, q, transpose = True)
            ADT_q = A * DT_q
            de_dS1_T  = 0.5 * vK_1.reshape(-1, 1) * left_multiply_dS_offdiagonal(pcmobj.surface, ADT_q, stream = None)
            de_dS1_T -= 0.5 * ADT_q.reshape(-1, 1) * left_multiply_dS_offdiagonal(pcmobj.surface, vK_1, transpose = True, stream = None)
            de_dS1_T  = cupy.asarray([cupy.sum(de_dS1_T[p0:p1], axis=0) for p0,p1 in gridslice])

            Sq = left_multiply_S(pcmobj.surface, q)
            ASq = A * Sq
            de_dD  = 0.5 * vK_1.reshape(-1, 1) * left_multiply_dD(pcmobj.surface, ASq, stream = None)
            de_dD -= 0.5 * ASq.reshape(-1, 1) * left_multiply_dD(pcmobj.surface, vK_1, transpose = True, stream = None)
            de_dD  = cupy.asarray([cupy.sum(de_dD[p0:p1], axis=0) for p0,p1 in gridslice])

            vK_1_S = left_multiply_S(pcmobj.surface, vK_1, transpose = True)
            vK_1_SA = vK_1_S * A
            # Attention on the transposed order here
            de_dD_T  = 0.5 * q.reshape(-1, 1) * left_multiply_dD(pcmobj.surface, vK_1_SA, stream = None)
            de_dD_T -= 0.5 * vK_1_SA.reshape(-1, 1) * left_multiply_dD(pcmobj.surface, q, transpose = True, stream = None)
            de_dD_T  = cupy.asarray([cupy.sum(de_dD_T[p0:p1], axis=0) for p0,p1 in gridslice])

            de_dK = de_dS0 - 0.5*fac * (de_dD + de_dS1 + de_dD_T + de_dS1_T)
            de_dK += 0.5 * contract_dSii(pcmobj.surface, vK_1 * q - 0.5*fac * (vK_1_DA * q + vK_1 * ADT_q), pcmobj.surface_discretization_method)

            de += de_dR - de_dK
            de += 0.5*fac * contract_dA(pcmobj.surface, vK_1_D * v_grids + 0.5 * (vK_1_D * Sq + vK_1_S * DT_q),
                                        pcmobj.surface_discretization_method) # First term from de_dR, second from de_dK

    else:
        raise RuntimeError(f"Unknown implicit solvent model: {pcmobj.method}")

    if pcmobj.frozen_dm0_for_finite_difference_without_response is not None:
        # Refer to the comments in gpu4pyscf/solvent/pcm.py::_get_vind()
        de *= 2

    de = gradient_nonghost_to_all(mol, de)

    t1 = log.timer_debug1('grad solver', *t1)
    return de.get()
