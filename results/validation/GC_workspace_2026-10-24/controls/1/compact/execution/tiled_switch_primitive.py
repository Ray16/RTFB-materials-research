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

def get_d2F_d2A(surface, surface_discretization_method = "SWIG"):
    '''
    Notations adopted from
    J. Chem. Phys. 133, 244111 (2010), Appendix C
    '''
    atom_coords = surface['atom_coords']
    grid_coords = surface['grid_coords']
    switch_fun  = surface['switch_fun']
    area        = surface['area']
    if surface_discretization_method.upper() == "SWIG":
        R_in_J = surface['R_in_J']
        R_sw_J = surface['R_sw_J']
    elif surface_discretization_method.upper() == "ISWIG":
        charge_exp = surface['charge_exp']
        R_J = surface['R_J']

    ngrids = grid_coords.shape[0]
    natom = atom_coords.shape[0]
    d2F = cupy.zeros([ngrids, natom, natom, 3, 3])
    d2A = cupy.zeros([ngrids, natom, natom, 3, 3])

    for i_grid_atom in range(natom):
        p0,p1 = surface['gslice_by_atom'][i_grid_atom]
        if p0 == p1:
            continue
        coords = grid_coords[p0:p1]
        si_rJ = cupy.expand_dims(coords, axis=1) - atom_coords
        norm_si_rJ = cupy.linalg.norm(si_rJ, axis=-1)
        si_rJ[:,i_grid_atom,:] = 0.0

        if surface_discretization_method.upper() == "SWIG":
            diJ = (norm_si_rJ - R_in_J) / R_sw_J
            diJ[:,i_grid_atom] = 1.0
            diJ[diJ < 1e-8] = 0.0
            si_rJ[diJ < 1e-8] = 0.0

            fiJ = switch_h(diJ)
            dfiJ = grad_switch_h(diJ)

            fiJK = fiJ[:, :, cupy.newaxis] * fiJ[:, cupy.newaxis, :]
            dfiJK = dfiJ[:, :, cupy.newaxis] * dfiJ[:, cupy.newaxis, :]
            R_sw_JK = R_sw_J[:, cupy.newaxis] * R_sw_J[cupy.newaxis, :]
            norm_si_rJK = norm_si_rJ[:, :, cupy.newaxis] * norm_si_rJ[:, cupy.newaxis, :]
            terms_size_ngrids_natm_natm = dfiJK / (fiJK * norm_si_rJK * R_sw_JK)
            si_rJK = si_rJ[:, :, cupy.newaxis, :, cupy.newaxis] * si_rJ[:, cupy.newaxis, :, cupy.newaxis, :]
            d2fiJK_offdiagonal = terms_size_ngrids_natm_natm[:, :, :, cupy.newaxis, cupy.newaxis] * si_rJK

            d2fiJ = gradgrad_switch_h(diJ)
            terms_size_ngrids_natm = d2fiJ / (norm_si_rJ**2 * R_sw_J) - dfiJ / (norm_si_rJ**3)
            si_rJJ = si_rJ[:, :, :, cupy.newaxis] * si_rJ[:, :, cupy.newaxis, :]
            d2fiJK_diagonal = cupy.einsum('qA,qAdD->qAdD', terms_size_ngrids_natm, si_rJJ)
            d2fiJK_diagonal += cupy.einsum('qA,dD->qAdD', dfiJ / norm_si_rJ, cupy.eye(3))
            d2fiJK_diagonal /= (fiJ * R_sw_J)[:, :, cupy.newaxis, cupy.newaxis]

            d2fiJK = d2fiJK_offdiagonal
            for i_atom in range(natom):
                d2fiJK[:, i_atom, i_atom, :, :] = d2fiJK_diagonal[:, i_atom, :, :]
        elif surface_discretization_method.upper() == "ISWIG":
            xi = charge_exp[p0:p1]
            erf_input_p = xi[:, None] * (R_J[None, :] + norm_si_rJ)
            erf_input_m = xi[:, None] * (R_J[None, :] - norm_si_rJ)
            fiJ = 1 - 0.5 * (erf(erf_input_p) + erf(erf_input_m))
            # fiJ[:,i_grid_atom] = 1.0
            dfiJ = 1/numpy.sqrt(numpy.pi) * xi[:, None] * (cupy.exp(-erf_input_m**2) - cupy.exp(-erf_input_p**2))
            ### It is necessary to zero out i \in I term in dfiJ, because the second term of d2fiJK_diagonal
            ### is not zeroed out otherwise.
            dfiJ[:,i_grid_atom] = 0

            fiJK = fiJ[:, :, cupy.newaxis] * fiJ[:, cupy.newaxis, :]
            dfiJK = dfiJ[:, :, cupy.newaxis] * dfiJ[:, cupy.newaxis, :]
            norm_si_rJK = norm_si_rJ[:, :, cupy.newaxis] * norm_si_rJ[:, cupy.newaxis, :]
            terms_size_ngrids_natm_natm = dfiJK / (fiJK * norm_si_rJK)
            si_rJK = si_rJ[:, :, cupy.newaxis, :, cupy.newaxis] * si_rJ[:, cupy.newaxis, :, cupy.newaxis, :]
            d2fiJK_offdiagonal = terms_size_ngrids_natm_natm[:, :, :, cupy.newaxis, cupy.newaxis] * si_rJK

            d2fiJ = 2.0/cupy.sqrt(cupy.pi) * (xi**2)[:, None] * (
                + erf_input_m * cupy.exp(-erf_input_m**2)
                + erf_input_p * cupy.exp(-erf_input_p**2)
            )
            # d2fiJ[:,i_grid_atom] = 0.0

            terms_size_ngrids_natm = d2fiJ / (norm_si_rJ**2) - dfiJ / (norm_si_rJ**3)
            si_rJJ = si_rJ[:, :, :, cupy.newaxis] * si_rJ[:, :, cupy.newaxis, :]
            d2fiJK_diagonal = cupy.einsum('qA,qAdD->qAdD', terms_size_ngrids_natm, si_rJJ)
            d2fiJK_diagonal += cupy.einsum('qA,dD->qAdD', dfiJ / norm_si_rJ, cupy.eye(3))
            d2fiJK_diagonal /= fiJ[:, :, cupy.newaxis, cupy.newaxis]

            d2fiJK = d2fiJK_offdiagonal
            for i_atom in range(natom):
                d2fiJK[:, i_atom, i_atom, :, :] = d2fiJK_diagonal[:, i_atom, :, :]
        else:
            raise NotImplementedError(f"surface_discretization_method = {surface_discretization_method} not recognized")

        Fi = switch_fun[p0:p1]
        Ai = area[p0:p1]

        d2F[p0:p1, :, :, :, :] += contract('q,qABdD->qABdD', Fi, d2fiJK)
        d2A[p0:p1, :, :, :, :] += contract('q,qABdD->qABdD', Ai, d2fiJK)

        d2fiJK_grid_atom_offdiagonal = -cupy.einsum('qABdD->qAdD', d2fiJK)
        d2F[p0:p1, i_grid_atom, :, :, :] = contract('q,qAdD->qAdD', Fi, d2fiJK_grid_atom_offdiagonal.transpose(0,1,3,2))
        d2F[p0:p1, :, i_grid_atom, :, :] = contract('q,qAdD->qAdD', Fi, d2fiJK_grid_atom_offdiagonal)
        d2A[p0:p1, i_grid_atom, :, :, :] = contract('q,qAdD->qAdD', Ai, d2fiJK_grid_atom_offdiagonal.transpose(0,1,3,2))
        d2A[p0:p1, :, i_grid_atom, :, :] = contract('q,qAdD->qAdD', Ai, d2fiJK_grid_atom_offdiagonal)

        d2fiJK_grid_atom_diagonal = -cupy.einsum('qAdD->qdD', d2fiJK_grid_atom_offdiagonal)
        d2F[p0:p1, i_grid_atom, i_grid_atom, :, :] = contract('q,qdD->qdD', Fi, d2fiJK_grid_atom_diagonal)
        d2A[p0:p1, i_grid_atom, i_grid_atom, :, :] = contract('q,qdD->qdD', Ai, d2fiJK_grid_atom_diagonal)

    d2F = d2F.transpose(1,2,3,4,0)
    d2A = d2A.transpose(1,2,3,4,0)
    return d2F, d2A
