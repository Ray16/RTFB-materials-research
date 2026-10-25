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

def get_dqsym_dx_fix_K_R(pcmobj, dm, atmlst, intopt_derivative):
    dV_on_charge_dx = get_dvgrids(pcmobj, dm, atmlst, intopt_derivative)

    f_epsilon = pcmobj._intermediates['f_epsilon']
    if pcmobj.if_method_in_CPCM_category:
        R = -f_epsilon
        R_dVdx = R * dV_on_charge_dx
    else:
        R_dVdx = _apply_R_batch(pcmobj, dV_on_charge_dx, False)

    K_1_R_dVdx = einsum_ij_Adj_Adi_inverseK(pcmobj, R_dVdx)
    K_1T_dVdx = einsum_ij_Adj_Adi_inverseK(pcmobj, dV_on_charge_dx, K_transpose = True)

    if pcmobj.if_method_in_CPCM_category:
        RT_K_1T_dVdx = R * K_1T_dVdx
    else:
        RT_K_1T_dVdx = _apply_R_batch(pcmobj, K_1T_dVdx, True)

    dqdx_fix_K_R = 0.5 * (K_1_R_dVdx + RT_K_1T_dVdx)

    return dqdx_fix_K_R
