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

    def build(self, ng=None):
        if hasattr(self, '_solvent'):
            self.solvent = self._solvent
        solvent_descriptors = self.solvent_descriptors or solvent_db[self.solvent]
        if self.radii_table is None:
            radii_table = smd_radii(solvent_descriptors[2])
        else:
            radii_table = cupy.asnumpy(self.radii_table)
        logger.debug(self, 'radii_table %s', radii_table*radii.BOHR)
        mol = self.mol
        if ng is None:
            ng = self.sasa_ng

        if natm_without_ghost(mol) != mol.natm:
            raise RuntimeError('SMD does not support ghost atoms')

        self.surface = pcm.gen_surface(mol, rad=radii_table, ng=ng)
        self._intermediates = {}
        F, A = pcm.get_F_A(self.surface)
        D, S = pcm.get_D_S(self.surface, with_S=True, with_D=True)

        epsilon = self.eps or solvent_descriptors[5]
        f_epsilon = (epsilon - 1.0)/(epsilon + 1.0) if epsilon != float('inf') else 1.
        DA = D*A
        DAS = cupy.dot(DA, S)
        K = S - f_epsilon/(2.0*np.pi) * DAS
        K_LU, K_LU_pivot = lu_factor(K, overwrite_a = True, check_finite = False)
        K = None

        intermediates = {
            'S': cupy.asarray(S),
            'D': cupy.asarray(D),
            'A': cupy.asarray(A),
            'K_LU': cupy.asarray(K_LU),
            'K_LU_pivot': cupy.asarray(K_LU_pivot),
            'f_epsilon': f_epsilon
        }
        self._intermediates.update(intermediates)

        charge_exp  = self.surface['charge_exp']
        grid_coords = self.surface['grid_coords']
        atom_coords = mol.atom_coords(unit='B')
        atom_charges = mol.atom_charges()

        # Move this to GPU
        intopt = int3c1e.VHFOpt(mol)
        intopt.build(1e-14)
        self.intopt = intopt

        int2c2e = mol._add_suffix('int2c2e')
        fakemol_charge = gto.fakemol_for_charges(grid_coords.get(), expnt=charge_exp.get()**2)
        fakemol_nuc = gto.fakemol_for_charges(atom_coords)
        v_ng = gto.mole.intor_cross(int2c2e, fakemol_nuc, fakemol_charge)
        v_grids_n = np.dot(atom_charges, v_ng)
        self.v_grids_n = cupy.asarray(v_grids_n)
