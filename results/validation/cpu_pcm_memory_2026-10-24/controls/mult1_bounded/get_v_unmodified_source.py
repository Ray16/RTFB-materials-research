    def _get_v(self, dms):
        '''
        return electrostatic potential on surface
        '''
        mol = self.mol
        nao = dms.shape[-1]
        grid_coords = self.surface['grid_coords']
        exponents   = self.surface['charge_exp']
        ngrids = grid_coords.shape[0]
        nset = dms.shape[0]
        v_grids_e = numpy.empty([nset, ngrids])
        max_memory = self.max_memory - lib.current_memory()[0]
        blksize = int(max(max_memory*.9e6/8/nao**2, 400))
        int3c2e = mol._add_suffix('int3c2e')
        cintopt = gto.moleintor.make_cintopt(mol._atm, mol._bas, mol._env, int3c2e)
        for p0, p1 in lib.prange(0, ngrids, blksize):
            fakemol = gto.fakemol_for_charges(grid_coords[p0:p1], expnt=exponents[p0:p1]**2)
            fakemol.cart = mol.cart
            v_nj = df.incore.aux_e2(mol, fakemol, intor=int3c2e, aosym='s1', cintopt=cintopt)
            for i in range(nset):
                v_grids_e[i,p0:p1] = numpy.einsum('ijL,ij->L',v_nj, dms[i])

        return v_grids_e
