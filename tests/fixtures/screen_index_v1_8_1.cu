/*
 * Copyright 2021-2024 The PySCF Developers. All Rights Reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

__global__
static void _screen_index(int8_t *non0shl_mask, double log_cutoff,
                          double *coords, int ngrids, int block_size,
                          int *atm, int natm, int *bas, int nbas, double *env)
{
    int grid_block_id = blockIdx.x;
    int grid_start = grid_block_id * block_size;
    int grid_stop = min(grid_start+block_size, ngrids);
    int shl_block_id = blockIdx.y;
    int thread_id = threadIdx.x;
    int ish = shl_block_id * blockDim.x + thread_id;
    if (ish >= nbas) {
        ish = 0;
    }

    int atm_id = bas[ish*BAS_SLOTS+ATOM_OF];
    int ang = bas[ish*BAS_SLOTS+ANG_OF];
    int nprim = bas[ish*BAS_SLOTS+NPRIM_OF];
    double *exps = env + bas[ish*BAS_SLOTS+PTR_EXP];
    double *coeffs = env + bas[ish*BAS_SLOTS+PTR_COEFF];
    double *ri = env + atm[atm_id*ATM_SLOTS+PTR_COORD];
    double atom_x = ri[0];
    double atom_y = ri[1];
    double atom_z = ri[2];

    __shared__ double gridx_cache[NG_PER_BLOCK*3];
    double *gridy_cache = gridx_cache + NG_PER_BLOCK;
    double *gridz_cache = gridy_cache + NG_PER_BLOCK;

    int8_t is_large = 0;
    for (int grid0 = grid_start; grid0 < grid_stop; grid0 += NG_PER_BLOCK) {
        if (grid0 + thread_id < ngrids) {
            gridx_cache[thread_id] = coords[         grid0+thread_id];
            gridy_cache[thread_id] = coords[ngrids  +grid0+thread_id];
            gridz_cache[thread_id] = coords[ngrids*2+grid0+thread_id];
        }
        __syncthreads();
        // check if any GTO values on grids are larger than threshold
        if (is_large) {
            continue;
        }
        int ng_in_tile = min(NG_PER_BLOCK, ngrids-grid0);
        for (int grid_id = 0; grid_id < ng_in_tile; ++grid_id) {
            double rx = gridx_cache[grid_id] - atom_x;
            double ry = gridy_cache[grid_id] - atom_y;
            double rz = gridz_cache[grid_id] - atom_z;
            double rr = rx * rx + ry * ry + rz * rz + 1e-300;
            double gto_sup = 1e-300;
            for (int ip = 0; ip < nprim; ++ip) {
                gto_sup += coeffs[ip] * exp(-exps[ip] * rr);
            }
            is_large |= log(fabs(gto_sup)) + ang*log(rr)/2 > log_cutoff;
        }
    }
    if (shl_block_id * blockDim.x + thread_id < nbas) {
        non0shl_mask[grid_block_id*nbas + ish] = is_large;
    }
}

