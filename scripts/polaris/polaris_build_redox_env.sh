#!/usr/bin/env bash
# Build a minimal gpu4pyscf DFT env on Polaris for the reorg recompute.
# Needs only: pyscf + gpu4pyscf (cu12x, matches Polaris cuda/12.9 on A100/sm_80) + rdkit + ase
# + geomeTRIC + pandas/scipy. No torch/fairchem/UMA (we warm-start the anion from the neutral
# SMD geometry) and no xtb (no thermal needed for inner-sphere lambda).
set -euo pipefail

CONDA=/eagle/FoundEpidem/ray/miniforge3/etc/profile.d/conda.sh
ENV=/grand/FRAME-IDP/rzhu/env/redox
PY=3.11

# /eagle/FoundEpidem (miniforge3 home) is over quota -> redirect ALL caches/tmp to /grand,
# which has the FRAME-IDP allocation and space. Otherwise conda's repodata/pkgs cache and
# pip's download cache hit "Disk quota exceeded".
export CONDA_PKGS_DIRS=/grand/FRAME-IDP/rzhu/conda_pkgs
export PIP_CACHE_DIR=/grand/FRAME-IDP/rzhu/pip_cache
export TMPDIR=/grand/FRAME-IDP/rzhu/tmp
mkdir -p "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR" "$TMPDIR"

echo ">> $(date) building $ENV (python $PY); caches -> /grand"
source "$CONDA"
conda create -y -p "$ENV" "python=$PY"

echo ">> pip installing DFT stack"
conda run -p "$ENV" pip install --no-input \
    pyscf gpu4pyscf-cuda12x cutensor-cu12 rdkit ase geometric pandas scipy numpy

# gpu4pyscf-cuda12x + cupy-cuda12x do NOT pull the CUDA runtime libs on their own, and Polaris
# has no system CUDA on PATH by default, so gpu4pyscf import fails (libcublasLt.so.12 /
# libnvJitLink.so.12 not found). Install the cu12 runtime wheels explicitly to make the env
# self-contained on A100 (sm_80). NOTE: the lambda "remove non-cu12 nvidia wheels" step is a
# Volta/sm_70-only workaround (CUDA 13 dropped Volta) — do NOT run it here; on A100 it strips
# the runtime and breaks the import.
echo ">> installing cu12 CUDA runtime wheels (self-contained)"
conda run -p "$ENV" pip install --no-input \
    nvidia-cuda-runtime-cu12 nvidia-cuda-nvrtc-cu12 nvidia-cublas-cu12 nvidia-cufft-cu12 \
    nvidia-curand-cu12 nvidia-cusolver-cu12 nvidia-cusparse-cu12 nvidia-nvjitlink-cu12 \
    nvidia-cuda-cupti-cu12 nvidia-nccl-cu12

echo ">> import check (CPU-only parts; gpu4pyscf import needs a GPU node)"
conda run -p "$ENV" python -c "import pyscf, rdkit, ase, geometric, pandas, scipy; print('pyscf', pyscf.__version__, 'rdkit', rdkit.__version__)"
echo ">> $(date) BUILD DONE (gpu4pyscf import + a GPU test still to run on a compute node)"
