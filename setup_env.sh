#!/usr/bin/env bash
# Bootstrap the `redox` conda environment; validated native-QM versions are pinned.
# See docs/DEVELOPMENT.md for the recorded working environment and validation gates.
#
#   ./setup_env.sh [ENV_NAME]
#
# IMPORTANT — CUDA build: torch must match your GPU driver's CUDA version. Check it with
#   nvidia-smi        # top-right "CUDA Version:"
# then set CUDA_TAG to a compatible PyTorch 2.8 wheel tag, or
# "cpu" for a CPU-only box. Installing a torch newer than your driver -> torch.cuda is False.
#
#   CUDA_TAG=cu128 ./setup_env.sh          # default on GPU hosts
#   CUDA_TAG=cpu   ./setup_env.sh          # CPU-only
set -euo pipefail

ENV_NAME="${1:-redox}"
PY_VER="${PY_VER:-3.11}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Auto-detect GPU vs CPU unless CUDA_TAG is set explicitly.
#   GPU present -> cu128 (requires a CUDA-12.8-capable driver).
#   No nvidia-smi -> cpu.
if [ -z "${CUDA_TAG:-}" ]; then
  if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
    CUDA_TAG="cu128"
    echo ">> GPU detected -> CUDA_TAG=$CUDA_TAG (requires a CUDA-12.8-capable driver)"
  else
    CUDA_TAG="cpu"
    echo ">> no GPU detected -> CUDA_TAG=cpu"
  fi
fi

# Prefer mamba if available (faster), else conda.
CONDA="$(command -v mamba || command -v conda)"
echo ">> using $CONDA ; env=$ENV_NAME python=$PY_VER cuda=$CUDA_TAG"

# 1. Create the environment.
"$CONDA" create -n "$ENV_NAME" "python=$PY_VER" -y

# 2. Install a CUDA-matched torch FIRST, from the CUDA-specific index.
if [ "$CUDA_TAG" = "cpu" ]; then
  TORCH_INDEX="https://download.pytorch.org/whl/cpu"
else
  TORCH_INDEX="https://download.pytorch.org/whl/${CUDA_TAG}"
fi
conda run -n "$ENV_NAME" pip install "torch==2.8.0" --index-url "$TORCH_INDEX"

# 3. Install the rest. torch is pinned in requirements.txt to the version above, so this
#    step leaves the CUDA build from step 2 untouched.
conda run -n "$ENV_NAME" pip install -r "$HERE/requirements.txt"

# 3a. xtb provides legacy GFN2 thermal corrections and qRRHO integration for the
#     independently derivative-checked native solution Hessian. These are distinct
#     thermal models; the selected policy determines which record can be admitted.
conda install -n "$ENV_NAME" -c conda-forge xtb -y

# 3b. GPU DFT backend (gpu4pyscf) — CUDA-tagged wheel, like torch. Only for cu12x GPUs;
#     hugely faster DFT+SMD than CPU pyscf. Skipped on CPU-only boxes.
case "$CUDA_TAG" in
  cu12*)
    conda run -n "$ENV_NAME" pip install gpu4pyscf-cuda12x==1.8.1 cutensor-cu12==2.7.0
    # The gpu4pyscf install can drag in CUDA-13 nvidia-* wheels (unsuffixed metapackages),
    # whose cu13/ headers+libs shadow the cu12 ones and break nvrtc on pre-Ada GPUs (e.g.
    # V100/sm_70 — CUDA 13 dropped Volta). Remove any nvidia-* that is NOT the -cu12 build
    # so cuda-pathfinder resolves the cu12 runtime cleanly.
    conda run -n "$ENV_NAME" python - <<'PYCLEAN'
import subprocess, sys
out = subprocess.check_output([sys.executable,"-m","pip","list","--format=freeze"]).decode()
bad = [l.split("==")[0] for l in out.splitlines()
       if l.startswith("nvidia-") and "-cu12" not in l.split("==")[0]]
if bad:
    print(">> removing stray non-cu12 nvidia packages:", bad)
    subprocess.check_call([sys.executable,"-m","pip","uninstall","-y",*bad])
import site, shutil, os
cu13 = os.path.join(site.getsitepackages()[0], "nvidia", "cu13")
shutil.rmtree(cu13, ignore_errors=True)
PYCLEAN
    ;;
  cpu)
    echo ">> CPU-only: skipping gpu4pyscf (DFT will use CPU pyscf)" ;;
  *)
    echo ">> CUDA_TAG=$CUDA_TAG: install a matching gpu4pyscf wheel manually if desired" ;;
esac

# 4. Verify (GPU visibility + imports).
conda run -n "$ENV_NAME" python "$HERE/scripts/check_env.py"

cat <<EOF

Done. Activate with:
  source \$(conda info --base)/etc/profile.d/conda.sh && conda activate $ENV_NAME

UMA weights are gated on HuggingFace. Request access at
  https://huggingface.co/facebook/UMA
then authenticate (keeps the token out of shell history):
  huggingface-cli login
EOF
