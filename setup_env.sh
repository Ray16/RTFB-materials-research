#!/usr/bin/env bash
# Reproducible setup for the `redox` conda environment.
#
#   ./setup_env.sh [ENV_NAME]
#
# IMPORTANT — CUDA build: torch must match your GPU driver's CUDA version. Check it with
#   nvidia-smi        # top-right "CUDA Version:", e.g. 12.4
# then set CUDA_TAG to the matching PyTorch wheel tag (cu121 / cu124 / cu126 / ...), or
# "cpu" for a CPU-only box. Installing a torch newer than your driver -> torch.cuda is False.
#
#   CUDA_TAG=cu124 ./setup_env.sh          # default
#   CUDA_TAG=cpu   ./setup_env.sh          # CPU-only
set -euo pipefail

ENV_NAME="${1:-redox}"
PY_VER="${PY_VER:-3.11}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Auto-detect GPU vs CPU unless CUDA_TAG is set explicitly.
#   GPU present -> cu124 (works with driver CUDA >= 12.4; override for older drivers).
#   No nvidia-smi -> cpu.
if [ -z "${CUDA_TAG:-}" ]; then
  if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
    CUDA_TAG="cu128"
    echo ">> GPU detected -> CUDA_TAG=$CUDA_TAG (override if your driver's CUDA < 12.4)"
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
conda run -n "$ENV_NAME" pip install "torch==2.6.0" --index-url "$TORCH_INDEX"

# 3. Install the rest. torch is pinned in requirements.txt to the version above, so this
#    step leaves the CUDA build from step 2 untouched.
conda run -n "$ENV_NAME" pip install -r "$HERE/requirements.txt"

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
