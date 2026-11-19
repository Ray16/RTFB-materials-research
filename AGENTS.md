# Project instructions

This repo at `/nfs/lambda_stor_01/homes/rzhu/0_redox_refactor` is the single home
for all active work: the redox library (`src/redox/`), campaigns (`results/`),
and the manuscript (`manuscript/`). The historical raw-record archive was missing
at the last check (2026-11-19). Do not assume it exists or replace missing records
with unrelated calculations. Use `REDOX_DATA_ROOT` only after verifying that the
selected data tree contains the requested campaign. Paper rendering uses the
committed inputs in `data/manuscript/` and `manuscript/evidence/`.

GPU work uses `gpu_reserve`. Do not set `CUDA_VISIBLE_DEVICES` directly.
Cap `OMP_NUM_THREADS`, `MKL_NUM_THREADS` and `OPENBLAS_NUM_THREADS` on this
shared machine. Do not launch DFT or Hessian jobs unless explicitly asked.

`src/redox/qm/qrrho.py`, `src/redox/qm/xtb_kernel.py` and
`src/redox/quality/sampling_snapshot.py` are certificate inputs. Do not reformat
them. A change that is not the single admitted xTB library-path transformation
invalidates the projected-thermal certificate.

Keep half-reactions and reorganization components as separate quantities.
Do not introduce a weighted-sum score. Use real chemical names.
