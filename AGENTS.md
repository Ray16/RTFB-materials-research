# Project instructions

This is the cleaned redox library at `/nfs/lambda_stor_01/homes/rzhu/0_redox_refactor`.
The manuscript, slides, dated campaign outputs and bulk calculations remain in
`/nfs/lambda_stor_01/homes/rzhu/0_redox`. Do not copy those trees back here.

GPU work uses `gpu_reserve`. Do not set `CUDA_VISIBLE_DEVICES` directly.
Cap `OMP_NUM_THREADS`, `MKL_NUM_THREADS` and `OPENBLAS_NUM_THREADS` on this
shared machine. Do not launch DFT or Hessian jobs unless explicitly asked.

`src/redox/qm/qrrho.py`, `src/redox/qm/xtb_kernel.py` and
`src/redox/quality/sampling_snapshot.py` are certificate inputs. Do not reformat
them. A change that is not the single admitted xTB library-path transformation
invalidates the projected-thermal certificate.

Keep half-reactions and reorganization components as separate quantities.
Do not introduce a weighted-sum score. Use real chemical names.
