# Exact site-3 fixed-density PCM response qualification

The [declared test](plan.json) uses the checked 115-atom site-3 geometry,
R2SCAN-D4/def2-SVP and SMD acetonitrile. A baseline SCF supplies a fixed AO density
and MO coefficients. The analytic electrostatic solvent-potential derivative is
compared with finite differences at atom 0/y and the first oxygen/z, at 0.0005
and 0.001 bohr. The displaced evaluations do not repeat SCF or change the density.
This tests the response component that exhausted memory in the full native job.

[The completed record](status.json) and [independent raw-matrix replay](independent_review.json)
pass all four one-percent comparisons; the largest relative error is 4.143e-6.
The analytic response has shape 115 × 3 × 1150 × 231. Observed per-process GPU
peak is 31,224 MiB on a 32 GiB V100. [Small full spin controls](../pcm_first_response_workspace_v2_2026-10-23/README.md)
separately establish original-versus-modified Hessian equivalence.

This is a fixed-density electrostatic response/resource test. It does not
establish a native Hessian minimum, total free energy, chemical accuracy or
candidate eligibility. The full native pilot has separate resource and physical
criteria; no native result is copied from this test.

Raw matrices, geometry, baseline SCF checkpoint and execution sources remain at
`calcs/validation/pmdi_site3_PCM_response_v2_2026-10-23/`. [plan.json](plan.json)
and [status.json](status.json) hash-bind these files. The large response tensor is
kept on shared storage. Reproduce the independent matrix arithmetic with
[scripts/validation/review_pcm_first_response.py](../../../scripts/validation/review_pcm_first_response.py).
The native [preparation](../../polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_3_native_first_response/README.md)
requires these exact-geometry receipts and retains unchanged Hessian gates.
