# Expanded conditional-polymer conformer searches

The first wider panel covers six occupancy graphs on a five-repeat phthalimide
control at two declared backbone configurations. These are conditional scenarios,
not an experimentally measured atactic sequence or loading distribution.

The ten missing graph/configuration inputs were built with 32 ETKDGv3 embeds,
seed 20261023, serial MMFF optimization and exact graph/stereochemistry checks.
The [input plan](inputs/plan.json) records their sources. Two additional full-
parallel empty-backbone repeats reuse the prior checked initial geometries;
[remaining_wave_plan.json](remaining_wave_plan.json) describes these explicitly.
The earlier two-core empty-backbone pilots remain separate, preserved searches.
No statistical RNG independence is assumed for the repeats.

The new searches use CREST 3.0.2 iMTD-GC, GFN2/ALPB(acetonitrile), an 8 kcal/mol
energy window, topology checks and the program's normal search settings. Actual
settings, trajectories and calculator inputs are retained on shared storage.
A dry run records the SHAKE/all-bond, 5-fs and trajectory settings. Sampling
trajectories are not treated as equilibrium populations or DFT free energies.

## Actual execution

One 14-core worker on lambda2 and eleven 14-core workers on lambda13 entered
metadynamics. OpenMP is bounded by the assigned 14 cores; BLAS, MKL, NumExpr and
external xTB pools are one thread. The previous serial-reference pilots' nested
BLAS warnings motivated this execution change; their code and running processes
were not altered. New start checks found zero OpenBLAS warnings.

[Local pilot verification](pilot_start_verification.json) and
[lambda13 verification](start_verification_lambda13.json) check actual parent and
predictor PID/start identity, inherited affinity, six thread caps and disjoint
CPU sets (154 cores on lambda13, 14 on lambda2). Node preflight and the resource
inventory accompany them. CREST is a CPU calculation; no GPU was reserved for
these searches. Raw outputs: `calcs/polymer_models/2026-10-23_conformer_revision/`.

## Remaining scientific steps

A completed low-level search still needs conformer graph/stereochemistry QC,
DFT reoptimization/re-ranking, duplicate-minimum handling, native thermochemistry
and a defined truncated-ensemble treatment. The target is 10 then 20 unique DFT
starts per graph/configuration, expanding when minimum/ensemble changes exceed
2 kJ/mol. This is a within-search/batch convergence diagnostic, not proof of a
global minimum. New DFT follow-up jobs have not yet been submitted by this plan.
The two configurations stay separate; no experimental population weights are
invented. Repeated searches must not double-count the same conformer in a
partition function.

Producers: `prepare_conformer_revision.py`, `sample_polymer_conformers.py` and
`launch_conformer_revision.py` under `scripts/pipeline/`. The worker preserves
long calculations, validates every output conformer, owns its child process group
and refuses incomplete XYZ frames. The shared CPU allocator now respects existing
scientific workers with more than eight assigned cores.

The two earlier two-core pilots were subsequently retired after their repetitive
warning logs reached about 108 GiB. Parent/predictor ownership was verified and
normal termination handlers preserved all raw trajectories and checkpoints;
[the retirement record](legacy_pilot_retirement.json) confirms both handles ended.
The twelve corrected source searches remain separate. Their final outputs feed
the [DFT/native continuation](../../polymer_models/campaigns/2026-10-23_crest_refinement/README.md),
with missing quantities retained explicitly until normal physical admission.
