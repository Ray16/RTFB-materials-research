# Project instructions

## Parallel scientific calculations

Use as many currently free, suitable GPUs as the authorized workload can use.
Do not impose an arbitrary small worker cap when more suitable GPUs are available.
This is the user's standing preference, reaffirmed on 2026-10-11.

Read the `lambda-fleet` skill before submitting, launching, restarting or resuming
fleet calculations. The canonical skill is
`/nfs/lambda_stor_01/homes/rzhu/.codex-personal/skills/lambda-fleet/SKILL.md`.
Use the system reservation policy and `gpu_reserve` for every GPU calculation.
Occupied GPUs are unavailable even when their utilization is low.

Suitability depends on measured per-job memory, architecture, dependencies and
node health. The current full molecular DFT/Hessian workflow requires the proven
32 GB class; do not send it to 16 GB cards without separate memory validation.
Use smaller cards for independently validated workloads when appropriate.
GPU memory is not pooled across independent workers.

Preflight every node, assign disjoint CPU sets, cap native thread pools, check
CPU load/RAM, and preserve resumable outputs and ownership records. Expand the
fleet within these constraints. Keep scientific protocols and validation gates
consistent across nodes; report actual starts, completed records and unresolved
physical validation separately.

## Candidate results and slides

Collect verified results for candidates that finish, with pending quantities
blank and field-specific QC. Keep source-parent measurements separate from the
properties of exact grafted structures. Update `slides/redox_workflow.tex` and
its checked PDF as verified candidate results arrive; preserve earlier sources
and label historical results whose reference or thermal treatment is superseded.

## Git checkpoints

Commit and push each completed, verified milestone to the existing GitHub remote.
Include code, configuration, relevant tests, small result artifacts and a record
of checks and remaining scientific limits. Preserve unrelated staged changes.
Keep bulk datasets, model weights and large calculation outputs on shared storage;
record their locations and hashes with the checked results. Report the pushed commit.
