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

Checkpoint long scientific calculations on durable shared storage. Save geometry
evaluations, optimizer-returned coordinates, numerical settings and ownership so
interruption or nonconvergence does not erase progress. Checkpoints are unqualified
until the normal convergence and physical QC pass. Keep oligomer construction,
optimization, energy evaluation and follow-up analysis in the versioned pipeline.

## Candidate results and slides

For a newly introduced candidate, screen the precursor and attachment route before
launching its full redox-property campaign. Record the reaction class, defined
conditions, precursor access, selectivity, redox-core retention and source/evidence
in `results/grafting/route_evidence.json`; run the prospective route gate in
`scripts/pipeline/screen_grafting_routes.py`. Use proceed/hold/exclude decisions.
Design future molecules as exact precursor/product pairs with the attachment
atom and new bond specified. Apply `docs/PRACTICAL_GRAFTING.md`; distinguish
exact-support demonstrations, molecular analogues and proposed reactions.
Do not equate a graph attachment, SA score or base-free reference deltaG with easy
grafting. Strongly unfavorable actual-condition routes should not enter the synthesis
shortlist. Preserve existing results. Explicitly scoped method controls may be
studied without implying synthesis eligibility.

Maintain `results/README.md` as the single human entry point and `results/INDEX.json`
as its machine-readable counterpart. Every new result category or campaign must
record its stable current outputs, status/QC, raw calculation location and source
snapshot or manifest. Update these pointers when publishing results, rather than
making readers search dated directories. Keep immutable evidence and historical
outputs at their original locations; label prepared, pending and verified results
explicitly. Check local index links before committing an update.

The confirmed synthesis route is post-functionalization: attach redox groups
onto an already formed polymer, such as a Merrifield-type polymer with reactive
benzyl sites (user correction, 2026-10-19). Existing 4-methylbenzyl structures are
local graft-site models; they do not include the full polymer environment.
Evaluate attachment feasibility, graft loading, swelling/accessibility, redox
retention and transport. The earlier vinyl-monomer-before-polymerization route
was superseded. Preserve its source-frozen calculations as a separate exploratory
branch; do not present vinyl polymerization compatibility as the current route's
main selection gate. Keep source-parent, capped-model, vinyl and polymer results
separate.
Keep a maintained collection of a few focused figures with stable filenames.
Separate computed descriptor leaders from synthesis priorities. Molecular and
redox-state diagrams must distinguish computed state assignments from proposed
ladders; illustrative resonance forms do not establish charge/spin localization.

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
