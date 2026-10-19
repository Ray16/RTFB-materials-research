# Development guide

The pipeline builds explicit molecular states, computes their geometries and
energies, applies field-specific physical QC, and derives battery screening
properties. Reusable implementations live in `src/redox/`. Start new orchestration
in `redox.workflows`; keep command wrappers small.

## Working environment

Use the existing environment without rerunning setup:

```bash
source /nfs/lambda_stor_01/homes/rzhu/miniforge3/etc/profile.d/conda.sh
conda activate redox
export PYTHONPATH="$PWD/src"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export FAIRCHEM_CACHE_DIR=/nfs/lambda_stor_01/homes/rzhu/.cache/fairchem
```

The reviewed NFS environment has Python 3.11, torch 2.8.0+cu128,
fairchem-core 2.21.0, PySCF 2.14.0, GPU4PySCF 1.8.1, cuTENSOR 2.7.0,
ASE 3.29.0, RDKit 2026.3.5, NumPy 2.2.6, SciPy 1.17.1, pandas 3.0.5,
openpyxl 3.1.5 and pytest 9.1.1. This records the working environment;
it is not a portable lock file. `setup_env.sh` bootstraps a new environment.
Native thermal records additionally require the library versions and derivative
controls specified by `config/native_thermal_validation.json`. Installing matching
packages alone does not prove a new node is ready for calculations.

## Supported tasks

| Task | Command or implementation |
| --- | --- |
| Run software regression checks | `python -m pytest -q` |
| Validate a candidate definition without writing structures | `python -m redox.build.candidates --config config/custom_candidates.py --dry-run` |
| Inspect every declared state | `python -m redox.validation.readiness` |
| Resolve original repair geometries without rewriting old evidence | `python -m redox status --include-lineage --out results/reviews/LINEAGE.csv` |
| Inspect current energy coverage | `python -m redox.qm.sp --audit` |
| Inspect authoritative GPU availability | `python -m redox.workflows.resources --all` |
| Plan all suitable local DFT/UMA workers | `python -m redox.workflows.launch --module dft --out logs/runs/NAME` |
| Run a reviewed full molecular QM command | `python -m redox.workflows.gpu --gpu INDEX --out logs/runs/NAME -- COMMAND` |
| Rebuild and audit tables using existing records | `python -m redox.workflows.tables --out results/reviews/NAME` |
| Prepare a candidate execution plan | `python scripts/pipeline/run_motifs.py --config config/custom_candidates.py --only ID --out results/campaigns/NAME` |
| Prepare matched-comparison hypotheses | `python -m redox expand --out results/discovery/NAME` (`--register` appends definitions) |
| Enumerate viologen/PMDI/NDI exact graphs | `python -m redox search --out results/discovery/NAME` |
| Select four matched contrasts per core | `python -m redox search --select-from POOL_DIR --out NEW_DIR` |
| Select matched OMe/Me ring contrasts | `python -m redox search --select-from POOL_DIR --panel donor_rings --out NEW_DIR` |
| Audit downloaded source quality flags | `python -m redox source-audit --source RAW_DIR --out results/discovery/NAME` |
| Revalidate historical SP adoption | `python -m redox.workflows.revalidate_sp prepare --out RUN`; `collect --run RUN --promote` preserves originals before an agreeing replacement |
| Collect verified candidate snapshots | `python scripts/pipeline/collect_verified_candidates.py --out results/reviews/NAME` (original cohort plus registered additional candidates) |
| Recompute corrected thermochemistry from existing qualified Hessians | `python -m redox thermal-postprocess --only ID --threads 2 --out results/validation/NEW_RUN` (CPU only) |
| Update candidate slides from a verified feed | `scripts/pipeline/update_verified_slides.py` |

`run_motifs.py` plans by default. Execution requires `--execute`, the selected GPU
and node suitability checks under the current fleet instructions. Worker counts
must follow suitable free resources and the authorized workload. Old dated plans
with fixed one/two-worker limits do not override `AGENTS.md`.

The explicit shared `FAIRCHEM_CACHE_DIR` matters across nodes: their home-directory
caches can differ even when Python itself is on NFS. The UMA cache checker reads
local checkpoint bytes without downloading. Verify the expected model hash before
launch; a missing node-local checkpoint is an environment failure, not a failed
molecular calculation. `--only` in the builder accepts exact comma-separated IDs.

Repair-source paths in older `result.json` files can refer to the replaced
canonical geometry. `status --include-lineage` resolves an original geometry only
when its bytes match the recorded SHA-256, using the explicitly named archive.
It never guesses by folder age or substitutes the current geometry. This view
preserves frozen record hashes. New repairs write archive pointers and capture
their promotion helper's source at import; their active SPs retain separate
geometry-addressed provenance. Lineage checks supplement per-state physical QC
and do not establish experimental accuracy.

`workflows.launch` plans by default; `--execute` submits every suitable free local
device. It checks total CPU/RAM requirements and assigns disjoint CPU sets before
starting resource controllers. Each controller rechecks resources under the shared
allocation lock, binds the worker, obtains a system GPU reservation and probes GPU
dependencies inside that reservation. These full-QM launch paths require the proven
32-GB-class devices. A smaller-card route requires separate workload validation.
For cross-node static shards, specify one global shard count and disjoint starting
indices for each node. Process exit status remains separate from scientific QC.

The compatibility command `python scripts/pipeline/review_tables.py --out ...`
uses the same table workflow. Each stage runs in a fresh child interpreter with
one native thread, a log and an explicit destination. A failed stage stops the
chain and records its name/exit status; an existing destination is never replaced.
The final independent audit determines arithmetic agreement. No review is
automatically published, and missing physical quantities remain incomplete.

## Source and evidence ownership

| Location | Ownership rule |
| --- | --- |
| `config/` | Explicit chemical definitions, scientific protocols and validation certificates |
| `library/manifest.csv` | All three builders use `build.manifest.locked_manifest` through read/build/merge; publish with `write_rows` |
| `calcs/dft/ID/STATE/result.json` and `opt.xyz` | Geometry evidence; preserve prior results when repairing |
| `calcs/dft/ID/STATE/sp/HASH.json` | Protocol-addressed energies; read through `core.protocol.load_record` |
| `calcs/dft/ID/STATE/native_thermal/` | Geometry-addressed thermal evidence; admit through the selected thermal policy |
| `calcs/dft/ID/STATE/projected_thermal/` | Corrected internal-mode partition, source hashes and isolated CPU-kernel receipts; raw Hessians remain preserved |
| `results/reviews/NAME/` | Isolated derived tables, stage logs and review status |
| `results/validation/` and dated campaign directories | Validation evidence and execution history; inspect manifests before reuse |
| `results/*.csv` and figures | Canonical derived products; update through checked publication paths |
| `slides/` | LaTeX source, generated verified includes and checked PDF |

Some library/calculation outputs are tracked, some are ignored, and substantial
new evidence is currently untracked. Check `git status` and `.gitignore`; never
assume a directory is disposable because it is generated. Paths and file hashes
are referenced by validation certificates and published snapshots. Do not move
them as a cosmetic cleanup.

Use `core.io.atomic_write_json` for new mutable records and status files. Atomic
replacement prevents a reader seeing a partial file; exclusive ownership or a
domain lock is still required when multiple processes might merge/write it.

## Scientific contracts

All charge states and Fc use `core.protocol.ACTIVE_SP`. Analysis reads the current
addressed energy; missing records do not fall back to `result.json`. New adoption
also requires the producer's `energy_protocol` and `energy_geom_sha1`, charge,
multiplicity, finite energies and convergence. Historical source results without
that proof require computation if their addressed energy is missing; existing
valid addressed records continue to load.

Free energy combines an admitted active energy and the explicitly selected
thermal definition. Candidate and Fc definitions must match. Exact graft
properties remain separate from experimental parent measurements. Gas Nelsen,
solution reorganization, outer-sphere approximations, nominal capacity and
accessible capacity have distinct meanings; preserve those labels.

Native Hessian execution evidence is checked separately from the physical level
of theory. The certified memory budgets and 1e-5/1e-6-Bohr internal CDS derivative
settings retain distinct receipts; a derivative variant needs its matching sealed
certificate and complete per-geometry contact checks. The default calculation
path is unchanged. [The qualification history](../results/validation/cds_finer_step_2026-10-18/README.md)
and [reviewed minimum correction](../results/validation/pmdi_trioxy_minimum_2026-10-18/README.md)
show how isolated calculations reach admission without changing physical gates.

Native primary free energies require `qm.projected_thermal` records backed by
the sealed `config/projected_thermal_validation.json` qualification. `qm.qrrho`
projects rigid modes using mass-weighted coordinates, retaining genuine soft
internal modes. `qm.xtb_kernel` runs the installed xTB thermodynamic kernel in a
dedicated CPU subprocess to avoid a PySCF/D4 shared-library conflict. The audit's
`native_partition.py` independently reconstructs the internal spectrum and thermal
equations. A raw native Hessian record alone is insufficient for primary admission.

Geometry optimization now passes an explicit full grid-response gradient scanner
to the PySCF/geomeTRIC kernel and records both geometry and SCF convergence.
Strict calls reject unconverged output before replacing an existing geometry;
explicit non-strict calls retain `optimized=False` when convergence fails.
Historical geometry metadata is not retroactively relabeled. Primary native
properties continue to require their independent stationarity and minimum checks.

Generic native minimum repairs prepare a complete replacement before changing
canonical files. Promotion includes `sp`, `native_thermal` and `projected_thermal`,
retains old addressed history, checks admission before and after installation and
rolls back failed installs. Checkpoint resumes rerun the bounded optimizer and
record convergence. Raw cross-point history remains preserved and its geometry
addresses prevent stale use after a minimum changes.

`scripts/validation/audit/recompute_axes.py` independently reimplements the
arithmetic and thermal admission. Keep this intentional duplication independent
of production property/ranking helpers. Reusing a production formula in its own
audit would weaken the test.

## Making changes

Read `AGENTS.md`, inspect the existing working tree and preserve unrelated work.
Change one ownership boundary at a time. Keep existing CLI paths working through
wrappers when moving reusable logic. Add regression tests for a real failure
mode, then run the relevant tests and the complete CPU suite. Geometry changes,
new QM methods and new candidate campaigns need their own physical validation
and resource preparation; a passing refactor suite does not establish those.

Spin publication predicates live in `redox.validation.spin_evidence`. The candidate
collector and scorecard both require completed, current source-bound helper
reports before treating alternate-spin energies as confidence evidence. Unresolved
helpers are labelled `DFT_check_required`, which blocks primary Pareto ranking.
A verified electronic gap close to zero still needs native spin thermochemistry;
matching near-degenerate native free energies require ensemble integration.
Numerical surface calculations remain distinct from primary spin admission.

## Equilibrium spin selections

`validation.spin_ensemble.prepare_policy(state_dir, helper_dir)` reads two existing,
qualified native minima without launching calculations. `select_policy` applies a
prepared policy under molecule/native/spin writer locks only when both the primary
and independent readers agree. The selection file binds current component inputs;
changed or invalid dependencies withhold the free energy rather than reverting to
one spin surface. Original helper seed identity is preserved separately from the
current selected minimum. Both component native terms already include spin degeneracy.

Equilibrium potentials do not qualify single-surface kinetics. Reorganization
values remain flagged diagnostics for selected mixtures, and integrity must pass
for every component pair before potential-window capacity is counted. See the
[implementation and evidence](../results/validation/spin_ensemble_integration_2026-10-19/README.md).

Publication dependency discovery is shared in `validation.provenance.spin_decision_paths`;
collectors and overlays must use that inventory when adding a new state policy.
Watch deadlines use `workflows.publication.watch_deadline`: the default closes at
23:55 Chicago on the current day, while `--max-hours` uses elapsed hours across
midnight/DST. The collector additionally honors its selected queue's terminal state.
CSV and workbook snapshots retain the union of optional QC fields, including
single-surface kinetic flags and explicit ensemble correction/population columns.
