# Refactor prompt — clean reconstruction of the redox pipeline

> Paste the block below into a fresh session (ideally in plan mode) to execute the refactor.

---

## Task: Port the redox screening pipeline into a clean architecture

You are refactoring a mature but messy research codebase into a clean, well-structured
reimplementation. This is a **port into a clean architecture**, NOT a blank-slate rewrite and
NOT an in-place shuffle: design a sound target structure, then carry the *validated scientific
logic* over module-by-module, rewriting the glue, interfaces and layering as you go. Every
scientific number the current code produces must be reproduced **bit-for-bit** by the new tree.

### Paths
- **Source (read-only reference):** `/nfs/lambda_stor_01/homes/rzhu/0_redox`
- **Destination (build here):** `/nfs/lambda_stor_01/homes/rzhu/0_redox_refactored`
- Do not modify the source tree except to read it. Initialize the destination as a fresh git repo.

### What this project is
Computational screening of redox-active motifs for **post-functionalization of Merrifield-type
polymers** in nonaqueous electrolyte. The pipeline: define candidates as exact precursor/product
graft pairs → screen attachment routes (proceed/hold/exclude) → build structures (capped graft,
soluble parent, finite oligomer) with explicit charge/spin/solvent/reference → run UMA and
GPU-DFT geometry + Hessian + energy jobs → compute redox potentials (vs Fc), reorganization
energies, attachment thermodynamics and finite-oligomer responses with per-state QC → publish
qualified tables/figures. Read `README.md`, `AGENTS.md`, and everything in `docs/`
(especially `MODELING.md`, `DESIGN_AXES.md`, `WORKFLOWS.md`, `PRACTICAL_GRAFTING.md`,
`PROPERTY_VALIDATION_AND_PUBLICATION.md`) before touching code.

### Scope — core pipeline only
Carry into the new tree:
- `src/redox/` (the library — this is the main job)
- `scripts/` (campaign / validation / release / plotting entry points) — curated, not wholesale
- `config/` (chemical definitions, protocols, execution policy)
- `tests/` — curated to a real contract suite (see below), not a 1:1 copy of ~250 files
- `data/reference/` and the `data/curated` inputs the pipeline actually reads; keep the
  `python -m redox seed-reference` seeding mechanism
- `docs/` — rewritten to match the new structure
- `examples/` — the portable cached replays (these are regression oracles, keep them working)
- `pyproject.toml`, `requirements*.txt`, `setup_env.sh`, `.gitignore`, `LICENSE`

Leave behind (do NOT copy): `manuscript/`, `slides/`, `validation/` (60+ dated campaign dirs),
`logs/`, `library/`, `calcs/`, `results/` bulk outputs, and `scripts/report/.browsers/` (a full
Chromium is checked in — never bring that over). These stay in the old repo; reference their
locations in docs if needed.

### The mess to fix (concrete)
- `src/redox/workflows/` is 75 files / ~11.8k LOC — a dumping ground of one-off campaign episodes.
  Collapse into a coherent, documented workflow layer.
- `src/redox/qm/` has duplicate-lineage modules accreted over time:
  `thermal` / `native_thermal` / `native_thermal_repair` / `thermal_source` / `projected_thermal`;
  `hessian_memory` / `hessian_resume` / `hessian_restart` / `saved_hessian`;
  `response_resume` / `response_diagnostic` / plus `*_checkpoints` variants. Identify the
  canonical behavior for each and unify; delete superseded copies.
- `src/redox/validation/` (45 files) and `tests/` (~250 files, many named after single campaign
  episodes like `test_counterion_T4_resource`, `test_archived_*`) need curation down to a
  maintainable, meaningful suite.
- The CLI is dispatched from a `COMMANDS` dict in `src/redox/cli.py`; `redox` and `redox-uma`
  are the console entry points. Preserve the user-facing command surface (or improve it
  deliberately and document the change).

### Method (follow in order)
1. **Map before moving.** Produce a written inventory: for each current subpackage, what it does,
   which modules are canonical vs. superseded duplicates, and what the public surface is
   (CLI commands, the `examples/` replays, the scientific-contract tests). Identify the true
   data-flow of the pipeline stages. Do not start writing the new tree until this map exists.
2. **Design the target architecture** and get it approved before bulk migration. Propose the new
   package layout (clear separation of: structure building · QM backends · thermochemistry/
   properties · screening/route gates · workflow orchestration · reporting · core utilities),
   the module consolidation plan for the `qm` duplicate clusters, and the curated test plan.
   Present it; wait for sign-off.
3. **Migrate incrementally, oracle-green at every step.** Port one subsystem at a time. After
   each, the regression oracle (below) must pass. Never let the tree sit in a state where the
   cached replays diverge.
4. **Curate, don't copy, the tests.** Keep every test that encodes a *scientific contract* or a
   *numerical value* (the Nelsen solver, QRRHO, reorg/potential math, reference conventions,
   the cached-thermochemistry and viologen/oligomer replays, UMA parity). Drop or merge tests
   that only pin a specific dated campaign's bookkeeping. The result should be a suite a new
   contributor can read.

### Regression oracle (the contract — this defines "correct")
The new code must reproduce these **exactly**:
- `python3 -I -S examples/cached_thermochemistry/reproduce.py` → capped PMDI-ethyl first
  reduction ≈ **−1.075791 V vs Fc** and the DBU separated-ion reaction ≈ **−103.068826 kJ/mol**.
- The other portable replays: `exact_viologen_replay`, `exact_viologen_corrected_replay`,
  `finite_viologen_replay`, `occupied_sampling_replay`.
- The scientific-contract subset of `tests/` (potentials, reorg, Nelsen, QRRHO, thermal,
  references, UMA parity, contracts).
Run these in the OLD tree first to capture ground-truth values, then assert the NEW tree matches.
If any number must change, STOP and flag it explicitly with the physical justification — do not
silently absorb a numerical difference.

### Hard constraints (from the user's standing instructions — non-negotiable)
- **GPU work only through the reservation gate** `~/bin/gpu_reserve` (shared NFS lock, cross-session).
  NEVER `CUDA_VISIBLE_DEVICES=X` directly. One task per GPU; a GPU is free only with zero resident
  compute processes. For this refactor you should need little/no GPU — the oracle runs on cached
  data/CPU. Do not launch DFT/Hessian jobs unless explicitly asked.
- **Cap CPU threads** per worker (`OMP_NUM_THREADS`/`MKL_NUM_THREADS`/`OPENBLAS_NUM_THREADS`) on
  this shared node.
- **Grounded facts only.** Every retained constant, reference value or literature claim must trace
  to a repo file, a fetched source, or a shown computation. Don't invent or guess numbers.
- **Real chemical names** (IUPAC/common) — never coined shorthand.
- **No scalarized figure-of-merit.** If any ranking/reporting logic survives, keep it as a
  σ-aware Pareto front over raw axes; do not reintroduce a weighted-sum score.
- Preserve the `gpu_reserve` / thread-cap / checkpoint discipline anywhere orchestration code is
  ported.

### Deliverables
- A clean, installable tree at `0_redox_refactored` (`pip install -e .` links only the `redox`
  package; runtime deps stay unpinned per the current `pyproject.toml` philosophy).
- Consolidated `src/redox/` with the duplicate `qm`/`workflows` sprawl resolved and documented.
- A curated, passing, readable test suite plus the working example replays.
- Rewritten `README.md` + `docs/` describing the new architecture and the migration (what moved
  where, what was merged, what was intentionally dropped).
- A `MIGRATION.md` mapping old modules → new locations (or → "dropped, reason").
- All regression-oracle values reproduced bit-for-bit, with the comparison shown.

### Git
- Make atomic commits per migrated subsystem with the oracle green.
- **Commit attribution: author as Ray Zhu only. Do NOT add any `Co-Authored-By` trailer.**
- Do not push unless asked.

### Do NOT
- Do not copy `manuscript/`, `slides/`, `validation/`, `logs/`, bulk `results/`/`calcs/`/`library/`,
  or the checked-in Chromium.
- Do not change scientific behavior or numerical outputs without flagging it.
- Do not create sibling export folders beyond the one destination path.
- Do not run GPU jobs or long campaigns as part of the refactor.

Start with step 1 (the map) and the target-architecture proposal. Do not begin bulk migration
until the architecture is approved.
