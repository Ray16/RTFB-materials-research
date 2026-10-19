# Repo map

A guide to maintained code and saved evidence. Start with [AGENTS.md](../AGENTS.md)
for operating rules, [Development](DEVELOPMENT.md) for supported commands, and the
[documentation index](README.md) to distinguish current guidance from historical reports.
`DESIGN_AXES.md` defines properties and QC; `DATASETS.md` locates external data.

## `src/redox/` — the `redox` package (import as `redox.<subpackage>.<module>`)

`pip install -e .` once (editable install via `pyproject.toml`), **or** prefix commands with
`PYTHONPATH=src`. Every CLI module runs as `python -m redox.<subpackage>.<module>`.

**`core/` — shared helpers and data contracts**
- `common.py` — repo paths, physical constants, `load_config`, `read_manifest`,
  `read_result` (returns ACTIVE-protocol energies only; `sp_status="missing"` otherwise),
  `state_names` (real redox states only — never `reorg/`, `sp/`, `_*`, `*.bak*`),
  `free_energy` (G = E_smd + g_thermal; None if either is missing — never a silent 0).
- `protocol.py` — **the energy protocol** (`ACTIVE_SP`: ωB97M-V/def2-TZVPD for EVERY charge
  state and Fc, SMD MeCN, RI-J, conv 1e-9) and the protocol hash that addresses every energy
  record `calcs/dft/<id>/<state>/sp/<hash>.json` (hash = geometry + charge + spin + level +
  solvent + SCF settings).
- `io.py` — atomic JSON/text replacement for evidence and mutable status records.

**`build/` — structure generation** (SMILES → 3D conformers into `library/`, stereo-checked)
- `structures.py` — shared builder (decorate, conformer ensemble, stereo check, XYZ).
- `candidates.py`, `standalone.py`, `validation.py` — the three library sets (grafted
  candidates, "before functionalization" references, known-E validation cores).
- `manifest.py` — shared lock, column schema and atomic CSV replacement used by all builders.
- `geometry.py` — declared SMILES to stored XYZ mapping, stereochemistry and redox-core
  indices. Established imports from `core.common` remain compatible.

**`qm/` — electronic structure**
- `uma.py` — UMA charge/spin-aware gas pre-opt.
- `dft.py` — DFT+SMD geometry optimization (r2SCAN-D4/def2-SVP(D) in SMD) + single points.
  `run_batch` refuses to silently accept a `result.json` from a different optimization
  protocol (`[STALE]`, needs `--force`). `torsion_scan.py` — rotated seeds for floppy species.
- `sp.py` — compute/adopt the active-protocol energy record for every state
  (`--audit`, `--all --shard n:i`, `--orop` for the benchmark tree).
- `native_thermal.py`, `thermal.py`, `vibrations.py` — addressed native solution
  thermochemistry, thermal identities and spectrum/qRRHO integration.
- `projected_thermal.py`, `qrrho.py`, `xtb_kernel.py` — corrected internal-mode
  partition caches, rigid-mode projection and an isolated installed xTB CPU kernel.
  `python -m redox thermal-postprocess` reuses existing qualified Hessians.
- `hessian_memory.py` — independently qualified execution budgets and derivative controls.
- `nonequilibrium_pcm.py`, `nonequilibrium_smd.py` — frozen solvent response implementations.

**`properties/` — per-molecule screening axes**
- `potentials.py` — E° vs a LIVE level-matched Fc/Fc⁺ (same G definition); INCOMPLETE couples
  flagged → `results/redox_potentials.csv`.
- `integrity.py` — bound + intact charged-state screen (`intact_bound`); a necessary
  condition, NOT electrochemical reversibility.
- `stability.py` — disproportionation ΔG (= F·ΔE₁₂).
- `reorg.py` — Nelsen 4-point λ (λ_O, λ_R, sum) with protocol-addressed cross-point caches.
  `nelsen.py` — on-the-fly 4-point for validation scripts.
- `lambda_outer.py` — outer-sphere λ_o (Born, molecular-cavity nonequilibrium PCM; sign-checked,
  point-charge sphere test). `solvated_reorg.py` — SASA radius + Born helper.
- `capacity_and_proxies.py` (capacity, SA score, solubility proxy). (Couple RMSD for λ QC lives in
  `reorg.py`; the separate RMSD-table module `descriptors.py` is archived.)

**`screening/` — ranking**
- `scorecard.py` — one row per (molecule, electrode pool): contiguous redox path from the
  resting state, explicit λ conventions (λ_het, λ_SE), INCOMPLETE propagation.
- `pareto.py` — σ-aware Pareto front on complete rows only (incomplete rows cannot dominate).
- `discovery.py` — systematic candidate identification: family patterns with a 2e⁻ rule backed
  by measured two-wave data, Merrifield grafting as RDKit reactions (O-benzyl ether, N-benzyl
  amine/imide, Menshutkin only when it completes a viologen, S-benzyl, benzyl ester) onto the
  4-methylbenzyl tether proxy, stereo/element rules, capacity + SA computed like the scorecard.
  `rediscover()` must regenerate every registered candidate (tested).
- `expansion.py` — matched-comparison hypotheses, exact-graph deduplication,
  explicit evidence scope and append-only registration. `python -m redox expand --out DIR`.
- `family_search.py` — systematic core/substituent/linker/ring enumeration with
  symmetry deduplication and recipe provenance. `python -m redox search --out DIR`
  prepares a pool without registration or QM. Search axes are declared in
  `config/family_search.py`; complementary sources in `config/discovery_sources.json`.
- `source_audit.py` — downloaded-source quality inventory, preserving ordinal
  optimizer grades, raw identity flags, source hashes and explicit label-admission scope.

**`validation/` — accuracy checks used by the pipeline**
- `stability.py` — ΔG_disp vs experimental wave spacing → σ_disp.
- `readiness.py` — all manifest states, including missing geometry, protocol/spin,
  molecular identity and selected thermal admission.
- `spins.py` — source-addressed spin comparisons and promotion evidence.
- `provenance.py` — referenced evidence file discovery.
- `orop.py` — external benchmark aggregation with current active energies and
  matched molecular/Fc thermal definitions. Historical/mismatched rows cannot
  calibrate the current scorecard.

**`workflows/` — supported orchestration**
- `tables.py` — isolated CPU table rebuild with child-process stages, logs and final
  independent audit. `python -m redox.workflows.tables --out results/reviews/NAME`.
- `resources.py` — shared CPU allocation, node suitability and system GPU discovery.
- `gpu.py` — single-command resource preflight, process ownership and reservation.
- `launch.py` — plans/runs resumable DFT/UMA shards on all suitable free local GPUs;
  explicit global shard ranges support multi-node allocation.
- `revalidate_sp.py` — frozen-input independent recomputation of historical SP
  adoptions; archived originals and checked comparisons precede canonical promotion.

## `scripts/` — operational entry points and campaign drivers

[Script guide](../scripts/README.md) describes supported paths. Reusable orchestration
is moving into `src/redox/workflows/`; historical campaign logic still lives here.
- `check_env.py`, `free_gpus.py` — env check and cluster GPU scan (paths used by `CLAUDE.md`).
- `pipeline/` — the production chain: `run_uma.sh`, `run_dft.sh` (compatibility wrappers
  through the resource-aware package launch workflow),
  `finalize_after_dft.sh` (fail-fast publication: tests → readiness → Fc → E° → saved λ →
  integrity → stability → λ_o → capacity → scorecard → Pareto → figures), `set_fc_reference.py`,
  `compute_lambda_outer.py` → `consolidate_lambda_outer.py`, `fill_candidates_xlsx.py`.
  `run_motifs.py` plans candidate stages by default; `review_tables.py` is the compatibility
  wrapper for the package workflow. Collection, queue and slide drivers remain here.
- `analysis/` — robustness checks on the candidates: `thermal_sensitivity.py` (RRHO cutoff
  spread), `reorg_conformer_matched.py` (tether-rotation-free λ).
- `validation/` — benchmark campaigns, one folder each:
  - `orop/` — `run_orop_benchmark.py` (external E° benchmark, same protocol as production),
    `orop_backfill_thermal.py`.
  - `solvation/` — `validate_solvation.py` (FreeSolv/MNSol ΔG_solv).
  - `reorg_d3tales/` — D3TaLES λ reproduction: set builder, workers at THEIR level
    (`validate_reorg_worker_d3tales.py`) and OUR production level (`…_prod.py`), aggregation +
    unbound-anion QC (`aggregate_d3tales_reorg_prod.py`, `finalize_reorg_qc.py`,
    `reorg_screen_unbound.py`, `reorg_recompute_guarded.py`).
  - `reorg_literature/` — literature λ anchors (oligoacenes, perfluoroacenes):
    `validate_reorg_worker_{generic,inner_acene}.py`, `preflight_reorg_set.py`.
  - `reference/` — `scf_level_check.py`: (1) SCF ground-state check — every state re-solved
    from several initial guesses (`--set reference` = anchors/Fc/OROP sample; `--set
    candidates` = every energy feeding a ranked candidate), (2) level crossing at OROP's
    B3LYP-D3/6-31G* (SMD and C-PCM) → `results/validation/{scf_ground_state_check,
    level_crossing_Efc,scf_candidates_check}.csv`; `scf_check_phase.py` re-runs one phase in a
    fresh process.
  - `audit/` — `recompute_axes.py`: independent from-raw recompute of every published axis
    (imports only `redox.core.protocol` + config) → `results/validation/audit_*.csv`.
- `mining/` — D3TaLES family coverage + candidate mining; `identify_candidates.py` runs
  `redox.screening.discovery` over D3TaLES → `results/discovery/` (pool, rediscovery check,
  per-family pre-filter front).
- `fleet/` — compatibility launch paths through `redox.workflows.gpu`/`launch`.
  Old claim-based workers and their launchers now stop with a migration message;
  original sources are preserved as text in `scripts/legacy/fleet/`.
  Historical roles were:
  (`dft_launch.sh` still supplies one global shard numbering across nodes;
  `sp_fleet_worker.sh` = SP records + λ cross points (`CROSS_IDS`) + OROP SPs;
  `scf_level_launch.sh`; d3level/reorg/OROP workers + launchers; `gpu_probe.sh` +
  `cluster.env` back `free_gpus.py`).
- `polaris/` — ALCF Polaris env build + PBS jobs.
- `plotting/` — **all figures**, one folder per `results/figures/` folder
  (`candidates/`, `validation/`, `reorg/`, `pipeline/`); every `plot_*.py` uses
  `plotting/plot_style.py` `apply_style()`.

## Adding candidates (the growing registry)
1. `python scripts/mining/identify_candidates.py [--pubchem]` — systematic identification
   (rules regenerate every registered candidate); writes `results/discovery/`. Already
   registered molecules are marked `parent_in_library` and never resurface as new.
2. Name the leads (real names from the PubChem parent) + chemistry flags in a
   `data/raw/candidates/discovered_names_<date>.csv`, then
   `python scripts/mining/register_candidates.py --from <front csv> --names <csv>` —
   APPEND-ONLY into `config/discovered_candidates.py` (provenance recorded, idempotent).
3. Prepare the complete selected-candidate pipeline with
   `python scripts/pipeline/run_motifs.py --config config/discovered_candidates.py --only ID --out results/campaigns/NAME`.
   Review the printed stages, then add `--execute --gpu INDEX` on a preflighted
   suitable node under the fleet rules. The runner handles building, UMA, DFT,
   active SPs, native thermochemistry, state QC and gas cross points. Use a new
   destination for a new run, or `--resume` for the same preserved inputs/methods.
4. Collect fields through the verified candidate feed and run
   `python -m redox.workflows.tables --out results/reviews/NAME` for an isolated
   table audit. Inspect its arithmetic and input-freshness results before publication.
   Registered discoveries enter the scorecard through `CANDIDATE_CONFIGS`; parent
   measurements never become exact-graft measurements automatically.

For a stratified family batch, start with `python -m redox search --out DIR` and
the [exploration guide](CANDIDATE_EXPLORATION.md). The pool and comparison selector
are preparation steps; registration and expensive calculations remain explicit.

## `tests/`
Tests cover thermodynamic identities, state ordering, protocol/cache invalidation,
contiguous capacity, Pareto completeness/σ, PCM sign + Born limit, native derivative
admission, molecular identity, spin decisions, publication freshness, manifest concurrency
and workflow isolation. `test_refactoring.py` exercises shared ownership and failure cases.
`PYTHONPATH=src python -m pytest -q tests` (also step T of the finalize chain).

## `config/` — parameters (no logic)
`project.json` + `electrolyte.py` (solvent, window, Fc record), `validation.py` (hand
anchors, each event flagged `grounded`), `benchmark_mecn.py` (sourced MeCN E1/E2/wave-spacing
benchmark from `data/raw/validation/two_wave_mecn/`, tiered A/B/excluded per wave),
`standalone.py`, `starting_candidates.py`, `merrifield_multielectron.py`, `redox_groups.py`,
`scorecard_config.py`.

## Data and outputs
- `data/raw/candidates/` — the input candidate sheets (`Candidates.xlsx`,
  `merrifield_multielectron_smiles.xlsx`). `data/raw/validation/` — external datasets
  (git-ignored) except `two_wave_mecn/`, our hand-curated literature table (tracked; every row
  has DOI + table/page + original reference electrode).
- `library/` — generated structures + `manifest.csv` (resting state = `n_e == 0`).
- `calcs/dft/<id>/<state>/` — `opt.xyz` + `result.json` (geometry, thermal) + `sp/<hash>.json`
  (energies) + `reorg/` cross points. `calcs/spincheck/<id>/<state>_m<k>/` — DFT checks of the
  alternative multiplicity. `calcs/orop/` — benchmark states.
- `results/*.csv` — axis tables (`redox_potentials`, `state_integrity`, `reorganization`,
  `stability_disproportionation`, `lambda_outer`, `scorecard`, `pareto_shortlist`,
  `orop_benchmark`, …), `results/Candidates_computed.xlsx`, `results/figures/{candidates,
  validation,reorg,pipeline}/` (PNG).
- `archive/` — historical location, empty in the reviewed working tree. Do not rely on
  old archive links without checking whether the referenced file exists.

## Git tracking

Some generated `calcs/`, `library/` and `results/` files are tracked. Many more are
untracked or ignored. `paper/`, downloaded external validation directories, logs,
model weights and tokens have explicit ignore rules. Inspect `.gitignore` and
`git status` before staging or removing generated evidence. This refactor preserves
scientific paths and existing local/staged work.
