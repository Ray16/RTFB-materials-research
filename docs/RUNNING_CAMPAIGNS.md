# Existing campaign instructions and recorded context

This preserves the detailed README present before the repository usability pass.
Use the root README for the quick start and docs/WORKFLOWS.md for entry-point roles.
Recorded timings and calibration claims below retain their original scope.

# Redox flow battery property pipeline

Computational screening of redox-active groups for post-functionalization of
preformed Merrifield-type chloromethyl-polystyrene, targeting **stability** and
**redox potential** in a nonaqueous (MeCN / PF₆⁻, Fc/Fc⁺-referenced) electrolyte.
Capped benzyl structures model individual graft sites; explicit finite oligomers
probe the local polymer environment. These model levels are reported separately.

This is the **refactored, self-contained** tree: maintained source, configuration,
curated input data, documentation, tests, and small cached examples. Bulk
calculation archives, campaign outputs, and manuscripts from before the refactor
were moved to `~/archive/0_redox` (git commit `2515f22`) and are intentionally not
carried over here. Provenance records that cite absolute paths under
`/nfs/lambda_stor_01/homes/rzhu/0_redox/{calcs,results}` refer to that archive.

## Layout

| Path | Contents |
| --- | --- |
| `src/redox/` | The installable package: build, QM, properties, screening, workflows, validation |
| `config/` | Candidate definitions, electrolyte/reference settings, execution policy |
| `data/raw/` | Curated inputs: candidate workbooks, two-wave MeCN benchmark |
| `docs/` | Modeling design, design axes, datasets, grafting notes |
| `scripts/` | Pipeline, mining, plotting and validation entry scripts |
| `tests/` | Pytest suite (self-contained via `pythonpath = ["src"]`) |
| `examples/` | Small cached thermochemistry example |
| `validation/nitro_prior/` | Source-bound Nitro force-repeatability evidence (see below) |

## Install

The package only links `redox`; the heavy environment (torch / fairchem /
gpu4pyscf / xtb) is built separately and lives in the `redox` conda env.

```bash
pip install -e .          # into the `redox` env, after the env is built
```

## Command-line entry points

`python -m redox --help` (or the installed `redox` command) lists every stage:

```
redox seed-reference    # install + verify the pinned Fc/benchmark calibration records
redox build             # build candidate structures (--dry-run validates only)
redox build-validation  # validation cores incl. the ferrocene Fc/Fc+ reference
redox uma               # gas-phase warm starts (explicit backend qualification)
redox dft               # DFT/SMD with checked warm-start provenance
redox sp                # active-protocol single points (adopts matching DFT energies)
redox native-thermal    # solution DFT Hessian thermochemistry (--select-policy)
redox thermal-postprocess  # projected qRRHO partitions from admitted Hessians (CPU)
redox status            # inspect all declared states and QC
redox review            # rebuild and audit tables in a fresh directory (no QM)
redox launch            # plan current DFT/UMA shards (--execute submits)
redox run-gpu           # preflight + reserve one GPU for a full-QM command
redox resources         # discover free GPUs through the system reservation gate
redox motifs / uma-check / energy-audit / expand / search / source-audit
```

## Running a new candidate end to end

`redox motifs` is the production orchestrator. For each registered candidate it runs
build -> UMA -> DFT/SMD (3 conformer seeds, `--no-freq`) -> active single points ->
native solution thermochemistry (ferrocene + candidate, `--select-policy --repair-minimum`)
-> Fc reference sync -> strict readiness -> reorganization cross points -> table review.
Every GPU stage goes through `gpu_reserve`; full DFT/Hessian work needs a 32 GB-class card.

```bash
python -m redox seed-reference                       # once per checkout (~7 MB, verified)
python -m redox motifs --config config/starting_candidates.py --only ID \
    --out results/runs/NAME --threads 8              # read-only plan
python -m redox motifs --config config/starting_candidates.py --only ID \
    --out results/runs/NAME --threads 8 --execute --gpu N
./scripts/pipeline/finalize_after_dft.sh             # tables, scorecard, Pareto, figures
```

`--gpu 0,1,2` (free GPUs on the launching host, see `gpu_reserve list`) optimizes a motif's
redox states concurrently, one state per GPU with its own `--threads` CPU block; every other
stage stays on the first GPU. DFT dominates the wall time and each state's three conformer
seeds still run in sequence, so the gain is bounded by the slowest state.

Measured on V100-32GB, `--threads 8` (2026-10-28):

| stage | 2-(4-benzyloxyphenyl)-1,4-benzoquinone, 3 states, 1 GPU | N-(4-methylbenzyl)phthalimide, 2 states, `--gpu 1,2,3` |
|---|---|---|
| build + UMA | 6.4 min | 0.8 min |
| DFT/SMD (3 seeds per state) | 74.7 min (neu/red1/red2 in sequence) | 42.1 min (set by red1) |
| active single points | adopted, <1 min | adopted, <1 min |
| native thermal | 51 min (incl. a minimum repair of neu) | 11.8 min |
| cross points | – | 2.8 min |

For the quinone, the per-state optimizer times (sum over seeds: neu 8.9, red1 16.9, red2
14.4 min) give an estimated ~31 min for the 3-GPU DFT stage instead of 75 (74.7 min scaled
by red1's 16.9/40.2 share; estimate, not measured). `finalize_after_dft.sh`
takes about 1.5 min on CPU plus its test step.

Finalization reports two campaign-level inputs explicitly rather than failing on a fresh
checkout: the OROP row is `NOT AUDITED` without the raw OROP dataset, and step 13 (campaign
decision publication) is `NOT PUBLISHED` without the campaign's publication benchmark and
verified candidate feed. Step 13 also needs the figure interpreter:
`python -m venv --system-site-packages .venv-figures && .venv-figures/bin/python -m pip install -r requirements-visualization.txt`.

The low-level stages (`redox uma`, `redox dft`, `redox sp`, `redox native-thermal`) remain
available for repair work; they do not apply the orchestrator's defaults (e.g. `redox dft`
uses one conformer seed unless `--nconf 3` is given).

### Calibration evidence (`data/reference/`)

Potentials are referenced to a level-matched Fc/Fc+ computed live from ferrocene
records, and the scorecard sigmas are re-derived from raw records of the grounded MeCN
benchmark molecules; serialized tables are never trusted for calibration. A fresh
checkout has none of these records, so `data/reference/validated_2026-10-19.*` pins them
exactly as read by the independently checked 2026-10-19 review
(`aryl_hydrogen_scope_tables`, 32,492 arithmetic checks): Fc_abs = 4.379195901594358 V,
family sigma_E imide 0.142 / viologen 0.134 / phenothiazine 0.322 V (pooled 0.246 V),
sigma_disp 0.575 eV, sigma_lambda 0.209 eV. `seed-reference` installs them (never
overwriting a differing file) and verifies the live recomputation; finalization then
requires `results/uncertainty_evidence.json` to equal the pinned evidence
(`seed-reference --compare-results`). Benchmark records added to the archive after that
review change these sigmas and are deliberately *not* included: adopting them is a
calibration update that needs its own review (rebuild with
`scripts/reference/build_reference_bundle.py`).

The scorecard refuses to run on config-prior sigmas for E, disproportionation or lambda
(`--allow-prior-sigmas` / `REDOX_ALLOW_PRIOR_SIGMAS=1` accepts them explicitly, and the
choice is recorded in `uncertainty_evidence.json`). Charge-class E sigmas remain priors
in the validated state (no OROP row passes the matched-reference gate) and are labelled.
The native thermal certificates cite the original checkout path
`/nfs/lambda_stor_01/homes/rzhu/0_redox`; a checkout elsewhere needs
`redox.qm.relocated_native.relocated_native_reader`.

## UMA backend policy

The backend is chosen from the **actual reserved GPU's** compute capability, not
from a config preference:

- **Stock (legacy) backend by default** — runs on a V100 (cc 7.0) and any
  non-Nitro-capable card. It is deterministic (force repeatability ≈ 10⁻⁶ eV/Å).
- **Nitro** is selected only on a Nitro-capable card (CUDA compute capability
  8.x–12.x, i.e. an H100). It requires Hopper/Ampere-class hardware and **cannot
  run on a V100**.
- Nitro additionally stays blocked by a production-qualification gate
  (`config/uma_nitro_qualification.json`) until its force repeatability passes the
  1×10⁻⁴ eV/Å tolerance. The retained H100 evidence under
  `validation/nitro_prior/` currently shows it does **not** (fused-kernel TF32 +
  atomic scatter nondeterminism), so production falls back to the stock backend.

GPU work must go through the system reservation gate (`gpu_reserve`); the pipeline
verifies lease ownership before creating a CUDA context.

## Tests

```bash
pytest            # self-contained; uses this tree's src via pyproject pythonpath
```

## Documentation

- [Modeling the system](MODELING.md)
- [Design axes and reliability](DESIGN_AXES.md)
- [Datasets](DATASETS.md)
- [Practical grafting](PRACTICAL_GRAFTING.md)
- Agents should read [AGENTS.md](../AGENTS.md) before changing or running the pipeline.
