# Validation and candidate discovery — 2026-10-10

The regenerated tables and workbook pass **8,306 checks with zero discrepancies**.
The software supports new motifs through validation, grafting, computation, and audited
export. Scientific accuracy remains limited by incomplete calculations and the model
limitations below. The authoritative snapshot for this review is
[`results/review_2026-10-10_final`](../results/review_2026-10-10_final/).
This is the preserved first audit snapshot. The subsequent parallel repair campaign and
39-proposal expansion are tracked in [`PARALLEL_PROGRESS_2026-10-10.md`](PARALLEL_PROGRESS_2026-10-10.md).
Its current status and eventual post-repair export have separate manifests.

## ✅ Verification

| Check | Outcome | Evidence |
|---|---|---|
| Independent arithmetic, membership and export audit | 8,306 checks; zero discrepancies | `validation/extended/audit_recompute.csv` in the review snapshot |
| Regression and chemistry contracts | 48 tests pass | `tests/test_contracts.py`, `tests/test_readiness.py` |
| New motif grafting | Real RDKit embedding; valid state ladder; rerun preserves all files; changed chemistry under the same ID is refused | `test_builder_real_embedding_resume_and_collision` |
| Pipeline integration | Build → UMA → DFT → active SP → coverage → cross points → audited tables passes for cached viologen | `results/pipeline_smoke_2026-10-10_final/run_manifest.json` |
| Environment | ENV OK, including xTB | `logs/env_review_2026-10-10.log` |
| Declared-state coverage | 197 states: 171 ready, 15 REVIEW, 11 INCOMPLETE | `validation/readiness.csv` in the review snapshot |

The integration run exercised reuse of existing quantum calculations. A complete new
motif was embedded in the chemistry tests; no new motif underwent a fresh full DFT run
during this review. Existing production outputs and pre-existing repository edits were
preserved. New exports are separate snapshots.

The initial audit found 11 discrepancies: a missing imaginary-mode flag, five missing
integrity rows, and five capacity/proxy entries that no longer matched the calculation
tree. Regeneration resolved them. Further code review found two semantic errors:
partially valid reorganization energies could qualify a two-electron molecule for the
primary Pareto front, and the multi-electron summary could combine oxidation and reduction
of an ambipolar molecule into a nominal two-electron anolyte. Both are fixed and tested.

The workbook now contains the full scorecard, Pareto, potential, and reorganization
tables alongside the original six-candidate comparison. Every exported table cell is
traced to its CSV and spreadsheet error cells are checked.

## 🔬 Scientific limits

Numerical agreement with the inputs establishes reproducibility, not agreement with
experiment. The revised experimental benchmark has 34 scored tier-A potential values:
MAE 0.223 V and bias +0.177 V. Eight tier-A reference-free spacings have MAE and positive
bias 0.258 V, approximately 24.9 kJ/mol. Absolute potentials and disproportionation
stability therefore retain systematic error; no fitted correction was applied.

Four quinone entries were checked directly against Wang et al., SI Table S1, page S4.
Their E1/E2 transcriptions match. The fluoro-benzoquinone row is internally inconsistent:
the printed potentials imply 0.94 V whereas the printed spacing is 0.86 V. Its second
wave and derived spacing are now tier B; its unambiguous first wave remains tier A.[^1]
Other literature entries retain their prior provenance and have not all been independently
rechecked against their primary PDFs in this review.

The 11 incomplete declared states are in the earlier eight-molecule discovery batch.
Their missing active energies/geometries remain blank. Imaginary xTB thermal modes in
15 otherwise complete states remain explicit review flags. Gas-phase dianion reorganization
and conformer jumps remain flagged; no gas-phase value is promoted to a validated solution
reorganization energy. PCM outer-sphere values remain annotations pending validation
against an independent nonequilibrium implementation. DFT spin validation for small UMA
gaps, decomposition chemistry, and polymer-phase behavior remain open.

Only **viologen** is on the primary complete Pareto front in both salt scenarios in this
snapshot. This reflects the available QC-clean calculations, not proof that the other
motifs are inferior. Partial kinetics cannot dominate complete candidates.

## ⚙️ Adding and running motifs

Add entries to `config/custom_candidates.py` with a unique ID, name, family, one terminal
`[*:1]` attachment in `frag`, and `(state, charge, multiplicity, n_e)` tuples. Exactly one
resting state has `n_e=0`; `n_e = charge - resting_charge`. Every state must have compatible
electron-count/spin parity, and the ladder must be contiguous. Unassigned stereochemistry,
disconnected structures, duplicate IDs/states, and changed chemistry under an existing ID
are refused before files are written.

```mermaid
flowchart LR
    accTitle: Motif Calculation and Validation Workflow
    accDescr: Motifs are validated and grafted before quantum calculations. Complete calculations are exported and independently audited; missing or flagged data remain explicit.
    motif_input["Motif and charge states"] --> validate_input["Validate and graft"]
    validate_input --> optimize_geometry["UMA then DFT in SMD"]
    optimize_geometry --> calculate_energy["Active energies and cross points"]
    calculate_energy --> assess_quality["Coverage and physical QC"]
    assess_quality --> export_tables["Scorecard and workbook"]
    export_tables --> audit_values["Independent audit"]
```

Activate the reproducible environment:

```bash
source /nfs/lambda_stor_01/homes/rzhu/miniforge3/etc/profile.d/conda.sh
conda activate redox
export PYTHONPATH=src
```

Validate and inspect the plan:

```bash
python -m redox.build.candidates --config config/custom_candidates.py --dry-run
python scripts/pipeline/run_motifs.py --config config/custom_candidates.py --out results/my_run
```

Execute the selected registered motifs on one reserved GPU:

```bash
python scripts/pipeline/run_motifs.py --config config/custom_candidates.py \
  --out results/my_run --execute --gpu 0
```

The runner uses the shared GPU reservation gate, thread caps and CPU affinity, per-motif
writer locks, per-stage logs and a run manifest. The default DFT search uses three seeds
per state. It stops on missing/stale geometries, failed SCF, absent thermal corrections,
invalid cross-point caches, or audit discrepancies. It never requests `--force`. Failed
partial calculations must be inspected before replacing them; choose a new run output
directory when resuming. Existing valid molecular calculations are reused.

Outer-sphere PCM is not automatically computed for new motifs by this entry point: that
axis remains missing until its separate calculation is run, and is not a primary ranking
objective. The exported snapshot can contain incomplete unrelated molecules; the selected
motifs must pass their own coverage gate. The legacy global finalizer now fails on any
incomplete declared state.

## 🧪 Candidate discovery

The refreshed D3TaLES search excludes both known parents and already-grafted structures.
It finds 104 new graftable quinone molecules; 50 have a parent D3TaLES electron
reorganization value within the preliminary 0–1.5 eV filter. Ten molecules are on the
preliminary capacity/SA/parent-reorganization front. These parent values are not production
properties of the grafted molecules. Results:

- `results/discovery/review_2026-10-10/d3tales_pool.csv`
- `results/discovery/review_2026-10-10/d3tales_prefilter_front.csv`
- `results/discovery/next_candidates_2026-10-10.csv` — ten quinones plus nine imide proposals

D3TaLES supplies few bis-imides with an available graft handle. The new imide generator
therefore produces four new PMDI and five NDI designs with methyl, ethyl, isopropyl, or
short ether-containing N-substituents. It removes the existing methyl-capped PMDI, validates
every charge ladder and regenerates every product from an N-H precursor using the same
grafting rules. PMDI designs are motivated by nonaqueous multi-electron imide studies;[^2]
short ether-bearing NDI derivatives are motivated by nonaqueous NDI functionalization
studies.[^3] These are proposed derivatives, not molecules with verified measurements.

Prioritize chemistry review of non-halogenated quinones `80NOPT`, `80LEVD`, `80GKVH`, and
`80XWVS`, and the imide proposals `design_pmdi_methoxyethyl`, `design_ndi_methyl`, and
`design_ndi_methoxyethyl`. This is a choice of informative chemistry to test, not a
prediction of the final Pareto front. Halo-quinones need checks for competing nucleophilic
substitution; aza-quinones need attachment selectivity checks; ester tethers and residual
protic groups need explicit stability assessment. None of the 19 new proposals has been
registered or launched for production quantum calculations.

Regenerate into fresh directories:

```bash
python scripts/mining/identify_candidates.py --out results/discovery/my_search \
  --scorecard results/review_2026-10-10_final/scorecard.csv
python scripts/mining/propose_imides.py --out results/discovery/my_imides
```

The next computational campaign should first repair the 11 incomplete states, address
solution reorganization and conformer matching for flagged candidate steps, and check
small spin gaps. Promote selected new proposals to `custom_candidates.py` once their
grafting chemistry is reviewed, then run the automated pipeline and compare them using
the same QC and uncertainty rules.

[^1]: Wang et al. (2020), Effect of Cycling Ion and Solvent on the Redox Chemistry of Substituted Quinones. [Primary paper](https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632); [SI Table S1](https://acs.figshare.com/articles/journal_contribution/12493136). Local downloaded SI: `results/validation/review_2026-10-10/jp0c03632_si_001.pdf`.
[^2]: Daub et al. (2021), [Imide-Based Multielectron Anolytes as High-Performance Materials in Nonaqueous Redox Flow Batteries](https://pubs.acs.org/doi/10.1021/acsaem.1c01490).
[^3]: Ahn et al. (2025), [Ethylene Glycol-Functionalization of Naphthalene Diimides in Nonaqueous Redox Flow Batteries](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cssc.202501011).
