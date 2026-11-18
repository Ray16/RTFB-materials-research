# Where the candidate data lives (read this before computing or evaluating candidates)

The code tree (`~/0_redox_refactor`) is **code-only**. The heavy QM lives in separate data
trees. Point the pipeline at a data tree with `REDOX_DATA_ROOT` (see `src/redox/core/common.py`).

## The two data trees

| Tree | `calcs/dft` molecules | Size | What it holds |
|------|----------------------|------|----------------|
| **`~/archive/0_redox`** | **160** | **~270 GB** | **The prior candidate campaign** — the pool candidates (`viologen`, `aq_benzyloxy`, `aq_benzylamino`, `nq_benzyloxy`, `bisviologen`, `d3_*`, imides, etc.) with their UMA/DFT/Hessian/continuation data. **This is where the current candidate pool was computed.** |
| `~/0_redox` (canonical) | 19 | ~52 GB | Benchmark/reference only (`bm_*`, `ferrocene`, `methyl_viologen`, `mophquinone`, controls) — the reliability-validation set, not the candidate pool. |

**To evaluate the current candidates, use the archive:**

```bash
REDOX_DATA_ROOT=/nfs/lambda_stor_01/homes/rzhu/archive/0_redox python -m redox.reporting.tables --out <dir>
```

With that, 86 of the 90 config candidates resolve to computed DFT (95%); the only 4 missing are
the 2026-11-16 seed additions (`exp_dmpz_mb_ethyl`, `hyp_dmpz_mb_methoxyethyl`,
`hyp_aq_2_benzyloxy`, `hyp_aq_2_amino`). Note `aq_benzyloxy`/`aq_benzylamino` already exist as
computed anthraquinone candidates in the archive — the `hyp_aq_*` seeds likely duplicate them;
reuse the existing QM rather than recomputing.

## Convention: new candidates go to the SAME place (the archive)

When the pool is expanded, new-candidate QM must be written into the archive tree so it sits
with the prior candidates and can be evaluated together:

```bash
REDOX_DATA_ROOT=/nfs/lambda_stor_01/homes/rzhu/archive/0_redox \
  python -m redox.workflows.motifs --config config/custom_candidates.py \
  --only <new_ids> --out <run_dir> --execute --gpu <idx>   # via gpu_reserve
```

## Write-side wiring — SOLVED via symlinks (2026-11-16)

The pipeline uses `ROOT` as the single root nearly everywhere, and the tests + the
`occupied_sampling_replay` ORACLE isolate by **monkeypatching `ROOT`**. Swapping reads/writes to
a `DATA_ROOT` constant broke that isolation (the oracle dropped to 4/5), so the env-override is
NOT the write-side mechanism. Instead, the refactor's data dirs are **symlinked to the archive**:

```
0_redox_refactor/calcs   -> ~/archive/0_redox/calcs
0_redox_refactor/library -> ~/archive/0_redox/library
```

Now every `ROOT/'calcs'|'library'` path (and the `common` constants, which default to `ROOT`)
resolves into the archive, so **build / UMA / DFT write new candidates into the archive**
alongside the prior ones, the manifest **merges** (file-locked `locked_manifest`, preserves
existing rows — verified 433 → 445 after adding the 4 seeds), and the monkeypatch tests are
unaffected (they patch to a tmp dir, not the real symlinked root). Tests 749 pass, oracle 5/5.

`REDOX_DATA_ROOT` still works for **read-only evaluation** from a different tree (readers use the
`common` constants); it is just not used for the write path.

## Compute recipe for new candidates (writes into the archive via the symlinks)

```bash
# on a GPU node (gpu_reserve-gated); lambda2 has no GPUs
ssh -n lambda5 "cd ~/0_redox_refactor && PYTHONPATH=src OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 \
  OPENBLAS_NUM_THREADS=2 <redox-python> -m redox.workflows.motifs \
  --config config/custom_candidates.py --only <new_ids> --out results/<run> \
  --execute --gpu 4,7,8 --thermal-definition legacy"   # UMA->DFT->SP potentials
```

The 4 seed structures are already built into the archive library; the dihydrophenazine seeds'
QM is in progress on lambda5.
