# Archive log

Records files relocated out of the canonical checkout to
`/nfs/lambda_stor_01/homes/rzhu/archive/0_redox/`. Content is **preserved, not
deleted**; tracked files also remain recoverable from Git history. Immutable
scientific evidence (`validation/`, `manuscript/evidence/`, release inventories)
stays at its original path, per `AGENTS.md`.

## 2026-11-15 — superseded Overleaf upload snapshots

Moved to `archive/0_redox/manuscript/`. Each is a draft/intermediate Overleaf
upload that a byte-distinct `_final` sibling in the repo supersedes; none is
referenced by any doc, script, manifest or the live `manuscript/overleaf/`.

| Archived (draft) | Superseded by (kept in repo) |
| --- | --- |
| `overleaf_upload_2026-10-31_exact_CPU_example.zip` | `..._exact_CPU_example_final.zip` |
| `overleaf_upload_2026-11-02_three_ionic_cycles.zip` | `..._three_ionic_cycles_final.zip` |
| `overleaf_upload_2026-11-02_direct_phenothiazine_pool.zip` | `..._direct_phenothiazine_pool_final.zip` |
| `overleaf_upload_2026-11-03_editorial_pass.zip` | `..._editorial_final.zip` |

## 2026-11-16 — stale `validation/` dirs (moved to archive/0_redox/validation)

Dated dirs superseded by a later one and not referenced outside themselves
(only this log named them). Moved on approval:

- `validation/ionic_electronic_cycles_2026-10-31` → `ionic_cycle_completion_2026-11-02`
- `validation/manuscript_polymer_replay_2026-10-30` → `publication_native_replay_2026-10-28` / later replays
- `validation/polymer_energy_continuation_2026-10-30` → `polymer_energy_classes_2026-10-30`
- `validation/polymer_native_continuation_2026-10-30` → `polymer_continuations_checked_2026-10-30`
- `validation/route_gate_identity_2026-10-28` → `candidate_grafting_audit_2026-10-28`
- `validation/viologen_native_2026-10-30` → `viologen_native_2026-11-02`

**Held back (do NOT move):** `charged_native_pilot_2026-10-31`,
`polymer_collection_2026-10-31`, `polymer_recovery_2026-10-31`,
`publication_viologen_2026-10-31` — each is referenced by the frozen publication
raw-inventory manifests under `results/redox_release_candidates/` and by sibling
`native_census.json` evidence. Relocating them would make the release inventory
inconsistent, which `AGENTS.md` prohibits.
