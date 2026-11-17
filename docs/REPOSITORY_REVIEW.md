# Repository release review

Initial assessment: GitHub main at e25cf5ff. At the user's instruction, all
work was consolidated into the canonical 0_redox directory. It is now on main;
pre-existing unrelated edits remain local. The older main
checkout and its review patch are preserved in Git metadata. The reconciliation preserves both committed histories and excludes unrelated
local edits from the publication checkpoint.

## Assessment and prioritized plan

The package constructs exact local graft models, prepares UMA/DFT calculations,
and combines qualified states into redox, reaction and oligomer comparisons.
Scientific responsibilities already have separate modules. Protocol-addressed
energies, thermal definitions, checkpoint ownership and admission thresholds are
essential behavior to preserve. The cached PMDI example provides a portable
starting point independent of the full scientific environment.

1. Make the README and installation instructions accurate for a clean clone.
   Explain editable installation, optional dependencies, inputs and outputs;
   foreground the runnable standard-library example and current manuscript.
2. Make entry points discoverable by role without importing calculation modules.
   Preserve existing command dispatch and legacy aliases.
3. Consolidate Python configuration loading and provide actionable errors for
   missing configurations without changing how valid configurations execute.
4. Add a minimal CPU CI suite, contributor instructions and regression coverage
   for the changed interfaces. Verify installation and documented commands.
5. Record actual test results and remaining portability/release limits.

No numerical kernels, scientific constants, random seeds, force backends,
reference conventions or physical QC thresholds are changed. Existing data and
historical results are retained. The existing MIT license is retained. This pass
adds no publication metadata or new license. Validation records describe the pre-publication worktree.

## Implemented in the canonical checkout

- Replaced the entry README with an executable cached quick start, workflow map,
  installation distinctions, output descriptions and scientific/release limits.
  Preserved its prior detailed campaign content in `docs/RUNNING_CAMPAIGNS.md`.
- Added `docs/INSTALLATION.md` and `CONTRIBUTING.md`. Marked the historical dataset
  inventory as dated context rather than a clean-clone guarantee.
- Consolidated trusted Python configuration loading between shared utilities and
  the candidate builder, preserving valid module execution and original errors.
  Missing files/required fields now have explicit diagnostics.
- Added a minimal GitHub CPU workflow for the existing standard-library cached
  example and lightweight tests. It has read-only repository permissions and
  launches no scientific jobs. Hosted CI subsequently passed; the current
  workflow includes 34 CPU interface, cache, relocation and release-gate tests.
- Updated dispatch test doubles to supply source provenance and isolate synthetic
  hosts from live campaign priorities. Production dispatch behavior is unchanged.
- Added conventional ignore patterns while preserving existing ignore rules.
  No raw records, checkpoints, models or license terms were removed or changed.

## Remaining high-value work

The seeded-reference failure was traced to the qRRHO library-location refactor:
the original certificate bound the previous literal path, while the maintained
source discovers the same hash-qualified library. A narrowly defined source
compatibility check now reconstructs the complete original source hash after
reversing only that locator replacement. Additional source changes are rejected;
the original certificate, actual library, worker, protocol and evidence checks
remain required. Thirty focused tests pass, including the live pinned-reference
integration and rejection of altered constants, cutoffs, imports and locator code.
The original failing run is preserved separately from the recovery evidence in
`calcs/validation/repository_reference_recovery.json`. This is source/runtime
compatibility, not new experimental validation or a new quantum calculation.

The canonical checkout retains pre-existing scientific edits. The publication
checkpoint joins the existing committed scientific history with main while
leaving those unrelated edits local. The running polymer calculations
must finish before their final relocation; temporary sibling links remain until
readers and namespaces are migrated. The full declared raw archive is now
verified, and the manuscript's method/reference clarification points are resolved.
Public deposition, data reuse terms and author declarations remain pending.

The machine-readable validation result is `validation/repository_review.json`.

The [current handoff](../validation/publication_release_2026-11-03/current_handoff.json)
records the completed local deliverables and the remaining input/state dependencies.
Earlier validation reports retain their original observation times and scope.
