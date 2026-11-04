Submission data companion — 3 November 2026

Read benchmark_endpoints.csv together with benchmark_context.json. The 39
endpoints cover 30 parent molecules. The model state charges and multiplicities,
solvent, temperature, reference, electronic/native protocols and thermal
convention are explicit. Experimental potentials retain their source DOI;
experimental conditions have not been homogenized by this export. Consult the
source records and SI for per-study handling and field-specific validation.
Parent measurements are not experimental validation of the exact grafts.

The 912 rows of design_SMILES.csv identify exact capped local graft designs.
Enumeration does not imply calculation completion, qualified thermochemistry,
synthetic accessibility or polymer performance. reported_evaluation_count also
includes historical evaluations; it is not a count of fully qualified properties.
Source-role and alias columns are JSON lists. Route holds/proceeds and material
decisions remain separate. Product and parent graph representations do not
establish charge or spin localization. The registry retains the full precursor,
attachment and property evidence; this CSV is an identity inventory, not a
replacement for that registry.

manifest.json binds the exported artifacts and the original source files by
SHA-256. Absolute source paths document provenance; they are not required to read
these tables. Sources are in the project repository on refactor/clean-workflow:
results/candidates/pool/2026-11-03_qualified_phenothiazine_context/registry.json
manuscript/evidence/2026-11-03_seven_chloride_outcomes/benchmark.json
manuscript/evidence/2026-11-03_seven_chloride_outcomes/screening_energy_inputs.json
manuscript/evidence/2026-11-03_seven_chloride_outcomes/primary_native_cohort_inputs.json
scripts/release/export_submission_tables.py

The export checks exact design identity and benchmark residual arithmetic. It
neither repeats quantum calculations nor promotes material eligibility. The
older 2026-11-03 directory is a preserved prototype; use this context-complete set.
