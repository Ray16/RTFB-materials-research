# Exact molecular CPU replay

The [versioned example](../../examples/exact_viologen_replay/README.md) provides
ten small cached states for the same exact viologen redox and quaternization
case. The raw numerical payload is about 1.2 MB; its deterministic archive is
651,404 bytes. It reconstructs nine projected Hessian mode sets/nuclear
partitions, six saved finite-difference panels and one monatomic partition.

The extracted archive passes in a fresh hash-locked NumPy 2.2.6 environment
without the project, PySCF, ASE or Torch installed. An audit hook blocks the
original checkout, archive, molecular calculation roots and polymer
calculation roots. It reproduces both reduction potentials, all three
precursor-specific quaternizations and monomer disproportionation at
25/50/100 cm-1. Maximum frequency and nuclear errors are below 1e-10 cm-1 and
5.4e-14 eV, respectively.

The main result remains −11.49 kJ/mol for the finite local-site separated-ion
quaternization in MeCN at 298.15 K. This replay verifies numerical records and
bookkeeping, not experimental accuracy, new SCFs, a polymer material or
complete historical executable/scratch closure. Historical Fc/local-site FD
summaries remain metadata-only. Public archival DOI and full deposition remain
open; the small example is versioned with the code.

`validation.json`, `replay.json`, `package_manifest.json` and `compression.json`
bind the extracted test and package. `tests.xml` records 24 example/checker
tests; `publication_tests.xml` records 51 checks including those 24. The current
21/28-page article/SI include the example and its scope. Revised pages are
visually checked; all nineteen deck frames match their previously checked text
and renders. The self-contained Overleaf upload compiles with matching text;
no Overleaf project is connected.

The bulk of the wider project remains on shared storage. This small example
does not replace the larger benchmark, PMDI/TEMPO, polymer-state or production
runtime archives. Their outstanding release and scientific requirements remain
listed in `docs/PUBLICATION_READINESS.md` and `manuscript/review_status.md`.
