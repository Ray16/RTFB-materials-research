# Working JCIM manuscript

Open [the article](manuscript.pdf) (19 pages) and
[Supporting Information](supporting_information.pdf) (32 pages).
Editable sources are [manuscript.tex](manuscript.tex) and
[supporting_information.tex](supporting_information.tex). Earlier article
versions remain preserved in Git and the recorded pre-rewrite archive; see
[restoration.json](restoration.json).

The [current immutable evidence pointer](evidence/current.json) selects the
source snapshot. The article presents a context-preserving assessment of exact
redox grafts on preformed polymers. Five quantitative interpretation-transfer
tests demonstrate consequences for reaction signs, polymer penalties, capacity
ceilings, sampled contacts and directional relaxation. These are applications
of established thermodynamic cycles; intrinsic DFT accuracy or superiority of
a new screening algorithm is not claimed.

The exact viologen molecular example now has a complete primary-model repeat:
nine new Hessians, one analytic chloride partition, 81 checked SCF frames and
independent occupied-orbital spin reconstructions. Six complete redox,
quaternization and disproportionation cycles reproduce the retained results
within 2.45e-8 V and 1.55e-6 kJ/mol across the three rotor settings. SI Table S15
reports component QC. [Calculation/review record](../validation/ionic_cycle_completion_2026-11-02/README.md).

The exact lead now also has one complete matched ionic-polymer neighboring-graft
G cycle: −3.94 kJ/mol for the adjacent all-same-backbone model, with a
−5.85 to −3.94 kJ/mol rotor-setting range. The article and SI distinguish this
finite r2SCAN-D4/SVP smooth-surface MeCN response from primary molecular
quaternization, direct motif binding, electronic coupling and an ensemble result.
Three other ionic cycles remain pending. SI Table S11 reports the four partners.

The retrospective benchmark still has 39 endpoints and MAE 0.247 V. Numerical
reproduction is distinct from experimental accuracy. Counterion/speciation
sensitivity, complete occupied-polymer cycles/sampling, material accessibility,
public deposition and final author metadata remain open. The newer alternate
polymer collection is separately scoped from the frozen primary molecular data.
See [review_status.md](review_status.md) and the
[publication acceptance criteria](../docs/PUBLICATION_READINESS.md).

Build using the redox interpreter, `PYTHONPATH=src`, and native thread pools
capped at one: `python manuscript/build.py`. The builder checks the frozen
source hashes, reconstructs claims and compiles both documents. Compiler facts
are in [build_check.json](build_check.json); PDF-hash-bound inspection is in
[visual_review.json](visual_review.json).

No Overleaf project is connected. The latest [self-contained export record](overleaf_export.json)
identifies a compiled upload archive with matching text, figures, bibliography
and numerical TeX inputs. It is a local reviewer upload, not a remote sync.

The [corrected CPU companion](../examples/exact_viologen_corrected_replay/README.md) supplies all nine new Hessians and gradient panels and passes extracted replay in a fresh NumPy-only environment with original paths blocked. It preserves the earlier example separately.
