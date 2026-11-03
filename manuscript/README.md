# Working JCIM manuscript

Open [the article](manuscript.pdf) (21 pages) and
[Supporting Information](supporting_information.pdf) (35 pages).
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
within 2.45e-8 V and 1.55e-6 kJ/mol across the three rotor settings. The SI
reports component QC. [Calculation/review record](../validation/ionic_cycle_completion_2026-11-02/README.md).

The exact lead now has all four declared matched ionic-polymer neighboring-graft
G cycles: −3.94 kJ/mol (all-same adjacent), +22.58 (all-same one spacer),
+8.28 (alternating adjacent), and −5.30 (alternating one spacer) at the primary
cutoff. Twelve checked native components supply these finite r2SCAN-D4/SVP
smooth-surface MeCN responses. Both signs reverse with backbone configuration
and remain unchanged across the tested rotor settings. Occupied ensembles remain pending.
Twelve frozen motif triplets additionally supply 72 checked gas/MeCN SCFs:
local motif repulsion and whole-chain electronic nonadditivity have opposite
signs in two of four parent frames. The article distinguishes these quantities
from primary molecular G, binding free energy, coupling and ensemble predictions.

The first complete corrected occupied phthalimide cohort now fails its declared
10-to-20 sampling diagnostic. Twenty inputs remain twenty numerical minima;
ensemble G drops by 2.40 kJ/mol and minimum E by 3.88 kJ/mol at the primary
rotor cutoff. All duplicate/rotor scenarios retain the failure. The abstract,
results, conclusion and SI Table S17 report this additional limit, with the
earlier empty-backbone cohorts and primary cycles kept separate. A shared
publication checker reconstructs the geometry/partition analysis for both paper
and slides; the [76 focused tests](../validation/corrected_sampling_2026-11-02/manuscript_shared_checks.xml)
include prefix relocation, numerical tampering and unsupported scope promotion.

The primary counterion test now supplies four checked chloride-contact outcomes
and finite 1-M association G of +15.42 to +22.36 kJ/mol. Geometry/energy duplicate
review places the two signed N0 repairs in one group at all tested thresholds,
leaving three groups. The SI table and supporting slide preserve all outcomes
without assigning population weights. This tests one part of the separated-ion
convention; full salt speciation and exact-graft experimental accuracy remain open.

The unified pool now has 912 exact graphs, preserving all 904 prior graphs and
evaluations. Eight new direct N-linked phenothiazine F/CN pairs remain held and
uncomputed. The unsubstituted baseline has a newly reviewed primary bead route;
three baseline routes pass in total. Route qualification supplies neither a
new exact-graft property nor a full material admission. The current article,
SI and slide 5 use checked frozen registry inputs and verified DOI metadata.

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
