# JCIM main-text rewrite, workflow contracts and viologen recovery

[Rewrite plan](PLAN.md) · [JCIM structure/build checks](submission_check.json)
· [Context-bound grafting check](quaternization_context_check.json)
· [Workflow catalogue](workflow_catalog.json) · [Focused tests](tests.xml)

The main manuscript was rewritten from a new outline for a full-length JCIM
Article. The title is nine words and the abstract has four sentences. The
scientific argument is context-preserving assessment across exact molecular,
reaction and polymer models, with the novelty distinguished from established
thermodynamic identities and multistep screening. The Methods define permitted
property combinations and the Results connect benchmark applicability to the
viologen, TEMPO, PMDI and phthalimide examples. The new conclusion states the
supported scope without declaring a complete material shortlist. All numerical
results remain tied to the same immutable scientific snapshot; the SI's
technical evidence is retained, with its title aligned to the new main text.

The checked PDFs contain 18 main-text pages including references/TOC and 30 SI
pages. Both compile without undefined references/citations or overfull boxes.
The main manuscript's 18-page layout was reviewed in a contact sheet, with pages
1, 3, 8 and 13 examined at readable resolution. The recompiled standalone
reviewer package contains 33 files and reproduces the text of both local PDFs.
The original main text, SI and build/pointer records are preserved under
`calcs/validation/jcim_rewrite_2026-11-01/before/`; rendered pages, compiler logs
and the standalone ZIP are also retained under the task's calcs directory.
`manuscript/build_check.json` and `manuscript/overleaf_export.json` bind the
published artifacts. No Overleaf remote is configured and the restoration
record says the draft is not connected; no Overleaf synchronization is claimed.
A cover-letter draft is supplied without inventing author approval or originality
statements. Authorship, full data/software archival release and remaining
scientific validation must be resolved before final submission.

## Workflow and property changes

`redox workflows` classifies maintained orchestration, design, analysis,
calculation stages and specialized controls/recovery. It lists the available
public registry rather than implementing another dispatcher. The shared
resource layer discovers the reservation executable and lock directory through
explicit configuration/runtime discovery, while CPU-only imports do not require
a GPU service. The reservation policy and physical admission thresholds are
unchanged. See `docs/WORKFLOWS.md` for configuration and scope.

The shared checked-arithmetic contract retains quantity/unit, exact component
identity and charge/spin, electronic/thermal method, solvent/model, temperature,
standard state, reference, QC and source binding. The quaternization publisher
uses this contract at every cutoff and precursor seed after normal raw-record
review. The new publication at
`results/grafting/viologen/2026-11-01_context_checked/` reproduces the previous
three cutoff values exactly, including −11.493329652 kJ/mol at 50 cm−1, and
retains the new metadata. Its source-digest fingerprint survives relocation of
producer code into unique frozen files. The old pure arithmetic API explicitly
reports an unverified context when callers omit context records; existing
scientific evidence is not silently requalified. This is one integration of a
shared boundary, not a completed migration of every historical property reader.

Sixty-five focused tests pass, covering context mismatches, nonfinite values,
missing spin/QC, equivalent checked arithmetic, authoritative runtime discovery,
CLI registry consistency, publication arithmetic and charged-native contracts.
The production code and paper were checked in a working tree with pre-existing
changes; those unrelated changes are preserved. The old product-review snapshot
and a numerical-only review were correctly rejected during reassembly; the
source-complete checked recovery receipt was used for the successful publication.

## Charged viologen calculation

The earlier site-2 native pilot had stopped at the shared-storage margin. Its
exact PID was absent, its checkpoints were preserved, and current free storage
was approximately 291 GiB. The recovery is a new namespace, not an overwrite:

`/nfs/lambda_stor_01/homes/rzhu/redox_candidate_continuations/2026-11-01_viologen_all_same_occ_2_native_storage_recovery/`

The preparation binds the independently checked 114-atom, charge +2, singlet
stationary geometry and frozen worker sources. Twenty charged-native protocol
checks passed before launch. Node dependency preflight passed on lambda4, and
`gpu_reserve run` owns GPU 0. The actual worker PID/start identity, CPU affinity
[0,1], six thread caps of two, one resident process on that GPU and the initial
SCF progress are recorded separately in `native_start_verification.json`.
The SCF progressed to cycle 28 and the calculation entered its first-derivative
stage. Bounds remain 80 GiB host RAM, four GiB durable output, eight hours and
a 32-GiB live shared-storage margin. Only the first charged native pilot is
admitted until its independent physical/resource review passes; other free GPUs
do not remove this dependency gate.

This is a verified calculation start, not a completed Hessian, admitted G,
occupied free-energy cycle or ensemble. The result remains absent from the paper.
The model is the separately identified smooth-surface MeCN branch and must not
be mixed into primary-model molecular/grafting sums. Inspect current ownership
and outputs before any future continuation rather than inferring liveness from
this dated record.
