# Figure 3 — the reliability-gating funnel

**File:** `fig3_qc_funnel.png` (300 dpi). PNG only, per the figure style contract.

**Figure 3. Every reported candidate passes explicit, automated QC; the identity
registry narrows to the scored pool through sequential gates.**
The funnel shows four tiers: registered graphs → state-complete → QC-clean →
scored candidates. A candidate advances only if it is state-complete (all declared
redox states built and converged), then QC-clean (declared-vs-realized charge and
spin, ⟨S²⟩/spin-contamination, geometry-protocol currency, full state coverage,
and cross-point consistency all pass; `quality/readiness.py --strict`,
`quality/spins.py`, `quality/stability.py`). **Candidates that fail a gate are
held, not dropped** (left annotation): they are carried on a separate
partial-information list and can be dominated but never silently removed or allowed
to dominate. The scored tier is then ranked on Pareto fronts (Figure 5).

**Honesty about the counts (critical).** The funnel spans **three grains**, as
annotated on the left:
- **Real / verified:** *registered graphs = 912* (**graph grain**; identity
  registry, counting molecular graphs incl. per-state and held/uncomputed graphs;
  `docs/PUBLICATION_READINESS.md`); *state-complete = 450 of 456 declared states*
  (**state grain**; from a pool-wide `redox.quality.readiness` run, 2026-11-17 —
  450 states have a converged geometry + active single-point, only 6 are
  `missing_geometry_or_result`); and *scored candidates = 91* (**candidate grain**;
  the frozen manuscript pool). The maintained
  `results/validation/full_pool_scored_legacyG_MeCN_TBAPF6.csv` currently carries
  **92** rows — the extra row is one phenazine catholyte appended by the live
  `expand_catholyte` campaign; the frozen pool is 91.
- **Placeholder (one tier, flagged):** *QC-clean* is drawn hatched and labelled
  "placeholder" because its final gate includes **validated native-thermal free
  energies**, whose certificates are path-locked to the original
  `/nfs/.../0_redox` checkout. Read from this archive-symlinked tree, the readiness
  run reports 351 *spurious* `missing_validated_native_thermal` states, so the naive
  "103 ready" is a **relocation artifact** and is deliberately not shown. A faithful
  QC-clean count needs `redox.quality.relocated_native` (or a checkout at the
  original path). Bar widths are a fixed decreasing schedule for legibility and do
  **not** encode counts proportionally.

**Provenance / reproducibility.** Scored count read live from the pool CSV; 912 and
the gate definitions are fixed from the cited sources; state-complete (450/456) is
from the dated `redox.quality.readiness` run documented in the generator. Generated
by `scripts/plotting/manuscript/fig3_qc_funnel.py` (CPU only).
