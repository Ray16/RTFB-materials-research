# Grafting reaction free energies

All result categories: [results index](../README.md).

**Open the [visual results browser](index.html)** for the maintained ΔG ranking,
reaction schemes, coverage and all 75 exact molecular structures. The same files
refresh as qualified calculations finish; reload the browser to see updates.

**Stable result table:** [current.csv](current.csv). **Progress:** [status.json](status.json).
The inventory covers all 75 current candidates. No grafting free energy is reported
until all reaction components have matched, qualified free energies.

The first screening solvent is acetonitrile at 298.15 K. The substrate is
4-methylbenzyl chloride, a capped local model of a Merrifield attachment site.
These calculations model local substitution thermodynamics, not whole-bead
chemical potentials, activation barriers, grafting yields or site accessibility.

The balanced reference reactions depend on the attachment chemistry:

- N–H/O–H precursors: `site–Cl + Nu–H -> site–Nu + HCl`.
- N quaternization: `site–Cl + Nu -> site–Nu(+) + Cl(-)` (with all other charges retained).

The complete atom/charge-balanced schemes and exact molecular identities are in
[reactions.json](reactions.json); shared reagents are in [species.json](species.json).
The neutral HCl reference is explicit. A base was not specified, so the free-energy
contribution for HCl capture/proton transfer is kept separate and pending. A common
base correction can later be applied to the N–H/O–H reactions under a defined
chemical model. Ionic strength, ion pairing and counterion exchange are additional
assumptions; cross-reaction-class rankings must retain these limitations.

The reference reactions have two reactant and two product solute species, so the
uniform ideal 1-atm-to-1-M standard-state correction cancels in their reaction sums.
Atomic chloride needs its monatomic translation/spin partition, not an artificial
vibrational Hessian. Molecular species use the existing matched DFT/SMD and native
internal-mode thermochemistry. A missing thermal correction is never zero-filled.

At preparation, 67 existing product free energies can be reused; eight products
need further qualified calculations. There are 72 unique precursor graphs and
three shared references. See the status file for subsequent progress. The
reference-reaction columns and any future base-coupled columns are distinct.

## Submitted calculation campaign

The [frozen campaign plan](campaigns/2026-10-19_mecn_reference/plan.json) contains
83 independent component calculations: 72 precursors, three shared references
and eight missing product states. The remaining 67 product free energies are
reused with source hashes. Eight H100 workers and two V100 workers have been
submitted through the system reservation gate. Their actual processes, GPU
residency, native thread caps and disjoint CPUs are recorded in
[start verification](campaigns/2026-10-19_mecn_reference/start_verification.json).
That is launch evidence, separate from qualified component and reaction counts.

Raw calculation outputs remain on shared storage at
`calcs/grafting/2026-10-19_mecn_reference/<component_id>/rest/` from the project root.
The campaign directory retains per-item claims, owner identities, logs, completion
receipts and preserved QC failures. One actual attempt is allowed per item;
failures require review before another attempt. The molecular protocol matches
the qualified production method. Atomic chloride has a separate, independently
checked monatomic partition and no artificial vibrational modes.

The finite GPU workers collect qualified components after each item. A reaction
gets a ΔG only when every component is available and unchanged. They refresh:

- [ΔG ranking](figures/deltaG_ranking.pdf), separated by reaction class.
- [Coverage](figures/coverage.pdf), including pending and QC-held reactions.
- [Balanced attachment schemes](figures/reaction_schemes.pdf), with vector molecular drawings.
- [Searchable structures and results](index.html), [JSON feed](current.json) and CSV.

Manual read-only refresh, from the project root:

```bash
env PYTHONPATH=src OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python \
  scripts/pipeline/collect_grafting_results.py
```

This refresh submits no calculations. The figures retain pending fields rather
than inventing values, and reference-reaction rankings do not include an
unspecified base contribution. The prepared oligomer inventory remains a separate
[polymer-model task](../polymer_models/README.md).
