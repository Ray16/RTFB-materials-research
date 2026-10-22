# Base-assisted grafting pilot

Open the maintained [comparison and status](index.html). The
[vector figure](figures/base_comparison.pdf) and [structured results](current.json)
refresh after each qualified component finishes. Pending quantities remain blank.

**Completed:** all six new components passed QC, with no failed components.
Five candidates have all four base-assisted values (DBU/DIPEA; ions/contact pair).
AQ-2-O/6-OMe remains unresolved because its upstream precursor failed QC.
The two BQ designs and AQ-1-O are **base required; retained for evaluation**.
Their positive HCl-reference values are not a reason to exclude them.
The [selection-update checks](../../reviews/2026-10-20/base_required_selection/summary.json)
record the independent arithmetic review, source hashes and checked slide revision.

This authorized finite pilot compares **DBU and DIPEA in MeCN at 298.15 K** on
the two highlighted BQ ether designs, two AQ ether designs, and NDI/PMDI ethyl
controls. Five reference reactions were qualified at preparation; the AQ
2-O/6-OMe reference remains held and receives no filled-in corrected value.

Only six new species are calculated: each neutral base, protonated base, and
protonated-base/chloride contact pair. Existing qualified candidate/product,
benzyl-chloride and HCl/chloride calculations are reused with hash verification.
The same production DFT/SMD geometry, active single-point and native internal-mode
qRRHO protocol and QC gates apply. Three geometry seeds are evaluated per new
species; thermal treatment applies to the selected electronic-energy minimum.

For an attachment precursor `H–Nu` and local capped site `P–CH2–Cl`, compare:

```text
Reference:       P–CH2–Cl + H–Nu     → P–CH2–Nu + HCl
Separated ions:  P–CH2–Cl + H–Nu + B → P–CH2–Nu + BH+ + Cl−
Contact pair:    P–CH2–Cl + H–Nu + B → P–CH2–Nu + [BH+···Cl−]
```

The correction to the reference is the free energy of base capture of HCl.
The paired model has one fewer independently translating solute. Its explicit
1-atm-to-1-M reaction correction is approximately −7.93 kJ/mol at 298.15 K;
that correction cancels for the separated-ion and original reference reactions.
The paired/free-ion thermodynamic cycle is checked for closure. Neither model
assumes precipitation of a solid salt or a measured ion-pair population.

**Scientific limits:** these are thermodynamic sensitivity calculations, not
barrier calculations, experimental yield predictions or synthesis admissions.
The base precedents involve AQ-carboxylate/DBU in NMP and
ferrocene-carboxylate/DIPEA in DMF; applying them to these exact phenolic motifs
in MeCN is a proposed test. Competing base alkylation, quinone side reactions,
polymer diffusion and attachment selectivity remain unresolved. Contact-pair
entropy and conformer sampling are approximate. Favorable modeled base-assisted
reactions support retaining these candidates; experimental selectivity and yield
are separate requirements. See
[practical chemistry and sources](../../../docs/PRACTICAL_GRAFTING.md).

The frozen [plan](../campaigns/2026-10-19_base_assisted_pilot/plan.json) records
inputs, producers, stopping tolerances and references. Raw calculations are at
`calcs/grafting/2026-10-19_base_assisted_pilot/<component>/rest/`.
Queue ownership, node preflight, launches, actual-start verification, preserved
failures and component receipts live in the campaign directory. The queue makes
one actual attempt per species and preserves existing outputs.

Entrypoint: `scripts/pipeline/base_assisted_grafting.py`; use `collect` to verify
receipts and refresh the report without launching calculations. This pilot remains a finite batch; subsequent validation runs are recorded separately.

A [scoped external acid-capture check](../../validation/acid_base_external_2026-10-22/README.md)
compares the DBU free-ion term with author-tabulated MeCN pKa values. A
fixed-reference hybrid sensitivity retains both BQ signs. It does not calibrate
the full reaction, and no exact DIPEA anchor was substituted from a related amine.
