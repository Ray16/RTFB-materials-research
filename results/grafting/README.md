# Grafting reaction free energies

All result categories: [results index](../README.md).

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
