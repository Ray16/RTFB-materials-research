# Source-bound native spin equilibrium

CN-PMDI's current singlet and triplet dianions are independently qualified own
minima. The two-spin equilibrium gives a **32.16% triplet population** at 298.15 K
and a **−0.009968 eV** correction relative to the singlet surface alone. Each
component native free energy already includes spin degeneracy; the mixture does
not multiply it again. This is a stated two-spin equilibrium model, not a kinetic
or experimental validation.

`redox.validation.spin_ensemble` reads current active energies and corrected native
partitions, enforces matching chemistry/charge/Hamiltonian and binds geometry,
energy, raw Hessian and partition artifacts by hash. The alternative's original
seed identity is preserved: the singlet was repaired later, so this is a new
comparison of two fully qualified own minima, not reuse of the earlier electronic
gap report. Wrong chemistry, unqualified minima or changed inputs are refused.

The independent reader in `scripts/validation/audit/spin_ensemble.py` reconstructs
identities and uses a separately factored partition equation. The prepared policy
passed **38 independent checks**, and an isolated preview passed **eight integration
checks**, including candidate, scorecard and independent-state admission. The
shared hook was activated only after both running methyl-NDI branches had passed
their frozen-source guards; exact before/candidate bytes, AST equivalence and live
process observations are retained. Only the analysis wrapper changed.

[Selection receipt](selection_receipt.json) records transactional selection under
writer locks and agreement between both real readers. Invalid selections restore
the exact previous policy bytes and retain failed candidates. The selected file
is `calcs/dft/search_pmdi_1864c2c74267/red2/spin_ensemble_policy.json`.

The full software suite passes **436 tests**. Integrity checks now cover every
included spin-component pair; both CN-PMDI pairs are bound with unchanged bond
graphs. Gas relaxation remains a flagged diagnostic, and solution kinetics stay
unvalidated. The global audit independently rebuilds distance graphs instead of
copying bond-change counts from the generated table. The current checked candidate export and slides include the equilibrium second
potential. The bound 25/50/100 cm⁻¹ diagnostic changes that potential by only
2.70 mV; 50 cm⁻¹ remains primary without fitting.

Raw arrays remain on shared storage, referenced by the policy and review hashes.
Historical reports and raw component records are preserved without relabeling.
