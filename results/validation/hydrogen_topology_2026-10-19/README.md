# Explicit hydrogen connectivity review

The new [shared checker](../../../src/redox/validation/protonation.py) compares
explicit hydrogen attachments while ignoring formal-charge and bond-order
inference from XYZ. It tolerates arbitrary atom order, equivalent sites and ASE
force columns. Six regression cases include acid-to-zwitterion proton transfer,
keto/enol tautomerization, unchanged reordered molecules and a missing H atom.
This is a covalent-graph check, with no claim about equilibrium speciation or
the quality of the electronic model.

The [read-only audit](current_geometry_audit.json) covers all 798 saved geometries
addressed by the 405-row current manifest. It verifies each geometry hash again
at review. These are geometry checks, independent of energy, Hessian or accuracy
admission.

| Geometry source | Matching full graph | H-attachment change | Heavy-graph warning |
| --- | ---: | ---: | ---: |
| DFT | 399 | 3 | 3 |
| UMA | 386 | 3 | 4 |

The three H changes occur in both saved stages for reduced hydroxyanthraquinone,
juglone and D3TaLES 05CZHI. The [case review](case_review.json) independently
cross-checks the published integrity table: all three second reductions already
have `bond_change`, two changed bonds and `intact_bound=False`. They cannot supply
a qualified second electron. This new screen therefore does not increase or
decrease their published usable capacities.

Ferrocene coordination and the explicit viologen ion-pair graphs do not pass
the simple covalent matcher. Their heavy-graph warnings are retained as
unsupported/unresolved here, rather than certified hydrogen matches or evidence
that a ferrocene calculation is wrong. Their existing domain checks remain
separate.

The checker is enforced on the new isolated diaryl-viologen screening outputs,
which include acidic and phenolic groups. Existing primary integrity exclusions
are preserved. No canonical geometry, energy, native thermal partition or
experimental benchmark value was overwritten.
