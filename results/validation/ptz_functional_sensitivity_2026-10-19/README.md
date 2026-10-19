# Predeclared fixed-geometry electronic-method sensitivity

Actual active-level relaxation changes DMeOEPT's electronic charging gap by
-0.027134 eV, insufficient to remove the original free-energy discrepancy.
This bounded follow-up uses all three same frozen private geometries for three
predeclared electronic approximations: current wB97M-V/VV10, PBE0-D4 and B3LYP-D4.
All use def2-TZVPD, RI-J, SMD/MeCN, grid level 3 and the same tighter SCF settings.
Neutral singlet, radical-cation doublet and dication singlet are held fixed.

The [plan](plan.json) records the choices before results. There are nine points,
including three reproductions of the prior WB energies. The alternatives differ
in functional and dispersion together; this is not an exact-exchange-only test.
They are not relaxed on their own surfaces, no thermal terms are computed, and
no matched alternative Fc reference is required for the electronic charging gap.
A gap is not a new standard potential or a complete free-energy prediction.

The [preflight](lambda13_preflight.json), [launch](launch.json) and actual-start
observations record atomic GPU reservations, disjoint CPU pairs and native pool
caps. All currently free suitable GPUs are used up to the nine independent points.
Each point gets one actual attempt; failures preserve outputs. No production
method is selected or canonical value changed on the basis of this single parent.

The installed [NLC parameter selection check](NLC_parameter_source_review.json)
confirms that energy and gradient branches prioritize the NLC-containing XC
functional. Standalone VV10 coefficients differ, but that does not show a
production parameter error. The initial dependency preflight asked for the wrong
module name (`dftd4`); the installed interface is `pyscf.dispersion.dftd4`. That
was corrected before any scientific attempt, without changing the environment.

After all points finish, `review.py` independently checks their input identities,
SCF/spin diagnostics, control reproduction and three-state energy combinations.
Method spread is sensitivity, not an empirically calibrated uncertainty or a
reason to choose whichever method matches one experiment most closely.

All nine points finished successfully. The [independent review](independent_review.json)
passes 71 checks, including the complete declared method/state matrix, executed
task identities, frozen-geometry identities, energy units, spin diagnostics and
three WB control reproductions. The [terminal verification](terminal_verification.json)
checks the actual original process identities and all eight controller exits;
the original reservations have ended.

| Electronic approximation | Charging gap (eV) | Change from WB (eV) |
| --- | ---: | ---: |
| wB97M-V | 0.886511 | 0.000000 |
| PBE0-D4 | 0.954183 | +0.067671 |
| B3LYP-D4 | 0.908061 | +0.021550 |

These sampled alternatives raise the electronic gap at the same geometries.
They do not supply an explanation for the original positive second-wave residual.
The 0.599-V measured wave spacing is a free-energy quantity; it is not a target
electronic gap. Neither an electronic-method winner nor a fitted correction is
inferred. Own-surface geometry, thermal, conformer and electrolyte effects remain
outside this experiment. The phenothiazine selection hold remains in force.

The [preservation review](preservation_review.json) compares canonical parent
records with the previously published frozen parent inputs, and published
candidate files and numeric property tables with the preceding Git checkpoint.
This private diagnostic does not promote any new property.
