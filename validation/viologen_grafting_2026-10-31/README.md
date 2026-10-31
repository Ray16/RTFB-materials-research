# Exact viologen quaternization

All three precursor states pass independent native/resource review. The local
cycle gives **−11.49 kJ/mol** at the primary setting, with a 1.50 kJ/mol
precursor-seed span and −13.30 to −11.49 kJ/mol rotor sensitivity. The
[summary](summary.json), [reaction](reaction.json) and [checked seed table](QC_table.tex)
record the completed finite-model result.

This campaign computes the missing mono-N-methyl-bipyridinium precursor for
the route-admitted baseline enumerated viologen. The recipe is
[`config/viologen_grafting_component.json`](../../config/viologen_grafting_component.json).
The precursor is a reaction component, not a new redox candidate.

The defined reaction is a monocationic precursor plus neutral capped benzyl
chloride giving the dicationic capped graft plus chloride. It is atom/charge
balanced, contains two separated species on each side, and has no HCl/base term
or net 1-atm-to-1-M correction. MeCN and 298.15 K are the computational
conditions. Counterion association, PF6 exchange, polymer accessibility,
conversion and the experimental reflux temperature are outside this model.

Three fixed RDKit/MMFF seeds (20261031) enter standard UMA S1.2.1 warm starts,
then r2SCAN-D4/SVP SMD geometry optimization, uniform wB97M-V/TZVPD energies,
and native DFT Hessians with projected qRRHO at 25/50/100 cm-1. No UMA energy
enters the grafting sum. The selected precursor is the lowest G at 50 cm-1
among these fixed seeds, with every seed and setting reported separately.
This is finite sampling and does not prove global conformer convergence.

Raw namespaces are on durable shared storage:

- `redox_candidate_continuations/2026-10-31_viologen_precursor/conf000/`:
  the completed V100 precursor, including its original frozen worker.
- `redox_candidate_continuations/2026-10-31_viologen_precursor_H100_checked/`:
  conformers 1 and 2. Before each warm start, three fixed-geometry standard-UMA
  energy/force evaluations must reproduce the V100 reference within 1e-4 eV
  and 1e-4 eV/angstrom. These are standard-model checks, not Nitro qualification.
- Original conformers 1 and 2 in the first preparation were never started:
  host RAM promises refused the placements before GPU reservation. They remain
  preserved and are superseded by the explicit H100 namespaces.

All paths above are relative to `/nfs/lambda_stor_01/homes/rzhu/`.
Placement logs and independent reviews are in
[`results/validation/viologen_grafting_2026-10-31/`](../../results/validation/viologen_grafting_2026-10-31).
Raw SCF checkpoints, optimizer returns, settings, forces and Hessians remain
with each private case. The execution uses shared GPU reservations, disjoint
CPU pairs, native thread caps, exclusive-residency checks, a prospective 24-GiB
host bound and a 32-GiB available-host-memory floor.

The independent review reconstructs the native Hessian's sampled columns,
projected modes and nuclear terms and checks density/spin and the complete
resource/runtime chain. The shared site/chloride records come from the existing
14-state portable neighbor package; `shared_reference_replay.json` reconstructs
all 13 cached Hessians and the one monatomic partition. The fresh precursor
uses its own thermal corrections throughout.

The versioned commands are:

```bash
python -m redox.workflows.reaction_component_namespace prepare --help
python -m redox.workflows.reaction_component_namespace probe --help
python -m redox.workflows.reaction_component_namespace submit --help
python scripts/pipeline/review_reaction_components.py --help
python scripts/pipeline/publish_quaternization.py --help
```

`publish_quaternization.py` requires independent native/resource review of all
three requested precursor seeds, the exact dicationic product, matching active
and native protocols, shared-reference numerical replay, and a balanced
quaternization. Missing thermochemistry cannot be filled with electronic
energies. Its output remains a local finite-sampling model and does not admit
the candidate as a complete polymer material.

The separate `polymer_derivative_failure.json` independently reconstructs a
failed occupied-polymer check from saved gradients. The 1.2138% column fails
the unchanged 1% limit. Approximately constant gradient-difference residuals
at other columns motivate checking quadrature/cavity switching; their origin
has not yet been established. This diagnostic does not qualify that state.
