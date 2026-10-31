# Exact viologen: cached redox and grafting replay

This example is maintained on the public `refactor/clean-workflow` branch.

This small CPU example uses the same exact capped graft for the reduction
ladder and local quaternization. It runs from a checkout with **NumPy only**:

```bash
python3.11 -m venv .venv-replay
.venv-replay/bin/python -m pip install --require-hashes --no-deps -r examples/exact_viologen_replay/data/requirements.lock.txt
.venv-replay/bin/python examples/exact_viologen_replay/replay.py replay
```

The locked binary requirement is tested on CPython 3.11/Linux x86_64. No redox
package installation, GPU, RDKit, PySCF, ASE or xTB library is required. Original
absolute paths in state records are provenance metadata; replay opens only
the supplied package and its numerical runtime.

At the 50 cm-1 rotor setting, the cached states reconstruct:

| Quantity | Value |
| --- | ---: |
| First reduction vs matched Fc | -0.6077915 V |
| Second reduction vs matched Fc | -1.2080017 V |
| Local separated-ion quaternization | -11.49333 kJ/mol |
| Monomer disproportionation | +57.91148 kJ/mol |

Ten states cover three declared graft redox states, three fixed precursor
starts, two Fc states, the local benzyl-chloride site and chloride. Nine
Cartesian Hessians are projected and their frequencies/nuclear terms
recomputed. Six recent states also retain arrays that independently reconstruct
two finite-difference columns at two steps. The inherited Fc/site FD summaries
remain metadata-only. The chloride partition is recomputed analytically.

The example preserves MeCN, 298.15 K, declared spin states, the matched spin-only
Fc convention and the lowest-precursor-G selection at 50 cm-1. All displayed
balanced reactions have zero net solute count, so the source 1-atm nuclear
convention has zero net conversion to a 1-M reaction convention. No ion-pairing,
counterion-exchange, reflux-temperature or polymer correction is inferred.
Rotor-setting sensitivity is not a statistical uncertainty interval.

The data directory is a hash-locked generated snapshot, including its as-used
runner/exporter. Preserve its bytes. The maintained runner is `replay.py`, and
the exporter is `scripts/release/molecular_example.py`. Re-export into a fresh
directory when changing its inputs or implementation.

This verifies cached numerical reconstruction. It does not rerun SCF, validate
all historical executable/scratch records, establish global conformer/spin
ordering, or admit a complete polymer material. Full archival deposition/DOI
and final deposition metadata remain separate publication tasks.
