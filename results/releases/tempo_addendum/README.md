# TEMPO local-graft raw-state addendum

[Archive and hashes](current.json) · [Relocated replay](relocation.json) ·
[Adverse-input tests](negative_tests.xml) · [Scientific result](../../grafting/tempo_exact_support_2026-10-26/README.md)

This local addendum contains the four cached states used in the TEMPO
local-graft oxidation potential: radical and oxidized TEMPO plus neutral and
oxidized ferrocene. It includes their geometries, active solution-energy
records, native Hessians and projected thermal records. The Fc states overlap
with earlier packages and must not be added as new unique molecules.

After extraction outside the original checkout, an isolated standard-library
replay reconstructed the potential at 25, 50 and 100 cm⁻¹. A separate locked
NumPy 2.2.6 replay reconstructed all four Hessian mode sets and nuclear sums.
Python audit hooks blocked opens under the original project; the package-relative
inventory and hashes were verified before replay. Six focused tests include
rejection of changed or rehashed energy/Hessian records, path escape, duplicate
identity and promotion of unsupported claims.

The archive is on shared storage at the path in [current.json](current.json).
After extraction, run `python3 -I -S /path/to/addendum/reproduce.py --data
/path/to/addendum`; the package README describes the optional `--full-native`
replay with the locked NumPy environment. Cached replay does not repeat SCF,
stationarity, spin or experimental measurements, establish grafting free energy
or polymer interaction, or demonstrate material eligibility. This addendum is
not yet part of the consolidated manuscript bundle and is not a public deposit;
data licensing, author metadata and DOI remain open.
