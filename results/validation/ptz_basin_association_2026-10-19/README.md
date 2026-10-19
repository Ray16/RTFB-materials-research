# Reciprocal local-basin association for the PTZ lead

Both charge-state controls optimize from the opposite state's geometry using the
same native geometry-level r2SCAN-D4/def2-SVP/SMD Hamiltonian and full grid response.
All five strict optimizer criteria, final SCF restoration, graph identity and
endpoint hashes are checked. Production coordinates are never replaced.

The neutral returns to the selected geometry with proper chemical-symmetry RMSD
0.000135 Å and energy difference −2.05 micro-eV. The radical cation returns within
0.00178 Å and +23.4 micro-eV. Both pass the predeclared local-pair criteria
(0.05 Å, 1 meV), with **10 checks and zero failures**. A reflection fit is recorded
separately and is not used to claim geometric return.

The initial tether angles alone differed by about 171° on one bond; chemical
symmetry makes angle proximity insufficient. Reciprocal optimization and proper
rotation/automorphism matching provide a stronger local association. This does
not prove a complete conformer ensemble or a polymer constraint.

[Derived descriptor qualification](descriptor_qualification.json) binds this review,
the selected native minima and the independently checked six cached SCF densities.
Directional electronic solution relaxations remain 1.1157 / 1.1373 eV (mean 1.1265).
They are composite active-SP descriptors at geometry-level minima. Active-SP
stationarity, finite-temperature barriers and experimental kinetic accuracy remain
unvalidated; primary rates are not promoted. Do not add separate gas/outer terms.
The raw earlier diagnostic and its unresolved-conformer flag remain preserved;
the new qualification is an explicit derived evidence layer.
