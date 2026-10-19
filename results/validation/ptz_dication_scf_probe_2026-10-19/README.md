# Frozen-geometry electronic diagnostic for the DMeOEPT dication

The restricted active-protocol energy is reproduced within 5.2e-9 eV. Four
independent zero-net-spin UKS runs start from opposite HOMO/LUMO rotations in the
alpha/beta occupied orbitals (0, 0.2, 0.5, 0.8 radians), preserving electron counts.
All recover the restricted energy within 2.4e-7 eV and have recorded S-squared 0.
No lower solution was found under these seeds. This is an orbital-seed diagnostic,
not a formal stability-eigenvalue or exhaustive electronic-state proof.

All calculations use the frozen selected geometry, active functional/basis/SMD
settings and reserved H100s with disjoint two-CPU sets. The baseline precedes the
four parallel seed probes. The initially attempted live-start observation ran
before every dependency probe had created its report; no job was restarted.
`verification.json` records outputs that finished before the later observation.

No geometry, thermal correction, reference or canonical calculation is changed.
Unrestricted zero-spin solutions are diagnostic and cannot be adopted as pure
singlets. The measured second-wave discrepancy is not explained or fitted away.
Raw orbitals remain in shared storage with hashes; small source/results are
tracked for reproducibility.
