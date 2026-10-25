# Checked raw numerical-control addendum

[Current package and hashes](current.json) describe a 1.13 MB compressed
addendum on shared storage, complementary to the v5 and first-site packages.
Its 510 regular files include eight four-atom geometry/Hessian records, six
fixed-geometry Fc electronic records, code and environment snapshots. The 392
stored SCF evaluations comprise eight baselines and 384 displaced evaluations.
These are calculation records with overlapping chemical identities.

[Archive extraction and relocated replay](relocation.json) verify every member
and run two independent modes while original project reads are blocked:
standard-library electronic/component arithmetic, and a newly installed,
hash-locked NumPy 2.2.6 CPU replay of forces, sampled Hessian columns, projected
modes and nuclear sums. No production redox, PySCF, CUDA or xTB module is loaded
by the replay. [Full numerical results](full_replay.json) retain per-state errors.

The small [self-contained example](../../../examples/control_replay/README.md)
ships the cached raw records with the repository. Bulk archives remain under
`calcs/releases/2026-10-25_numerical_control_addendum_checked{,.tar.zst}`.
The frozen source revision, code hashes and original data locations are retained
in the package. Cached replay does not rerun quantum calculations, validate
solvation accuracy, establish polymer G or show ensemble convergence. Formal
public archival deposition, data license and DOI remain pending.
