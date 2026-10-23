# SMD temporary-buffer execution controls

The 139-atom PMDI resource pilot exhausted H100 memory during dense SMD matrix
assembly, before geometry optimization. This resource error does not reject its
geometry. The original failed pilot remains intact.

A run-local adapter releases the DA scratch matrix, reuses DAS storage for K,
and releases the unused geometry scanner/SCF before the larger single point.
Installed libraries, surface/grid, solvent descriptors, functionals, basis sets,
precision and physical acceptance thresholds are unchanged. Source snapshots
and upstream license notices accompany each execution.

Both [singlet](controls/1/status.json) and [doublet](controls/2/status.json) controls
pass full original/compact comparisons. Maximum energy difference is 1.4e-11 eV;
coordinates are identical; gradient and Hessian differences are below 4.0e-14
and 9.4e-12 in their respective atomic units. This qualifies numerical equivalence
on the controls, not memory capacity at the larger geometry.

The initial launch-observation timeout did not stop the detached controls.
[Actual PID/start verification](actual_start_verification.json) confirms their
owned starts; they were not resubmitted. Full-size resource qualification runs
in a separate recovery case linked from the PMDI campaign, preserving the
original failure and unchanged 80-GiB pilot safety budget.
