# Complete, transactional spin-helper repair adoption

The alternate-spin minimum repair previously replaced helper geometry/result files
before validating the replacement at its final location, and copied only `sp` and
`native_thermal`, omitting the corrected `projected_thermal` partition.

It now calls the already qualified `native_thermal_repair.install_repaired_state`
transaction. The adapter binds the original helper geometry before optimization,
carries all three caches and old addressed histories, verifies the private install,
and rolls back any failed post-install admission. Original files remain archived
with hashes. The installer producer and optimizer producer have separate records.

Five regression cases exercise the real helper adapter: complete partition and
old evidence preservation; missing partition; private validation failure; final
validation failure with rollback; and a source changed during refinement. All
466 software tests pass (`software_checks_final.log`). No physical thresholds,
Hamiltonian or numerical producer for the running primary panel changed.

`producer_before.py` preserves the e910d65 helper implementation. Primary workers
had already passed their source guards. Subsequent alternate-spin jobs use the
new helper source epoch, with source hashes frozen by their own run plan.
