# Alternate-spin helper optimizer correction

The helper minimum repair used geomeTRIC `optimize()`, which returns a geometry
without retaining the convergence flag. It now uses the same `kernel()` interface
as the main geometry driver and requires both optimizer and final-SCF convergence
before writing a replacement geometry. Full grid-response gradients and existing
optimization/physical acceptance thresholds are unchanged.

Refined metadata records the actual optimized path and source geometry/result
hashes. Inherited direct-energy geometry/protocol declarations are cleared;
active energies still require their own addressed SP records.

Five new regression cases cover successful convergence and separate optimizer/SCF
failures, including the actual helper failure path preserving its input files.
The full suite passes 461 tests (`software_checks.log`). These tests establish the
implementation contract; they do not constitute an independent quantum replay.

`source_epoch.json` and `producer_before.py` preserve the earlier producer. All
12 primary parent workers had passed their startup source guards before this
change. Their geometry, SP and native thermal numerical producers are unchanged.
The new helper will be used for subsequent source-addressed spin calculations.
