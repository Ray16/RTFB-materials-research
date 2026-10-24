# CPU electronic-structure equivalence controls

[Both stored spin controls](comparison.json) pass CPU/GPU comparisons for active gas/SMD energies and geometry-level r2SCAN-D4 energies/gradients, with the same basis, solvent and model settings. The [execution adapter](../../../scripts/pipeline/cpu_quantum_execution.py) changes resource limits, SCF checkpoints and metadata inspection, and prevents CUDA imports in the CPU-only worker. Physical models and thresholds are unchanged.

Raw SCF checkpoints and outputs are at `calcs/validation/CPU_quantum_2026-10-23/`. The independently reserved CPU worker used sixteen disjoint cores on lambda0. These two rigid controls do not qualify the 139-atom polymer resource footprint, its optimization, active parent or native free energy. A large exact-geometry pilot is still required before a production CPU fallback. The existing GPU calculations remain separate.
