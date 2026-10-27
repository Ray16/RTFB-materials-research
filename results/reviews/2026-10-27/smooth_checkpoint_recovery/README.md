# Four preserved geometry restarts

Three source attempts stopped on HDF5 checkpoint locking; the fourth was interrupted by foreign GPU residency. All source producer and gate owners were confirmed stopped before isolated restart preparation. Three inputs use the last SCF-converged, explicitly unqualified geometry evaluation; the first source has no such frame and retains its original input. Each exact graph and stereochemistry was audited again.

The new worker retries only EAGAIN/EWOULDBLOCK errors whose messages identify a checkpoint lock, with six bounded delays totaling 7.8 seconds. Locking remains enabled; other errors and exhausted retries remain failures. No old output is replaced. A real CPU HDF5 round trip and adverse-input tests pass. Every geometry, force, resource and model check remains in the original scientific worker.

The [startup snapshot](startup_snapshot.json) distinguishes dispatch from scientific-stage status. The separate independent reviewer is declared in [the launch record](launch.json). Original and new calculation paths remain bound in each case under `results/polymer_models/campaigns/2026-10-27_smooth_checkpoint_recovery/`; raw outputs use the matching `calcs/polymer_models/2026-10-27_smooth_checkpoint_recovery/` namespace. No converged geometry, energy or ensemble result is inferred from starting these recoveries.
