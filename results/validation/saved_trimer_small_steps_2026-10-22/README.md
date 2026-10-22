# Fine-step validation of the saved trimer force

The fixed 0.0005/0.001 bohr energy differences at the optimizer-returned trimer disagree with its analytic gradient by about 1.47e-4 Ha/bohr. They remain preserved as a failed check. The earlier evaluation-42 diagnostic showed strong deterministic step dependence and fine-step agreement; it did not qualify the final geometry.

This separate experiment tests the actual optimizer-returned geometry at 1e-4, 1e-5 and 5e-6 bohr, with three repeated baseline solves. [plan.json](plan.json) predeclares the same 1e-6 Ha/bohr analytic-versus-FD bound, an additional 1e-6 finest-pair agreement bound, and a 1e-12 Ha repeat-energy spread limit. The minimum geometry/SP identities remain fixed. This is a numerical step-convergence test; an acceptable result will still require explicit independent review and pipeline integration rather than silently changing the old check.

Raw evaluation records and durable source/ownership status: `calcs/validation/saved_trimer_small_steps_2026-10-22/`. The existing checkpointed diagnostic is launched through a fresh reserved H100. Running jobs and emitted files do not establish scientific acceptance; no Hessian, G, electron-transfer rate or synthesis admission is claimed.

The diagnostic completed, but the finest two steps still disagree with analytic gradients by about 1.0e-4 to 1.25e-4 Ha/bohr, above the unchanged 1e-6 criterion. The repeated baseline energy spread is 3.64e-12 Ha, also above the predeclared 1e-12 requirement. No geometry was promoted. Further numerical-surface work is required, or this original trimer case must remain excluded.
