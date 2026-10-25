# Checked first-site polymer milestone, 25 October

[Status snapshot](summary.json) · [Chemical results](../pmdi_first_site_grafting_checked/README.md) · [Software checks](../../../publication_readiness/software_checks.json) · [Current paper](../../../../manuscript/manuscript.pdf)

The exact PMDI example now connects checked single-site polymer free energies to balanced first-grafting reactions and a base-canceling local-versus-polymer comparison. The paper, SI, slide table and self-contained upload archive contain these verified results. Sampling limits and the missing paired-state free energy remain explicit.

The latest large V100 calculation saved durable partial Hessian and H1 checkpoints, then stopped when the orbital-response Krylov solver diverged. It remains unqualified and published no free energy. Complete small singlet/doublet controls validate the new bounded derivatives and compatible tensor runtime. Broader DFT conformer refinement has completed passing structures and energy records, while affected structures remain excluded after final force-consistency failures. New geometry dispatch is paused pending a separately gated repair pilot. Confirmed foreign-residency interruptions are preserved and retried only in new source-bound directories.

Validation: 1,158 regression tests pass; the one figure-environment skip passes separately. Both manuscript documents have resolved references and zero overfull boxes, with separate hash-bound visual reviews. The 74-page slide deck was inspected. Local entry-point checks cover 613 targets, excluding two explicitly parameterized raw-path templates.

This checkpoint does not establish a converged polymer ensemble, neighboring-site free energy, electronic coupling, reaction yield or full-axis material eligibility. Public raw-data deposition and final manuscript metadata remain unresolved. The broad repository reorganization remains deferred by the requested scientific work order.
