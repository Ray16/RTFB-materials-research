# Stopped-search conformer warm starts

[Preparation/QC](preparation_status.json) · [Declared inputs](request.json) · [Versioned builder](../../../../scripts/pipeline/stage_interrupted_conformers.py)

Eight original low-level searches ended at their 48-hour limits. Six have a recorded completed metadynamics phase and now supply 20 checked DFT starting structures each, including inherited reference poses where specified. Two remain unprepared under the declared phase requirement. Together with the four completed original searches, this makes DFT inputs available for ten of the twelve graph/configuration cases.

These are **warm starts from incomplete searches**, not completed CREST ensembles, DFT minima or population weights. Every selected frame is checked against the exact open-chain, occupancy and stereochemical graph. The source ensemble/log/terminal-owner hashes and checked source indices are retained; original failures and raw files remain unchanged. No quantum jobs are launched by this preparation step. Full individual geometry/energy checks, native thermochemistry and convergence of the 10/20-start batches remain required.
