# Final-search DFT re-ranking and native ensemble follow-up

The [request](request.json) declares twelve neutral phthalimide graph/configuration
cases and nested 10/20-start batches. A final CREST source must finish normally
and pass independent frame, graph, stereochemistry and physical-input checks
before its DFT inputs are frozen. MD frames and intermediate ensembles are excluded.

The all-same configuration retains its two previously checked DFT poses before
adding final CREST energy-ranked starts. The alternate configuration uses final
CREST starts. These references prevent a previously known low-energy backbone
pose from being lost when switching searches. DFT minima are deduplicated within
each prefix; rediscovery never adds an artificial entropy weight.

- [Preparation status](preparation_status.json): actual frozen graph/start counts.
- [Current results and matched cycles](current.json): pending G stays blank.
- [Finite controller](coordinator.json): actual owner and source hashes.
- [Node preflight](lambda13_preflight.json) and [RAM-promise check](host_RAM_promise_check.json).
- [Baseline continuation](../2026-10-23_baseline_two_pose_ensemble_v2/continuation_status.json).

New scientific calculations remain subject to real energy/native resource pilots,
the shared GPU reservation gate, disjoint two-core workers, conservative CDS
forces, total-energy derivative checks and native Hessian admission. RAM admission
counts the promised growth of already launched workers, subtracting only their
observed unique allocations. A free GPU alone does not establish host capacity.

Launch observation is asynchronous and retains PID/start ownership. A live launch
observer is polled rather than duplicated. Temporarily unverified devices remain
unavailable and do not consume molecule attempts. An ended observer without a
matching outcome requires review. The original baseline monitor stopped after a
transient device-identity query timeout; all twelve native results subsequently
completed without restarting scientific jobs. Observer history and source snapshots
remain here. The source upgrade changed no numerical/physical acceptance threshold.

The initial finite partitions, once present, still need review of minimum and
ensemble changes at 10 versus 20 starts, rotor-cutoff sensitivity and the two
backbone configurations. More sampling is required when the declared diagnostic
fails. No global-minimum, experimental sequence-population or synthesis-admission
claim follows from completing this queue.

The directory also provides a shared resource root for explicitly named candidate
baseline cases. Their own declarations/results remain with each candidate campaign;
they are not counted as completed phthalimide search cases.
