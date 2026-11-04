# Upstream hybrid-energy continuations for missing backbone conformers

[Summary](summary.json) · [Input and launch bindings](inputs_and_launches.json) ·
[Actual worker inspection](actual_workers.json) · [Tests](tests.xml)

Two hybrid electronic-energy calculations are active for the missing alternating
empty-backbone inputs 10 and 19. Each runs on its own V100-32GB GPU, with disjoint
two-CPU allocations, two-thread pools and the original 128 GiB prospective host
RAM ceiling plus 32 GiB free-memory floor. The first has passed the prior 80 GiB
RSS ceiling while remaining below its predeclared 128 GiB bound. Resource bounds
do not themselves establish successful electronic or thermochemical results.

Both plans were prepared before the refactor but had no scientific stage or raw
output. A nine-host process inventory found no surviving predictor for the four
missing empty-backbone inputs. Exact plans, scored geometries and original
scientific source bytes are retained unchanged. Each node-local namespace probe
checks 1,507 combined source/input bindings and the declared runtime distributions
before submission. The actual worker checks its cuTENSOR FP64 runtime under the
shared GPU gate. Live RAM promises from the twelve native workers constrain
placement; memory is not pooled across GPUs.

These previously unstarted raw-batch directories had no archived mountpoint. The
energy namespace therefore uses a private temporary parent and rebinds every
existing archived raw-batch sibling read-only before exposing the new writable
raw directory. Nothing is created inside the source archive. Input/source hashes
are checked again through the resulting original-path view.

The first input-19 attempt stopped because a missing CuPy kernel tried to write
its cache under the read-only home directory. That failed attempt and its closed
guard remain preserved in [the failure binding](preserved_cache_failure.json).
The fresh retry redirects only CuPy's compilation cache to private namespace
temporary storage; the method, geometry, precision and physical gates are
unchanged. Input 10 retains its actual original cache environment.

The independent observer runs the unchanged archived electronic reviewer only
after a successful stage and closed ownership. Its live report is
`results/validation/polymer_energy_namespace_2026-10-30/independent_reviews/status.json`.
Sixteen focused software tests pass. Neither electronic energy nor native free
energy is admitted by this startup milestone. Native Hessians and thermochemistry
for these exact geometries require separate follow-up after the electronic
review; the two other missing empty-backbone inputs and the wider occupancy
campaign remain open.

The versioned entry points are
`python -m redox.workflows.archived_energy_namespace prepare|probe|launch` and
`scripts/pipeline/observe_archived_energy_namespace.py`. Read the fleet skill,
inspect exact previous ownership, and use a fresh destination before submission.
The source archive is still an explicit dependency; a public standalone release
remains a separate requirement.
