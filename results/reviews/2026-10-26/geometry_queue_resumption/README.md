# Occupied-conformer geometry queue recovery

The [preflight audit](preflight_audit.json) checks the previous dispatcher process
identity, availability of its dispatch lock, every declared source hash and the
current per-case geometry statuses. Nine batches have a stopped dispatcher and
unfinished geometry; the tenth has all nineteen expansion geometries checked.
These geometry checks do not qualify hybrid electronic energies, native free
energies or ensemble convergence.

The new source-bound geometry adapter preserves the scientific producer and its
acceptance criteria. It isolates per-host reservation inventory failures, shares
short-lived discovery records across controllers, and repeats live reservation
checks immediately before placement. CPU-only node preflight defers GPU imports
until the reserved worker. Unknown remote ownership and failed or partially
written scientific attempts remain preserved rather than being reset. Tests cover
cache expiry, policy changes, unreachable hosts and adaptation of the existing
queue contract.

Existing per-case intent, launch and scientific records remain authoritative;
these checks alone do not assert that a worker has started or completed.

The [launch verification](launch_verification.json) records nine live dispatchers and a separate live independent reviewer. At that snapshot the nine batches had 58 checked geometries, 101 pending inputs and 12 unresolved ownership cases. No GPU slot was available in the verified inventory; no new scientific start was established. Live batch status is in each audited batch directory; controller and observer launch records and logs are under `calcs/validation/2026-10-26_geometry_queue_recovery/`.

A subsequent [ownership reconciliation](../geometry_intent_reconciliation/README.md) verified all twelve unresolved intents as unstarted using process checks across every permitted host, archived their original records and reopened them. The launch snapshot above remains unchanged as a historical record.

The [latest recovery verification](../geometry_timeout_reconciliation/verification.json) includes three later remote timeouts and their separately audited recovery: nine live controllers, 77 independently checked geometries and 113 pending inputs, with disjoint CPU allocations and no unresolved ownership at the recorded time.
