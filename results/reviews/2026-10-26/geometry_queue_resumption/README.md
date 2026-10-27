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
