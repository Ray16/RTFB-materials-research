# Occupied electronic gate-contention recovery

[Startup snapshot](summary.json) · [Transition audit](transition_audit.json) ·
[Tests](tests_final.xml) · [Earlier working-directory recovery](../electronic_workdir_recovery/README.md)

Two occupied electronic placements reached the atomic GPU gate after the GPUs
had become busy. Their preserved logs state that the gate did not run the worker;
neither case had a scientific stage, worker guard or raw calculation output.
The old controller treated the exited gate process as an unreviewed stopped
attempt, leaving those inputs held.

The revised controller reopens an entry only when the recorded owner is verified
dead, the gate log exactly reports a reservation refusal, and no scientific
output or guard exists. It archives the original launch, intent and log before
marking the attempt unstarted. Unknown, failed or partially started cases remain
held. Before transitioning the owned CPU controllers, all nine permitted hosts
were checked for matching workers and the old gate owners were checked; no
scientific worker was stopped.

At the [recorded snapshot](summary.json), both proven gate refusals had been
archived and returned to pending, and all three replacement CPU controllers were
live. This is not a completed electronic calculation. Fourteen focused tests
cover explicit contention, ambiguous/started cases, controller-template binding,
inventory and resource preflight. The live batch dispatch files can change after
the snapshot.
