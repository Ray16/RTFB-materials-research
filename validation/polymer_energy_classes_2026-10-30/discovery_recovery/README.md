# Availability discovery and contention recovery

[Summary](summary.json) · [Stopped controller](stopped_controller.json) ·
[Resume receipt](resume_launch.json) · [Actual starts](actual_starts.json) ·
[Tests](tests.xml)

A cluster-wide GPU availability query exceeded its timeout and terminated the
CPU dispatcher. Scientific workers and the independent collector retained their
ownership and outputs. The dispatcher is now explicitly resumed after checking
the prior PID/start identity. Availability is queried independently per node
with a bounded timeout; an unverified node contributes no free GPUs while
verified peer observations remain usable.

One placement had been refused by the GPU reservation gate before any scientific
stage started. Its controller and namespace child were both verified stopped;
there was no residency guard, electronic stage or raw output. The original
claim was preserved in the contention history and released without consuming a
scientific attempt. Unreachable owners, live producers and any scientific
outputs prevent this release path.

Twenty focused checks cover the queue and source bindings, including peer-node
timeouts, malformed availability, known gate refusals and preservation of
uncertain/live scientific ownership. Four resumed energy predictors were
verified by exact host/PID/start identity. No existing scientific output or
quantum setting was replaced.
