# Remote working-directory recovery

[Checked snapshot](summary.json) · [Retirement audit](retirement_audit.json) ·
[Old-adapter dry run](remote_dryrun_v1.log) · [Corrected-adapter dry run](remote_dryrun_v2.log) · [Tests](tests.xml)

Two electronic placement responses became ownership-unverified when a GPU
briefly appeared free. Neither case had a launch, attempt, guard or raw output.
The previous controller did not persist the SSH exception, so its exact cause
cannot be established retrospectively.

A separate dry run on the same host and case reproduced a concrete defect:
SSH started the old remote adapter in the account's home directory, where a
repository-relative source path could not be found. The corrected adapter first
changes to the repository root. With a deliberately invalid GPU index, the
corrected dry run passed source verification and reached GPU-index validation;
neither dry run reserved a GPU or started a scientific worker. The original
exceptions are consistent with this defect, but their cause remains inferred.

After inspecting all nine currently permitted hosts and confirming that no
worker referenced the two cases, the three idle, owned CPU controllers were
retired. Their two unstarted intents were archived and reopened, and three
replacement coordinators started. Placement exceptions are now written to a
case-local record while uncertain ownership remains held. At the recorded
snapshot the replacements were live, both cases were pending, and no new
scientific worker had started. The native coordinator was not replaced.

The scientific worker, geometry, model, numerical settings, reservation gate
and physical QC were not changed. Fourteen focused adapter, inventory and
recovery tests passed. This is a placement fix, not evidence of a completed
polymer calculation; current batch dispatch files give live status.
