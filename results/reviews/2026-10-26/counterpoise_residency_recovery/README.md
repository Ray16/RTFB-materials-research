# Adjacent-contact GPU residency recovery

[Stop review](review.json) · [Identical geometry review](prepared_geometry.json) · [Summary](summary.json) · [Tests](tests.xml)

The first attempt was interrupted when a foreign process appeared on the selected GPU before the scientific calculation began. The residency guard signaled only its own worker. The nested signal handler reported a host-RAM message, but the saved residency record identifies the actual interruption and the recorded RAM samples did not breach the configured limits. The original attempt and logs remain unchanged.

The original coordinator, gate owner and worker were verified stopped, and no scientific output or raw calculation directory was created. A new separate case has the same geometry, task and numerical settings and has passed the independent 38-atom/ghost-basis construction review. The [retry coordinator](retry_launch_verification.json) is live and waiting for a suitable uncontended device; no retry scientific calculation has started at the recorded snapshot. No counterpoise interaction or MeCN correction is yet established.

The [prelaunch recheck](prelaunch_recheck.json) again verified that the original owners were stopped and no scientific output had been created before the new coordinator was started. Live retry status is stored with the new case.
