# Consolidated manuscript reviewer bundle

The [current record](current.json) identifies a local archive of the current manuscript source and evidence plus five overlapping raw-state packages, including the TEMPO addendum. Its exact source commit, archive hash and shared-storage locations are recorded there. The packages overlap, so their state counts must not be added.

The [relocated standard-library replay](stdlib_replay.json) checks all six arithmetic stages with Python opens of the original checkout blocked. Independent native replays reconstructed [242 primary Hessians](primary_native_isolation.json), [13 neighboring-cycle Hessians](neighbor_native_isolation.json), [8 numerical-control Hessians](controls_native_isolation.json), and [4 TEMPO/Fc Hessians](tempo_native_isolation.json). The primary replay used a [fresh hash-locked CPU environment](primary_replay_environment.json); its lock is [here](../../../environments/reviewer_primary_cpu.lock.txt).

These checks reproduce saved numerical records; they do not rerun SCF calculations, repeat stationarity or spin checks, demonstrate experimental accuracy or ensemble convergence, or qualify a production GPU environment. The archive is local, not a public deposit. Data licensing, author metadata, a DOI and wider deployment portability remain open.
