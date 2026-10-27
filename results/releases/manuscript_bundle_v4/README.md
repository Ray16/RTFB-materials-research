# Consolidated manuscript reviewer bundle

The [current record](current.json) identifies a local archive of the current manuscript source and evidence plus six overlapping raw-state packages, including the TEMPO and local-motif addenda. Its source commit, archive hash and shared-storage locations are recorded there. The packages overlap, so state counts must not be added.

The [relocated replay](stdlib_replay.json) checks seven stages with Python opens of the original checkout blocked. It includes [12 local-motif interaction triplets](local_motif_replay.json) from 36 reviewed electronic geometries. Independent native replays reconstruct [242 primary Hessians](primary_native_isolation.json), [13 neighboring-cycle Hessians](neighbor_native_isolation.json), [8 numerical-control Hessians](controls_native_isolation.json), and [4 TEMPO/Fc Hessians](tempo_native_isolation.json). The primary replay uses a [hash-locked CPU environment](../../../environments/reviewer_primary_cpu.lock.txt).

These checks reproduce saved records; they do not rerun SCF calculations, repeat physical QC, demonstrate experimental accuracy or ensemble convergence, or qualify a production GPU environment. The archive is local, not a public deposit. Data licensing, author metadata, a DOI and wider deployment portability remain open.
