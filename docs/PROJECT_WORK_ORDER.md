# Project work order

User instruction, 21 October 2026: complete and review the four scientific tasks
before the repository refactor. In particular, properly address the oligomer work
before reorganizing files. Current result/status entry points are in
[results/README.md](../results/README.md).

The four tasks are property validation, systematic candidate expansion, defined
grafting thermodynamics, and neighboring-motif/polymer modeling. Closing a task
requires reporting qualified findings, explicit exclusions and remaining physical
limits; prepared code or an optimizer checkpoint alone is not a scientific result.

Current priority is the fourth task: continue the checked oligomer geometries,
complete the matched spacing panel, distinguish electronic nonadditivity from
qualified conditional grafting G, and preserve geometry/provenance checks.

After that scientific review, perform the requested repository refactor:

- Put maintained figures in one consistent location with stable entry points.
- Keep plotting scripts in a separate, organized plotting directory.
- Inventory redundant, one-off and obsolete files and move them to an archive.
- Preserve immutable evidence, source snapshots, raw calculation locations and
  hashes. Keep historical provenance resolvable rather than silently rewriting it.
- Update imports, commands, links and result indexes; verify regeneration and
  relevant tests before committing and pushing the refactor.

No repository-wide moves or archival cleanup are authorized to start before the
scientific work above is properly addressed. Active calculations retain their
frozen numerical producers and inputs throughout the current continuation.
