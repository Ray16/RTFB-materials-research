# Cached raw-state reproduction

[summary.json](summary.json) records a relocated standard-library-only run with Python isolated mode and site initialization disabled. [reconstruction.json](reconstruction.json) contains the potential and both grafting reaction sums reconstructed from ten cached raw states. The [example](../../../examples/cached_thermochemistry/README.md) includes 58 hashed raw artifacts and the original provenance. This closes the small cached arithmetic example, not a full fresh quantum environment or public archival release.

The complete suite passes 782 tests, with one missing plotting-library dependency in the scientific environment. That test passes separately in the figure environment. The earlier full run invoked pytest through stdin, which prevented two multiprocessing spawn tests from locating __main__; its failed XML is preserved. Rerunning through the normal `python -m pytest` entry point resolves both without a code change.

A separate export from committed Git revision `836a01d` also reproduces exactly after extraction away from the checkout; [fresh_git_archive.json](fresh_git_archive.json) records the archive hash and isolated command.
