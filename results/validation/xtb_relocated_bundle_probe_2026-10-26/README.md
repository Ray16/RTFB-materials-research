# Relocated reviewer-bundle xTB kernel probe

[Evidence and hashes](summary.json) · [Probe output](probe.json) ·
[Tests](tests.xml) · [Reviewer bundle](../../releases/manuscript_bundle_v2/README.md)

The portable xTB kernel entry point from the extracted version 2 reviewer bundle
ran using an external, hash-qualified xTB prefix. The atomic masses and five
thermodynamic values matched the previously qualified reference exactly. The
helper's process map confirmed the staged library path. A Python audit hook
blocked opens under the original project, while an `strace -f` open-system-call
trace covered the parent and child processes; the final trace contained no
absolute open of the original project path.

An initial run reproduced the numbers but Python startup read an editable-install
`.pth` path into the original checkout before the audit hook was installed. That
run is retained as a failed isolation attempt. The final run used `python -S`,
explicitly added the existing environment's site-packages directory, and placed
the relocated source first on the import path.

The external prefix is a local test artifact, assembled from the existing qualified
installation using hardlinks where possible, with system-library dependencies. It
is not a distributable lockfile. This is a kernel probe using an existing Python
environment: no quantum calculation was run, and native/projected production
entrypoints, a fresh platform installation, GPU portability and conformer sampling
are not qualified by this result. The raw trace, staged prefix and relocated
checkout are retained on shared storage at the paths in [summary.json](summary.json).
