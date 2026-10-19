# Fresh-gas public API forwarding

After actual owner/worker termination was verified for the parent spin campaigns,
`properties/nelsen.py` was updated to use the qualified shared gas solver. Its
Atoms/eV interface, functional/basis arguments and force-field embedding remain
compatible. Five tests cover finite/converged energy admission, eV conversion,
and refusal to rescore a stopped geometry. Prior producer bytes are preserved.
No saved reference output is overwritten or retrospectively certified.
