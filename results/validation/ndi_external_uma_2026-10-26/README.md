# Superseded NDI UMA adapter draft

This unstarted draft is preserved. Its adapter could replace a live process via
`exec`, which would discard a process-local GPU residency guard. The
[corrected campaign](../ndi_external_uma_2026-10-26_v2/README.md) instead requires
the scheduler to select the matching interpreter before reservation. The exact
draft source is retained at `calcs/validation/2026-10-26_ndi_external_control_seeds_v2/uma_adapter_v1.py`.
No scientific worker started from this draft.
