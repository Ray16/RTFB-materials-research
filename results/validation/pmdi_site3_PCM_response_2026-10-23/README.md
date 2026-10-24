# Preserved first response/resource attempt

The exact saved site-3 geometry completed its R2SCAN/SVP baseline SCF and the
blocked first-derivative assembly, but failed during construction of a dense R
operator in the fixed-grid-potential response. [status.json](status.json) retains
the failure; raw records remain under `calcs/validation/pmdi_site3_PCM_response_2026-10-23/`.
The derivative stage used bounded rectangles, but the remaining dense operator
still exceeded device memory. No native Hessian or G was admitted.

The [fresh successor](../pmdi_site3_PCM_response_v2_2026-10-23/README.md) uses
separately controlled contracted R applications. The original output is neither
overwritten nor treated as a completed response benchmark.
