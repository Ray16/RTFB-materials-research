"""Recognize the exact runtime-location refactor of the qualified qRRHO source.

Only one literal source transformation is admitted. Restoring its historical
library assignment must reproduce the original certificate's full source hash.
Numerical code, constants, imports elsewhere and even other formatting changes
remain source-bound. The caller must still verify the actual library, worker,
protocol and all certificate evidence independently.
"""
import hashlib

_QUALIFIED_LOCATION = b"import numpy as np\n\nLIBRARY = Path('/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/lib/libxtb.so.6.7.1')"
_DISCOVERED_LOCATION = b"import numpy as np\nfrom redox.qm.xtb_runtime import configured_library\n\nLIBRARY = configured_library()"


def matches_qualified_source(source: bytes, expected_sha256: str) -> bool:
    """Check exact source identity or the single reviewed location-only change."""
    if hashlib.sha256(source).hexdigest() == expected_sha256:
        return True
    if source.count(_DISCOVERED_LOCATION) != 1:
        return False
    restored = source.replace(_DISCOVERED_LOCATION, _QUALIFIED_LOCATION, 1)
    return hashlib.sha256(restored).hexdigest() == expected_sha256
