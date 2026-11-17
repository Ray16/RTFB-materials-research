"""Historical import path for the certified sampling snapshot."""
import sys

from redox.quality import sampling_receipts as _implementation

sys.modules[__name__] = _implementation
