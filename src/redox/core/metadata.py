"""CPU-only validation of diagnostic metadata before expensive qualifiers start."""
import json
from pathlib import Path
import numpy as np


def plain_metadata(value):
    if isinstance(value,np.generic):return plain_metadata(value.item())
    if isinstance(value,Path):return str(value)
    if isinstance(value,np.ndarray):
        if value.size>4096:raise ValueError('Store large scientific arrays separately from diagnostic JSON')
        return plain_metadata(value.tolist())
    if isinstance(value,dict):return {str(k):plain_metadata(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [plain_metadata(v) for v in value]
    if value is None or isinstance(value,(str,int,float,bool)):return value
    raise TypeError(f'Unsupported diagnostic metadata type: {type(value).__name__}')


def validate_metadata(value):
    result=plain_metadata(value)
    json.dumps(result,allow_nan=False)
    return result


def planned_allocation_metadata_check():
    example=dict(shape=[np.int64(632),np.int32(632),np.uint64(714)],bytes=np.int64(2281509888),
                 device_free_bytes=np.uint64(30519001088),pool_used_bytes=np.int64(167668224),
                 relative_error=np.float64(1e-10),elapsed_s=np.float32(.125),native_gates_passed=np.bool_(True))
    return validate_metadata(example)

if __name__=='__main__':
    print('DIAGNOSTIC_METADATA_PREFLIGHT_OK',json.dumps(planned_allocation_metadata_check(),sort_keys=True))
