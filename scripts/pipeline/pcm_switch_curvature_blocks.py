"""Host-backed switch curvature from unchanged all-atom grid tiles.

All atoms remain in every primitive evaluation. Only the surface-point axis is
partitioned; original per-atom ownership, switching weights and area are retained.
"""
import numpy as np


def local_surface(surface, begin, end):
    local = dict(surface)
    for key in ('grid_coords','switch_fun','area','charge_exp','norm_vec'):
        if key in surface:
            local[key] = surface[key][begin:end]
    size = end-begin
    local['gslice_by_atom'] = [(min(max(a-begin,0),size),min(max(b-begin,0),size))
                              for a,b in surface['gslice_by_atom']]
    return local


class HostSwitchTensor:
    def __init__(self, grid_major, xp, metrics, maximum_upload_bytes=512*1024**2):
        self.storage = grid_major
        self.array = grid_major.transpose(1,2,3,4,0)
        self.shape, self.ndim, self.dtype = self.array.shape, self.array.ndim, self.array.dtype
        self.xp, self.metrics, self.limit = xp, metrics, maximum_upload_bytes

    def upload(self, value):
        if value.nbytes > self.limit:
            raise MemoryError('Switch-curvature upload exceeds its bounded block size')
        self.metrics['maximum_uploaded_bytes'] = max(self.metrics['maximum_uploaded_bytes'],value.nbytes)
        self.metrics['uploads'] += 1
        return self.xp.asarray(np.ascontiguousarray(value))

    def __getitem__(self, key):
        return self.upload(self.array[key])

    def flatten(self):
        # The original dense diagonal primitive is useful for small controls.
        # Large native jobs must use the separately qualified diagonal operator.
        return self.upload(self.array).reshape(-1)

    def __matmul__(self, vector):
        if vector.shape != (self.shape[-1],):
            raise ValueError('Expected a surface-vector curvature contraction')
        result = self.xp.empty(self.shape[:-1],dtype=self.dtype)
        for atom in range(self.shape[0]):
            result[atom] = self[atom] @ vector
        return result


def assemble(surface, primitive, xp, metrics, tile_points=32, maximum_host_bytes=None):
    n, atoms = len(surface['grid_coords']),len(surface['atom_coords'])
    required = 2*n*atoms*atoms*9*8
    if maximum_host_bytes is not None and required > maximum_host_bytes:
        raise MemoryError('Host switch-curvature tensors exceed the declared bound')
    storage = [np.empty((n,atoms,atoms,3,3),dtype=np.float64) for _ in range(2)]
    metrics['host_tensor_bytes'] = max(metrics['host_tensor_bytes'],required)
    for begin in range(0,n,tile_points):
        end = min(n,begin+tile_points)
        values = primitive(local_surface(surface,begin,end))
        for target,value in zip(storage,values):
            # The installed primitive returns atom/atom/xyz/xyz/grid order.
            current = value.transpose(4,0,1,2,3)
            host = current.get() if hasattr(current,'get') else np.asarray(current)
            if not np.isfinite(host).all():
                raise ValueError('Nonfinite analytic switch-curvature tile')
            target[begin:end] = host
        metrics['primitive_calls'] += 1
        metrics['maximum_tile_tensor_bytes'] = max(metrics['maximum_tile_tensor_bytes'],2*(end-begin)*atoms*atoms*9*8)
        del values,current,host
    for array in storage:
        array.flags.writeable = False
    return tuple(HostSwitchTensor(array,xp,metrics) for array in storage)
