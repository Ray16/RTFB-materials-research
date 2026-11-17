"""Bounded execution of the installed analytic PCM curvature equations.

No surface points, electronic settings or derivative terms are removed. Dense
diagonal tensors are contracted by atom; off-diagonal tensors use the unchanged
installed CPU primitive on unions of requested surface blocks. Full-Hessian
replays and target finite-gradient checks gate production use.
"""

from collections import OrderedDict
import math
import numpy as np


def diagonal_contraction(surface, dF, d2F, vector, xp=np):
    """Exactly get_d2Sii(surface,dF,d2F) @ vector without its dense tensor."""
    derivative = dF.transpose(2, 0, 1)
    atoms, components, points = derivative.shape
    if components != 3 or d2F.shape != (atoms, atoms, 3, 3, points):
        raise ValueError("Unexpected PCM diagonal derivative layout")
    if vector.shape != (points,):
        raise ValueError("Unexpected PCM contraction vector")
    switch = surface["switch_fun"]
    charge = surface["charge_exp"]
    result = xp.empty((atoms, atoms, 3, 3), dtype=d2F.dtype)
    for atom in range(atoms):
        block = 2 * derivative[atom][None, :, None, :] * derivative[:, None, :, :]
        block /= switch**3
        block -= d2F[atom] / switch**2
        block *= math.sqrt(2 / math.pi) * charge
        result[atom] = block @ vector
    return result


class DiagonalOperator:
    def __init__(self, surface, first, second, xp, metrics):
        self.surface, self.first, self.second = surface, first, second
        self.xp, self.metrics = xp, metrics

    def __matmul__(self, vector):
        self.metrics["diagonal_contractions"] += 1
        return diagonal_contraction(
            self.surface, self.first, self.second, vector, self.xp
        )


class CurvatureBlocks:
    """Cache CPU primitive cross-blocks, never a full GPU curvature tensor."""

    def __init__(
        self, surface, primitive, transfer, metrics, cache_bytes=128 * 1024**3
    ):
        self.surface, self.primitive, self.transfer = surface, primitive, transfer
        self.metrics, self.limit = metrics, cache_bytes
        self.points = len(surface["charge_exp"])
        self.cache, self.bytes = OrderedDict(), 0

    def block(self, rows, columns):
        key = (tuple(rows), tuple(columns))
        if key not in self.cache:
            combined = np.concatenate([rows, columns])
            indices, inverse = np.unique(combined, return_inverse=True)
            local = {
                name: self.surface[name][indices]
                for name in ("charge_exp", "grid_coords", "norm_vec")
            }
            D, S = self.primitive(local, with_S=True, with_D=True)
            left, right = inverse[: len(rows)], inverse[len(rows) :]
            entries = [(key, left, right)]
            reverse = (key[1], key[0])
            if reverse != key:
                entries.append((reverse, right, left))
            for target, i, j in entries:
                values = tuple(
                    np.ascontiguousarray(np.take(np.take(array, i, axis=2), j, axis=3))
                    for array in (D, S)
                )
                size = sum(array.nbytes for array in values)
                if size > min(self.limit, 512 * 1024**2):
                    raise MemoryError(
                        "Requested PCM block exceeds bounded transfer/cache budget"
                    )
                if target != key and self.bytes + size > self.limit:
                    continue
                if target in self.cache:
                    old = self.cache.pop(target)
                    self.bytes -= sum(array.nbytes for array in old)
                while self.cache and self.bytes + size > self.limit:
                    _, old = self.cache.popitem(last=False)
                    self.bytes -= sum(array.nbytes for array in old)
                self.cache[target] = values
                self.bytes += size
                self.metrics["peak_cached_host_bytes"] = max(
                    self.metrics["peak_cached_host_bytes"], self.bytes
                )
            self.metrics["primitive_calls"] += 1
        self.cache.move_to_end(key)
        return self.cache[key]


class OffDiagonalTensor:
    def __init__(self, blocks, kind, transposed=False):
        self.blocks, self.kind, self.transposed = blocks, kind, transposed

    def transpose(self, *axes):
        if axes != (0, 1, 3, 2):
            raise ValueError("Unsupported PCM tensor transpose")
        return OffDiagonalTensor(self.blocks, self.kind, not self.transposed)

    def __getitem__(self, key):
        if len(key) != 4:
            raise ValueError("Unexpected PCM curvature slice")
        rows = np.arange(self.blocks.points)[key[2]]
        columns = np.arange(self.blocks.points)[key[3]]
        if np.ndim(rows) != 1 or np.ndim(columns) != 1:
            raise ValueError("Use surface blocks, not scalar PCM indexing")
        original_rows, original_columns = (
            (columns, rows) if self.transposed else (rows, columns)
        )
        result = self.blocks.block(original_rows, original_columns)[self.kind]
        if self.transposed:
            result = result.transpose(0, 1, 3, 2)
        result = np.ascontiguousarray(result[key[0], key[1]])
        self.blocks.metrics["uploaded_blocks"] += 1
        self.blocks.metrics["maximum_uploaded_block_bytes"] = max(
            self.blocks.metrics["maximum_uploaded_block_bytes"], result.nbytes
        )
        return self.blocks.transfer(result)
