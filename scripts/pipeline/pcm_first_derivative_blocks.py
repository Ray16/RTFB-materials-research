"""Lazy PCM first derivatives using the unchanged installed pairwise primitive.

Only requested rectangles are assembled. Each tile is obtained by evaluating
the original primitive on the union of its row/column points, preserving point
identity, diagonal terms, normal vectors, precision and transpose signs.
"""
import numpy as np


class FirstDerivativeBlocks:
    def __init__(self, surface, primitive, xp, metrics, tile_points=512,
                 maximum_output_bytes=512*1024**2):
        if type(tile_points) is not int or tile_points <= 0:
            raise ValueError("Positive integer PCM tile size required")
        self.surface, self.primitive, self.xp, self.metrics = surface, primitive, xp, metrics
        self.points = len(surface["charge_exp"])
        self.tile = tile_points
        self.maximum_output_bytes = maximum_output_bytes

    def rectangle(self, rows, columns, kind):
        shape = (3, len(rows), len(columns))
        size = int(np.prod(shape))*8
        if size > self.maximum_output_bytes:
            raise MemoryError("Requested PCM first-derivative rectangle exceeds the explicit bound")
        self.metrics["maximum_output_bytes"] = max(self.metrics["maximum_output_bytes"], size)
        output = self.xp.empty(shape, dtype=self.xp.float64)
        for i in range(0, len(rows), self.tile):
            for j in range(0, len(columns), self.tile):
                r, c = rows[i:i+self.tile], columns[j:j+self.tile]
                indices, inverse = np.unique(np.concatenate((r, c)), return_inverse=True)
                local = {name: self.surface[name][indices]
                         for name in ("charge_exp", "grid_coords", "norm_vec")}
                D, S = self.primitive(local, with_D=(kind == 0), with_S=True)
                tensor = D if kind == 0 else S
                left, right = inverse[:len(r)], inverse[len(r):]
                block = self.xp.take(self.xp.take(tensor, left, axis=1), right, axis=2)
                output[:, i:i+len(r), j:j+len(c)] = block
                self.metrics["primitive_calls"] += 1
                self.metrics["maximum_primitive_points"] = max(
                    self.metrics["maximum_primitive_points"], len(indices))
                self.metrics["maximum_primitive_tensor_bytes"] = max(
                    self.metrics["maximum_primitive_tensor_bytes"],
                    3*len(indices)**2*8*(2 if kind == 0 else 1))
                del D, S, tensor, local, block
        return output


class FirstDerivativeTensor:
    def __init__(self, blocks, kind, transposed=False, sign=1):
        self.blocks, self.kind, self.transposed, self.sign = blocks, kind, transposed, sign
        self.shape = (3, blocks.points, blocks.points)

    def transpose(self, *axes):
        if axes != (0, 2, 1):
            raise ValueError("Unsupported PCM first-derivative transpose")
        return FirstDerivativeTensor(self.blocks, self.kind, not self.transposed, self.sign)

    def __neg__(self):
        return FirstDerivativeTensor(self.blocks, self.kind, self.transposed, -self.sign)

    def __getitem__(self, key):
        if not isinstance(key, tuple) or len(key) != 3:
            raise ValueError("Expected a three-axis PCM first-derivative slice")
        rows = np.arange(self.blocks.points)[key[1]]
        columns = np.arange(self.blocks.points)[key[2]]
        if np.ndim(rows) != 1 or np.ndim(columns) != 1:
            raise ValueError("Scalar PCM surface indexing is outside this execution profile")
        r, c = (columns, rows) if self.transposed else (rows, columns)
        result = self.blocks.rectangle(r, c, self.kind)
        if self.transposed:
            result = result.transpose(0, 2, 1)
        if self.sign != 1:
            result *= self.sign
        self.blocks.metrics["rectangles"] += 1
        return result[key[0]]
