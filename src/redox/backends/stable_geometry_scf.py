"""Process-local SCF recovery without loosening geometry or electronic tolerances."""

import hashlib
import json
import time


def stabilize(mean_field, log_path):
    """Install a scanner-safe class method, so every call uses its own geometry."""
    parent = type(mean_field)
    kernel = parent.kernel
    if getattr(parent, "_redox_stable_geometry", False):
        return mean_field

    def record(method, stage, error=None):
        import numpy as np

        value = getattr(method, "e_tot", None)
        event = dict(
            stage=stage,
            converged=bool(method.converged),
            energy_Ha=float(value) if value is not None else None,
            geometry_sha256=hashlib.sha256(
                np.asarray(method.mol.atom_coords()).tobytes()
            ).hexdigest(),
            conv_tol=method.conv_tol,
            conv_tol_grad=method.conv_tol_grad,
            updated=time.time(),
            error=error,
        )
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a") as stream:
            stream.write(json.dumps(event) + "\n")

    def recovered(self, *args, **kwargs):
        shift, damping = self.level_shift, self.damp
        self.max_cycle = max(self.max_cycle, 300)
        try:
            energy = kernel(self, *args, **kwargs)
            record(self, "ordinary")
            if self.converged:
                return energy
            self.level_shift, self.damp = 0.3, 0.25
            energy = kernel(self, dm0=self.make_rdm1())
            record(self, "level_shift_damping")
            if self.converged:
                return energy
            self.level_shift, self.damp = 0.0, 0.0
            second = self.newton()
            second.conv_tol, second.conv_tol_grad = self.conv_tol, self.conv_tol_grad
            second.max_cycle = 150
            energy = second.kernel(dm0=self.make_rdm1())
            for name in ("mo_coeff", "mo_occ", "mo_energy", "e_tot", "converged"):
                setattr(self, name, getattr(second, name))
            record(self, "second_order")
            return energy
        except Exception as error:
            record(self, "exception", repr(error))
            raise
        finally:
            self.level_shift, self.damp = shift, damping

    mean_field.__class__ = type(
        "StableGeometry" + parent.__name__,
        (parent,),
        {"kernel": recovered, "_redox_stable_geometry": True},
    )
    return mean_field
