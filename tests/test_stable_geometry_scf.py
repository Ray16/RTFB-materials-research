"""SCF recovery must use the scanner's geometry and preserve requested tolerances."""

import copy
import json
import numpy as np
import pytest
from redox.qm.stable_geometry_scf import stabilize


class Molecule:
    def __init__(self, value):
        self.value = value

    def atom_coords(self):
        return np.array([[self.value, 0.0, 0.0]])


class MeanField:
    def __init__(self, outcomes):
        self.mol = Molecule(1.0)
        self.outcomes = iter(outcomes)
        self.seen = []
        self.max_cycle = 200
        self.level_shift = self.damp = 0.0
        self.conv_tol, self.conv_tol_grad = 1e-12, 1e-8
        self.converged, self.e_tot = False, None

    def kernel(self, dm0=None):
        self.seen.append((self.mol.value, self.level_shift, self.damp))
        self.converged = next(self.outcomes)
        self.e_tot = -10.0
        return self.e_tot

    def make_rdm1(self):
        return np.eye(1)


def test_scanner_copy_uses_current_geometry_and_recovers_without_looser_tolerances(
    tmp_path,
):
    original = stabilize(MeanField([False, True]), tmp_path / "events.jsonl")
    scanner = copy.copy(original)
    scanner.mol = Molecule(2.0)
    scanner.kernel()
    assert scanner.converged
    assert scanner.seen == [(2.0, 0.0, 0.0), (2.0, 0.3, 0.25)]
    assert scanner.conv_tol == 1e-12 and scanner.conv_tol_grad == 1e-8
    assert scanner.level_shift == scanner.damp == 0.0
    events = [
        json.loads(line)
        for line in (tmp_path / "events.jsonl").read_text().splitlines()
    ]
    assert [row["stage"] for row in events] == ["ordinary", "level_shift_damping"]


@pytest.mark.parametrize("converged", [True, False])
def test_second_order_state_is_copied_and_nonconvergence_stays_explicit(
    tmp_path, converged
):
    mean_field = MeanField([False, False])

    class Newton:
        def kernel(self, dm0=None):
            self.converged = converged
            self.e_tot = -11.0
            self.mo_coeff = self.mo_occ = self.mo_energy = np.ones(1)
            return self.e_tot

    mean_field.newton = lambda: Newton()
    recovered = stabilize(mean_field, tmp_path / "events.jsonl")
    recovered.kernel()
    assert recovered.converged is converged
    assert recovered.e_tot == -11.0
    assert recovered.conv_tol == 1e-12 and recovered.conv_tol_grad == 1e-8
