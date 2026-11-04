"""Finite conformer partitions, after independent state and identity checks.

These functions do not certify input quantum states. The caller must supply
checked free energies of one chemical species, backbone configuration and
protocol. Search hit counts are never statistical degeneracies. Each unique
minimum contributes once; electronic spin and rotational symmetry are already
included in each state's nuclear/electronic partition.
"""
from __future__ import annotations

import math

import numpy as np

KB_EV_K = 1.380649e-23 / 1.602176634e-19
EV_KJ_MOL = 1.602176634e-19 * 6.02214076e23 / 1000


def finite_number(value, name):
    if type(value) not in (float, int) or not math.isfinite(value):
        raise ValueError(f"Finite {name} required")
    return float(value)


def partition(free_energies_eV, temperature_K=298.15):
    """A stable log-sum over unique, equally counted conformational wells.

    Populations are conditional on this finite set and its thermodynamic model,
    not measured bead populations or the frequency of CREST discoveries.
    """
    temperature_K = finite_number(temperature_K, "temperature")
    if temperature_K <= 0 or not free_energies_eV:
        raise ValueError("A positive temperature and nonempty minimum set required")
    values = [finite_number(value, "conformer G") for value in free_energies_eV]
    minimum = min(values)
    kt = KB_EV_K * temperature_K
    weights = [math.exp(-(value - minimum) / kt) for value in values]
    total = math.fsum(weights)
    return dict(
        G_eV=minimum - kt * math.log(total),
        minimum_G_eV=minimum,
        configurational_lowering_kJ_mol=kt * math.log(total) * EV_KJ_MOL,
        conditional_populations=[value / total for value in weights],
        temperature_K=temperature_K,
        unique_minima=len(values),
    )


def cluster_minima(records, rmsd_A, rmsd_cutoff_A=0.1, energy_tolerance_kJ_mol=0.2):
    """Deterministic complete-link grouping of nearly identical DFT minima.

    Every pair within a group must meet both thresholds. This avoids transitive
    chains merging distinct endpoints. Geometry and energy criteria are numerical
    duplicate heuristics; sensitivity to their values must be reported separately.
    The input matrix must use a proper rotation and symmetry-aware heavy atoms.
    """
    cutoff = finite_number(rmsd_cutoff_A, "RMSD cutoff")
    tolerance = finite_number(energy_tolerance_kJ_mol, "energy tolerance")
    if cutoff <= 0 or tolerance < 0 or not records:
        raise ValueError("Positive RMSD cutoff and nonnegative energy tolerance required")
    keys = [row["key"] for row in records]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate input conformer identity")
    energies = [finite_number(row["E_eV"], "electronic energy") for row in records]
    distances = np.asarray(rmsd_A, dtype=float)
    n = len(records)
    if (
        distances.shape != (n, n)
        or not np.isfinite(distances).all()
        or np.any(distances < 0)
        or not np.allclose(distances, distances.T, atol=1e-10, rtol=0)
        or not np.allclose(np.diag(distances), 0, atol=1e-10, rtol=0)
    ):
        raise ValueError("Finite symmetric nonnegative RMSD matrix required")
    groups = []
    for i in sorted(range(n), key=lambda index: (energies[index], keys[index])):
        for group in groups:
            if all(
                distances[i, j] <= cutoff
                and abs(energies[i] - energies[j]) * EV_KJ_MOL <= tolerance
                for j in group
            ):
                group.append(i)
                break
        else:
            groups.append([i])
    return [[keys[i] for i in group] for group in groups]


def summarize(records, groups, expected_keys, temperature_K=298.15):
    """Partition only a complete declared set, with explicit duplicate groups.

    The minimum-G member represents a duplicate group; adding another search hit
    of the same minimum cannot create an extra entropy term. Incomplete sets
    retain their missing keys and never substitute electronic energies for G.
    """
    expected = list(expected_keys)
    lookup = {row["key"]: row for row in records}
    if len(lookup) != len(records) or len(set(expected)) != len(expected):
        raise ValueError("Repeated conformer identity")
    if not expected or not set(lookup) <= set(expected):
        raise ValueError("Empty or unexpected conformer set")
    members = [key for group in groups for key in group]
    if any(not group for group in groups) or sorted(members) != sorted(expected):
        raise ValueError("Duplicate groups must partition the declared set exactly")
    missing = [key for key in expected if key not in lookup or lookup[key].get("G_eV") is None]
    if missing:
        return dict(status="pending_native_G", missing_keys=missing, G_eV=None)
    for row in records:
        finite_number(row["G_eV"], "conformer G")
        finite_number(row["E_eV"], "electronic energy")
    representatives = [
        min(group, key=lambda key: (lookup[key]["G_eV"], key)) for group in groups
    ]
    result = partition([lookup[key]["G_eV"] for key in representatives], temperature_K)
    result.update(
        status="finite_set_partition_computed",
        representative_keys=representatives,
        duplicate_groups=groups,
        minimum_E_eV=min(row["E_eV"] for row in records),
        search_hits_are_weights=False,
        scope="Conditional finite conformer set; not a globally complete ensemble",
    )
    return result


def batch_change(first, expanded, tolerance_kJ_mol=2):
    """Report a within-protocol convergence diagnostic, never a global proof."""
    tolerance = finite_number(tolerance_kJ_mol, "convergence tolerance")
    if tolerance <= 0:
        raise ValueError("Positive convergence tolerance required")
    if any(row.get("status") != "finite_set_partition_computed" for row in (first, expanded)):
        return dict(status="pending_complete_batches", within_batch_tolerance=None)
    if not set(sum(first["duplicate_groups"], [])) < set(sum(expanded["duplicate_groups"], [])):
        raise ValueError("Convergence comparison requires strictly nested input sets")
    changes = {
        "minimum_E_change_kJ_mol": (expanded["minimum_E_eV"] - first["minimum_E_eV"]) * EV_KJ_MOL,
        "minimum_G_change_kJ_mol": (expanded["minimum_G_eV"] - first["minimum_G_eV"]) * EV_KJ_MOL,
        "ensemble_G_change_kJ_mol": (expanded["G_eV"] - first["G_eV"]) * EV_KJ_MOL,
    }
    return dict(
        status="within_protocol_batch_comparison",
        **changes,
        within_batch_tolerance=all(abs(value) <= tolerance for value in changes.values()),
        tolerance_kJ_mol=tolerance,
        global_ensemble_convergence_proven=False,
    )


def neighbor_cycle(empty, first, second, both):
    """Balanced conditional nonadditivity; no pure motif-binding interpretation."""
    values = [finite_number(value, "ensemble G") for value in (empty, first, second, both)]
    return math.fsum((values[3], values[0], -values[1], -values[2])) * EV_KJ_MOL
