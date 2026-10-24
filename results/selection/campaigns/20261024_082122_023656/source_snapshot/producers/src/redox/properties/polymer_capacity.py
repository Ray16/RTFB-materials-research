"""Composition-defined capacity bounds for neutral Merrifield substitution.

The bound assumes a declared electron budget per graft and full utilization.
It cannot establish an accessible capacity, reaction yield or loading distribution.
The mass convention is an H-capped linear fragment; optional added inactive mass
must be explicit. Charged attachment routes require their own counterion model.
"""
import math

from rdkit import Chem
from rdkit.Chem import Descriptors

from redox.build.grafting import precursor

FARADAY_C_MOL = 1.602176634e-19 * 6.02214076e23


def integer(value, name, minimum):
    if type(value) is not int or value < minimum:
        raise ValueError(f"Integer {name} >= {minimum} required")
    return value


def substitution_masses(fragment):
    molecule, reaction = precursor(fragment)
    if reaction != "NH_OH_substitution" or Chem.GetFormalCharge(molecule) != 0:
        raise ValueError("This bound requires neutral N-H/O-H substitution; define an ionic mass convention separately")
    table = Chem.GetPeriodicTable()
    carbon, hydrogen, chloride = [table.GetAtomicWeight(number) for number in (6, 1, 17)]
    return dict(
        precursor_g_mol=Descriptors.MolWt(molecule),
        styrene_repeat_g_mol=8 * carbon + 8 * hydrogen,
        chloromethyl_increment_g_mol=carbon + hydrogen + chloride,
        substitution_increment_g_mol=Descriptors.MolWt(molecule) - hydrogen - chloride,
        two_H_end_caps_g_mol=2 * hydrogen,
    )


def fragment_capacity_bound(fragment, repeat_units, reactive_handles, grafts, electrons_per_graft,
                            additional_inactive_mass_g_mol=0.0):
    repeat_units = integer(repeat_units, "repeat count", 1)
    reactive_handles = integer(reactive_handles, "reactive-handle count", 0)
    grafts = integer(grafts, "graft count", 0)
    electrons_per_graft = integer(electrons_per_graft, "electron budget", 1)
    if not grafts <= reactive_handles <= repeat_units:
        raise ValueError("Grafts must occupy declared reactive handles")
    if type(additional_inactive_mass_g_mol) not in (float, int) or not math.isfinite(additional_inactive_mass_g_mol) or additional_inactive_mass_g_mol < 0:
        raise ValueError("Additional inactive mass must be finite and nonnegative")
    mass = substitution_masses(fragment)
    total = (repeat_units * mass["styrene_repeat_g_mol"] + reactive_handles * mass["chloromethyl_increment_g_mol"]
        + grafts * mass["substitution_increment_g_mol"] + mass["two_H_end_caps_g_mol"] + additional_inactive_mass_g_mol)
    if not math.isfinite(total) or total <= 0:
        raise ValueError("Positive finite fragment mass required")
    capacity = grafts * electrons_per_graft * FARADAY_C_MOL / (3.6 * total)
    return dict(
        mass_g_mol=total, capacity_upper_bound_mAh_g=capacity,
        capacity_lower_bound_mAh_g=0.0, repeat_units=repeat_units,
        reactive_handles=reactive_handles, grafts=grafts,
        reactive_handle_density=reactive_handles / repeat_units,
        grafts_per_repeat=grafts / repeat_units,
        reactive_handle_conversion=grafts / reactive_handles if reactive_handles else None,
        electrons_per_graft=electrons_per_graft,
        additional_inactive_mass_g_mol=float(additional_inactive_mass_g_mol),
        assumptions="Full utilization of the declared electron budget; neutral H-capped linear fragment; no unlisted crosslinker/electrolyte mass",
        experimental_loading_inferred=False, accessible_capacity_established=False,
    )


def target_exceeds_upper_bound(result, target_mAh_g):
    if type(target_mAh_g) not in (float, int) or not math.isfinite(target_mAh_g) or target_mAh_g < 0:
        raise ValueError("A finite nonnegative capacity target required")
    # False means only that this necessary condition does not reject the
    # scenario. It never means all-axis admission or achieved capacity.
    return target_mAh_g > result["capacity_upper_bound_mAh_g"]
