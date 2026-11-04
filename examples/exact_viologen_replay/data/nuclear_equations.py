"""Frozen independent CPU replay function; source recorded in manifest."""
import math
import numpy as np

def nuclear_equations(frequencies, coordinates, mass, symmetry, cutoff):
    """Independent Python equations from the pinned upstream thermodyn source."""
    temperature, gas_constant, boltzmann = 298.15, 1.98726, 3.166808578545117e-6
    beta = 1 / (boltzmann * temperature)
    omega = np.abs(frequencies) / 219474.63067
    centered = coordinates - np.average(coordinates, axis=0, weights=mass)
    inertia = sum(m * (np.eye(3) * np.dot(x, x) - np.outer(x, x)) for m, x in zip(mass, centered))
    inertia += np.array([[1., 2., 4.], [2., 3., 5.], [4., 5., 6.]]) * 1e-10
    moments = np.linalg.eigvalsh(inertia)
    rotational = np.array([16.8576522 / v if v >= 3e-4 else 0. for v in moments])[::-1]
    average_SI = 1e-47 * 1.66053 * np.mean(moments)
    average_au = average_SI / 9.10938356e-31 * (1 / .52917726)**2 * 1e20
    mu = .5 / (omega + 1e-14)
    reduced = mu * average_au / (mu + average_au)
    harmonic_weight = 1 / (1 + (cutoff / np.abs(frequencies))**4)
    x = omega * beta
    sv_ho = x / np.expm1(x) - np.log(-np.expm1(-x))
    sv_fr = .5 + .5 * np.log(2 * np.pi * reduced / beta)
    entropy_vib = gas_constant * np.sum(harmonic_weight * sv_ho + (1 - harmonic_weight) * sv_fr)
    enthalpy_vib = np.sum(omega / np.expm1(x)) * gas_constant / boltzmann
    if len(mass) == 1:
        entropy_rot, enthalpy_rot = 0., 0.
    elif rotational[-1] < 1e-10:
        entropy_rot = gas_constant * (math.log(1 / (beta * rotational[0] / 219474.63067 * symmetry)) + 1)
        enthalpy_rot = gas_constant * temperature
    else:
        entropy_rot = gas_constant * (.5 * (math.log(np.pi / np.prod(rotational / 219474.63067) / beta**3)
                                                    - 2 * math.log(symmetry)) + 1.5)
        enthalpy_rot = 1.5 * gas_constant * temperature
    entropy_trans = 2.2868 * (5 * math.log10(temperature) + 3 * math.log10(sum(mass))) - 2.3135
    enthalpy_trans = 2.5 * gas_constant * temperature
    return (.5 * sum(omega) + ((enthalpy_vib + enthalpy_rot + enthalpy_trans)
            - temperature * (entropy_vib + entropy_rot + entropy_trans)) / 1000 / 627.50947428) * 27.211386245988
