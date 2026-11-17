"""Provenance-addressed solution thermal diagnostics; never hide unstable modes.

The geometry/energy protocol and thermal model are separate. A stable GFN2-ALPB
minimum is evidence for this low-level thermal model, not a DFT Hessian validation.
Production adoption additionally requires a paired reference and transfer evidence.
"""
from __future__ import annotations
import hashlib
import json
import math
import re
from pathlib import Path

SOLUTION_THERMAL = dict(schema=1,model='GFN2-xTB',solvation='ALPB',solvent='acetonitrile',
    temperature_K=298.15,rotor_cutoff_cm=50,optimization='verytight',electronic_accuracy=0.1,
    imaginary_acceptance_cm=-5.0,geometry='stationary model geometry; DFT source preserved',
    branch_search='both signs of lowest unstable mode, maximum two follow-up levels')


def thermal_identity(geometry: Path,charge:int,mult:int,protocol=None):
    return dict(geometry_sha256=hashlib.sha256(Path(geometry).read_bytes()).hexdigest(),
        charge=int(charge),mult=int(mult),protocol=dict(SOLUTION_THERMAL if protocol is None else protocol))


def thermal_hash(identity):
    return hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:16]


def spectrum_values(text):
    values=[]
    for line in text.splitlines():
        tokens=line.split()
        if not tokens or not tokens[0].isdigit():continue
        for token in tokens[1:]:
            try:
                values.append(float(token));break
            except ValueError:pass
    if not values or not all(math.isfinite(x) for x in values):
        raise ValueError('finite vibrational spectrum required')
    return values


def stationary_verdict(record,expected_identity=None):
    """Reject deep modes even when xTB's RRHO imaginary count excluded them."""
    reasons=[]
    if expected_identity is not None and record.get('identity')!=expected_identity:reasons.append('source_or_protocol_mismatch')
    if record.get('returncode')!=0:reasons.append('failed_execution')
    if record.get('optimization_converged') is not True:reasons.append('nonstationary')
    if record.get('connectivity_intact') is not True:reasons.append('changed_chemistry')
    frequencies=record.get('frequencies_cm')
    if not isinstance(frequencies,list) or not frequencies or not all(isinstance(x,(int,float)) and math.isfinite(x) for x in frequencies):
        reasons.append('missing_or_invalid_spectrum')
    elif min(frequencies)<SOLUTION_THERMAL['imaginary_acceptance_cm']:reasons.append('unstable_mode')
    if record.get('n_imag')!=0:reasons.append('imaginary_modes')
    g=record.get('g_thermal_eV')
    if not isinstance(g,(int,float)) or not math.isfinite(g):reasons.append('missing_thermal')
    return not reasons,reasons


def paired_potential(e_ox,e_red,g_ox,g_red,e_fc_ox,e_fc_neu,g_fc_ox,g_fc_neu):
    """One-electron potential with identical thermal model on all four species."""
    values=(e_ox,e_red,g_ox,g_red,e_fc_ox,e_fc_neu,g_fc_ox,g_fc_neu)
    if any(not isinstance(x,(int,float)) or not math.isfinite(x) for x in values):return None
    return (e_ox+g_ox)-(e_red+g_red)-((e_fc_ox+g_fc_ox)-(e_fc_neu+g_fc_neu))


def spin_degeneracy_eV(mult,temperature_K=298.15):
    """Ground-multiplet contribution absent from xTB's nuclear qRRHO routine.

    Uses an unsplit spin multiplet (no spin-orbit/Zeeman splitting or additional
    orbital degeneracy). Keep this assumption explicit in any thermal protocol.
    """
    if not isinstance(mult,int) or mult<1 or not math.isfinite(temperature_K) or temperature_K<=0:
        raise ValueError('positive integer multiplicity and temperature required')
    return -8.617333262145e-5*temperature_K*math.log(mult)


def load_solution_thermal_record(path,geometry,charge,mult,require_transfer_validation=False):
    """Fail closed on stale spin/geometry, failed branches and unsupported transfer claims.

    A returned model record may be used for explicitly labeled model diagnostics.
    Primary thermochemical predictions require paired transfer evidence as well.
    """
    try:
        record=json.loads(Path(path).read_text())
        identity=thermal_identity(geometry,charge,mult)
        selected=record.get('selected')
        if not isinstance(selected,dict) or record.get('identity')!=identity:return None
        valid,_=stationary_verdict(selected,identity)
        if not valid:return None
        n_atoms=int(Path(geometry).read_text().splitlines()[0])
        if len(selected['frequencies_cm'])!=3*n_atoms:return None
        by_cutoff=record.get('g_thermal_by_cutoff_eV') or {}
        term=by_cutoff.get('50')
        if not isinstance(term,(int,float)) or not math.isfinite(term):return None
        if abs(term-selected['g_thermal_eV'])>1e-5:return None
        if require_transfer_validation:
            evidence=record.get('transfer_evidence') or {}
            if (evidence.get('status')!='paired_transfer_validated'
                    or evidence.get('thermal_protocol')!=identity['protocol']
                    or evidence.get('geometry_sha256')!=identity['geometry_sha256']
                    or evidence.get('charge')!=charge or evidence.get('mult')!=mult):return None
        return record
    except (OSError,ValueError,TypeError,KeyError):return None
