"""Matched arithmetic with explicit physical context and retained provenance.

This checks the declarations supplied by a reader. It does not replace that
reader's verification of raw energies, Hessians, convergence or source hashes.
"""
from copy import deepcopy
import math

CONTEXT_FIELDS = ('quantity','unit','electronic_protocol','thermal_definition',
                  'solvent','solvent_model','temperature_K','standard_state','reference')


def combine_properties(records, coefficients):
    if not records or set(records) != set(coefficients):
        raise ValueError('Exactly one coefficient per component required')
    context = None
    retained = {}
    terms = []
    for role, record in records.items():
        declared = record.get('context', {})
        if any(key not in declared or declared[key] is None for key in CONTEXT_FIELDS):
            raise ValueError('Complete method, solvent, temperature, standard state and reference required')
        if type(declared['temperature_K']) not in (int,float) or not math.isfinite(declared['temperature_K']) or declared['temperature_K'] <= 0:
            raise ValueError('A finite positive temperature is required')
        if context is None: context = declared
        elif context != declared: raise ValueError('Property contexts differ; declare a conversion before combination')
        state = record.get('state', {})
        if (type(state.get('charge')) is not int or type(state.get('multiplicity')) is not int
                or state['multiplicity'] < 1 or not state.get('identity')):
            raise ValueError('Exact component identity, charge and spin required')
        if record.get('qc') != 'independently_checked' or not record.get('sources'):
            raise ValueError('Checked component QC and source bindings required')
        value, weight = record.get('value'), coefficients[role]
        if any(type(v) not in (int,float) or not math.isfinite(v) for v in (value,weight)):
            raise ValueError('Finite values and coefficients required')
        terms.append(value*weight)
        retained[role] = dict(coefficient=weight, state=deepcopy(state),
                              qc=record['qc'], sources=deepcopy(record['sources']))
    return dict(value=math.fsum(terms), context=deepcopy(context),
                components=retained, context_verified=True, empirical_accuracy_validated=False)
